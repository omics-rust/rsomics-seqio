use std::fs::File;
use std::ops::Range;
use std::path::{Path, PathBuf};

use noodles_core::{Position, Region};
use noodles_fasta as fasta;
use rsomics_common::{Context, Result, RsomicsError};

const DEFAULT_CACHE_CAPACITY: usize = 1024 * 1024;

pub struct IndexedFasta {
    reader: fasta::io::IndexedReader<fasta::io::BufReader<File>>,
    path: PathBuf,
    references: Box<[(Box<[u8]>, usize)]>,
    cache_capacity: usize,
    cached_name: Box<[u8]>,
    cached_start: usize,
    cached_sequence: Vec<u8>,
}

impl IndexedFasta {
    pub fn open(path: &Path) -> Result<Self> {
        Self::with_cache_capacity(path, DEFAULT_CACHE_CAPACITY)
    }

    pub fn with_cache_capacity(path: &Path, cache_capacity: usize) -> Result<Self> {
        if cache_capacity == 0 {
            return Err(RsomicsError::ConfigError(
                "indexed FASTA cache capacity must be greater than zero".to_owned(),
            ));
        }
        let reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .rs_with_context(|| format!("opening indexed FASTA {}", path.display()))?;
        let mut references = Vec::new();
        for record in reader.index().as_ref() {
            if references
                .iter()
                .any(|(name, _): &(Box<[u8]>, usize)| name.as_ref() == record.name())
            {
                return Err(RsomicsError::InvalidInput(format!(
                    "indexed FASTA {} has duplicate reference {}",
                    path.display(),
                    String::from_utf8_lossy(record.name())
                )));
            }
            let length = usize::try_from(record.length()).map_err(|_| {
                RsomicsError::InvalidInput(format!(
                    "indexed FASTA {} reference {} exceeds addressable length",
                    path.display(),
                    String::from_utf8_lossy(record.name())
                ))
            })?;
            references.push((record.name().to_vec().into_boxed_slice(), length));
        }
        Ok(Self {
            reader,
            path: path.to_path_buf(),
            references: references.into_boxed_slice(),
            cache_capacity,
            cached_name: Box::default(),
            cached_start: 0,
            cached_sequence: Vec::new(),
        })
    }

    pub fn len(&self, name: &[u8]) -> Result<usize> {
        self.references
            .iter()
            .find_map(|(candidate, length)| (candidate.as_ref() == name).then_some(*length))
            .ok_or_else(|| {
                RsomicsError::InvalidInput(format!(
                    "indexed FASTA {} has no reference named {}",
                    self.path.display(),
                    String::from_utf8_lossy(name)
                ))
            })
    }

    pub fn is_empty(&self, name: &[u8]) -> Result<bool> {
        self.len(name).map(|length| length == 0)
    }

    pub fn fetch(&mut self, name: &[u8], range: Range<usize>) -> Result<&[u8]> {
        let length = self.len(name)?;
        if range.start > range.end || range.end > length {
            return Err(RsomicsError::InvalidInput(format!(
                "indexed FASTA {} reference {} range {}..{} is outside length {length}",
                self.path.display(),
                String::from_utf8_lossy(name),
                range.start,
                range.end
            )));
        }
        if range.is_empty() {
            return Ok(&self.cached_sequence[..0]);
        }
        if self.cached_name.as_ref() != name
            || range.start < self.cached_start
            || range.end > self.cached_start + self.cached_sequence.len()
        {
            self.load(name, range.clone(), length)?;
        }
        let start = range.start - self.cached_start;
        let end = range.end - self.cached_start;
        Ok(&self.cached_sequence[start..end])
    }

    fn load(&mut self, name: &[u8], range: Range<usize>, length: usize) -> Result<()> {
        let start = range.start / self.cache_capacity * self.cache_capacity;
        let end = start
            .checked_add(self.cache_capacity)
            .map_or(length, |end| end.min(length))
            .max(range.end);
        let interval_start = Position::try_from(start + 1).map_err(|_| {
            RsomicsError::InvalidInput(format!(
                "indexed FASTA {} reference position exceeds supported range",
                self.path.display()
            ))
        })?;
        let interval_end = Position::try_from(end).map_err(|_| {
            RsomicsError::InvalidInput(format!(
                "indexed FASTA {} reference position exceeds supported range",
                self.path.display()
            ))
        })?;
        let record = self
            .reader
            .query(&Region::new(name.to_vec(), interval_start..=interval_end))
            .rs_with_context(|| {
                format!(
                    "reading indexed FASTA {} reference {} range {start}..{end}",
                    self.path.display(),
                    String::from_utf8_lossy(name)
                )
            })?;
        self.cached_sequence.clear();
        self.cached_sequence
            .extend_from_slice(record.sequence().as_ref());
        self.cached_name = name.into();
        self.cached_start = start;
        Ok(())
    }
}
