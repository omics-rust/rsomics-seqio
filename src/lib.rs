#![forbid(unsafe_code)]
#![allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::must_use_candidate
)]

mod detect;
mod reader;
mod reader_gz;
mod record;
mod writer;

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use rsomics_common::{Result, RsomicsError};

use detect::{CompressionKind, ReplayReader};
use reader_gz::GzipStream;

pub use reader::Reader;
pub use record::{LegacyFastqRecord, OwnedRecord, Record};
pub use writer::{Writer, create_path};

const INPUT_BUFFER: usize = 256 * 1024;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Fasta,
    Fastq,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Compression {
    Plain,
    Gzip { level: u32 },
}

enum GenericInput<R: Read> {
    Plain(ReplayReader<R>),
    Gzip(Box<MultiGzDecoder<ReplayReader<R>>>),
}

impl<R: Read> Read for GenericInput<R> {
    fn read(&mut self, output: &mut [u8]) -> std::io::Result<usize> {
        match self {
            Self::Plain(reader) => reader.read(output),
            Self::Gzip(reader) => reader.read(output),
        }
    }
}

fn decode_reader<R: Read>(source: R) -> Result<BufReader<GenericInput<R>>> {
    let (compression, replayed) = detect::probe(source)?;
    let input = match compression {
        CompressionKind::Plain => GenericInput::Plain(replayed),
        CompressionKind::Gzip => GenericInput::Gzip(Box::new(MultiGzDecoder::new(replayed))),
    };
    Ok(BufReader::with_capacity(INPUT_BUFFER, input))
}

pub fn open_reader<R: Read>(source: R) -> Result<Reader<impl BufRead>> {
    Reader::detect(decode_reader(source)?)
}

pub fn open_reader_with_format<R: Read>(source: R, format: Format) -> Result<Reader<impl BufRead>> {
    Ok(Reader::new(decode_reader(source)?, format))
}

enum PathInput {
    Plain(ReplayReader<File>),
    Gzip(GzipStream),
}

impl Read for PathInput {
    fn read(&mut self, output: &mut [u8]) -> std::io::Result<usize> {
        match self {
            Self::Plain(reader) => reader.read(output),
            Self::Gzip(reader) => reader.read(output),
        }
    }
}

pub struct PathReader {
    inner: Reader<BufReader<PathInput>>,
}

impl PathReader {
    #[must_use]
    pub fn format(&self) -> Format {
        self.inner.format()
    }

    pub fn read_record(&mut self) -> Result<Option<Record<'_>>> {
        self.inner.read_record()
    }
}

pub fn open_path(path: &Path) -> Result<PathReader> {
    let file = File::open(path).map_err(|error| {
        RsomicsError::Io(std::io::Error::new(
            error.kind(),
            format!("opening {}: {error}", path.display()),
        ))
    })?;
    let (compression, replayed) = detect::probe(file).map_err(|error| match error {
        RsomicsError::Io(error) => RsomicsError::Io(std::io::Error::new(
            error.kind(),
            format!("reading header of {}: {error}", path.display()),
        )),
        other => other,
    })?;
    let input = match compression {
        CompressionKind::Plain => PathInput::Plain(replayed),
        CompressionKind::Gzip => PathInput::Gzip(GzipStream::new(replayed)?),
    };
    let inner = Reader::detect(BufReader::with_capacity(INPUT_BUFFER, input))?;
    Ok(PathReader { inner })
}

pub struct LegacyFastqSource {
    inner: PathReader,
    done: bool,
}

impl Iterator for LegacyFastqSource {
    type Item = Result<LegacyFastqRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        match self.inner.read_record() {
            Ok(Some(record)) => Some(Ok(LegacyFastqRecord {
                id: record.id.to_vec(),
                seq: record.seq.to_vec(),
                qual: record
                    .qual
                    .expect("LegacyFastqSource only wraps FASTQ readers")
                    .to_vec(),
            })),
            Ok(None) => {
                self.done = true;
                None
            }
            Err(error) => {
                self.done = true;
                Some(Err(error))
            }
        }
    }
}

pub fn open_fastq_legacy(path: &Path) -> Result<LegacyFastqSource> {
    let inner = open_path(path)?;
    if inner.format() != Format::Fastq {
        return Err(RsomicsError::InvalidInput(format!(
            "{} contains FASTA, not FASTQ",
            path.display()
        )));
    }
    Ok(LegacyFastqSource { inner, done: false })
}

#[cfg(test)]
mod tests {
    use std::io::{Cursor, Read, Write};

    use super::*;

    fn gzip_member(bytes: &[u8]) -> Vec<u8> {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(bytes).unwrap();
        encoder.finish().unwrap()
    }

    struct OneByte<R>(R);

    impl<R: Read> Read for OneByte<R> {
        fn read(&mut self, output: &mut [u8]) -> std::io::Result<usize> {
            let limit = output.len().min(1);
            self.0.read(&mut output[..limit])
        }
    }

    #[test]
    fn generic_reader_detects_multimember_gzip() {
        let mut encoded = gzip_member(b">one\nAC");
        encoded.extend_from_slice(&gzip_member(b"GT\n"));
        let mut reader = open_reader(Cursor::new(encoded)).unwrap();
        let record = reader.read_record().unwrap().unwrap();
        assert_eq!(record.id, b"one");
        assert_eq!(record.seq, b"ACGT");
    }

    #[test]
    fn generic_reader_handles_one_byte_plain_and_gzip_sources() {
        let mut plain = open_reader(OneByte(Cursor::new(b">one\nACGT\n"))).unwrap();
        assert_eq!(plain.read_record().unwrap().unwrap().seq, b"ACGT");

        let encoded = gzip_member(b"@one\nACGT\n+\nIIII\n");
        let mut gzip = open_reader(OneByte(Cursor::new(encoded))).unwrap();
        assert_eq!(
            gzip.read_record().unwrap().unwrap().qual,
            Some(b"IIII".as_slice())
        );
    }

    #[test]
    fn explicit_format_supports_empty_and_compressed_readers() {
        let mut empty =
            open_reader_with_format(Cursor::new(Vec::<u8>::new()), Format::Fastq).unwrap();
        assert!(empty.read_record().unwrap().is_none());

        let encoded = gzip_member(b">one\nACGT\n");
        let mut gzip = open_reader_with_format(Cursor::new(encoded), Format::Fasta).unwrap();
        assert_eq!(gzip.read_record().unwrap().unwrap().seq, b"ACGT");
    }

    #[test]
    fn path_reader_detects_content_not_extension() {
        let mut file = tempfile::Builder::new()
            .suffix(".fastq")
            .tempfile()
            .unwrap();
        let encoded = gzip_member(b">one\nACGT\n");
        file.write_all(&encoded).unwrap();
        file.flush().unwrap();

        let mut reader = open_path(file.path()).unwrap();
        assert_eq!(reader.format(), Format::Fasta);
        assert_eq!(reader.read_record().unwrap().unwrap().seq, b"ACGT");
    }

    #[test]
    fn corrupt_gzip_is_reported_by_generic_and_path_openers() {
        let corrupt = [0x1f, 0x8b, 0x08, 0x00, 0xff, 0xff];
        assert!(open_reader(Cursor::new(corrupt)).is_err());

        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(&corrupt).unwrap();
        file.flush().unwrap();
        assert!(open_path(file.path()).is_err());
    }

    #[test]
    fn truncated_gzip_is_reported_after_any_decoded_records() {
        let mut encoded = gzip_member(b"@one\nACGT\n+\nIIII\n");
        encoded.truncate(encoded.len() - 6);

        let generic_result = open_reader(Cursor::new(encoded.clone()));
        match generic_result {
            Err(_) => {}
            Ok(mut reader) => {
                assert!(reader.read_record().unwrap().is_some());
                assert!(reader.read_record().is_err());
            }
        }

        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(&encoded).unwrap();
        file.flush().unwrap();
        let path_result = open_path(file.path());
        match path_result {
            Err(_) => {}
            Ok(mut reader) => {
                assert!(reader.read_record().unwrap().is_some());
                assert!(reader.read_record().is_err());
            }
        }
    }

    #[test]
    fn legacy_fastq_adapter_has_an_explicit_name_and_shape() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(b"@one\nACGT\n+\nIIII\n").unwrap();
        file.flush().unwrap();

        let records = open_fastq_legacy(file.path())
            .unwrap()
            .collect::<Result<Vec<_>>>()
            .unwrap();
        assert_eq!(
            records,
            vec![LegacyFastqRecord {
                id: b"one".to_vec(),
                seq: b"ACGT".to_vec(),
                qual: b"IIII".to_vec(),
            }]
        );
    }

    #[test]
    fn legacy_fastq_adapter_rejects_fasta() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(b">one\nACGT\n").unwrap();
        file.flush().unwrap();
        assert!(matches!(
            open_fastq_legacy(file.path()),
            Err(RsomicsError::InvalidInput(_))
        ));
    }
}
