#![forbid(unsafe_code)]
#![allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::must_use_candidate
)]

mod detect;
mod output_writer;
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

pub use output_writer::{Compression, OutputWriter};
pub use reader::Reader;
pub use record::{OwnedRecord, Record};
pub use writer::Writer;

const INPUT_BUFFER: usize = 256 * 1024;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Fasta,
    Fastq,
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

#[cfg(test)]
mod tests {
    use std::io::{Cursor, Read, Write};

    use super::*;

    fn gzip_member(bytes: &[u8]) -> Vec<u8> {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(bytes).unwrap();
        encoder.finish().unwrap()
    }

    fn bgzf_block(bytes: &[u8]) -> Vec<u8> {
        let mut encoder =
            flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(bytes).unwrap();
        let deflated = encoder.finish().unwrap();

        let block_size = 18 + deflated.len() + 8;
        let bsize = u16::try_from(block_size - 1).unwrap();
        let mut block = vec![
            0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0,
        ];
        block.extend_from_slice(&bsize.to_le_bytes());
        block.extend_from_slice(&deflated);

        let mut crc = flate2::Crc::new();
        crc.update(bytes);
        block.extend_from_slice(&crc.sum().to_le_bytes());
        block.extend_from_slice(&(bytes.len() as u32).to_le_bytes());
        assert_eq!(block.len(), block_size);
        block
    }

    fn bgzf_stream(bytes: &[u8]) -> Vec<u8> {
        let mut encoded = Vec::new();
        for chunk in bytes.chunks(32 * 1024) {
            encoded.extend_from_slice(&bgzf_block(chunk));
        }
        encoded.extend_from_slice(&bgzf_block(&[]));
        encoded
    }

    fn generic_reader_fails(bytes: &[u8]) -> bool {
        let Ok(mut reader) = open_reader(Cursor::new(bytes)) else {
            return true;
        };
        loop {
            match reader.read_record() {
                Ok(Some(_)) => {}
                Ok(None) => return false,
                Err(_) => return true,
            }
        }
    }

    fn path_reader_fails(bytes: &[u8]) -> bool {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(bytes).unwrap();
        file.flush().unwrap();
        let Ok(mut reader) = open_path(file.path()) else {
            return true;
        };
        loop {
            match reader.read_record() {
                Ok(Some(_)) => {}
                Ok(None) => return false,
                Err(_) => return true,
            }
        }
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
    fn real_multiblock_bgzf_decodes_through_generic_and_path_readers() {
        let mut raw = Vec::new();
        for index in 0..2_000 {
            writeln!(raw, "@read{index}\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII").unwrap();
        }
        let encoded = bgzf_stream(&raw);
        assert_eq!(
            &encoded[..16],
            &[
                0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0
            ]
        );
        let first_block_size = usize::from(u16::from_le_bytes([encoded[16], encoded[17]])) + 1;
        assert!(first_block_size < encoded.len());

        let mut generic = open_reader(Cursor::new(&encoded)).unwrap();
        let mut generic_count = 0;
        while generic.read_record().unwrap().is_some() {
            generic_count += 1;
        }
        assert_eq!(generic_count, 2_000);

        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(&encoded).unwrap();
        file.flush().unwrap();
        let mut path = open_path(file.path()).unwrap();
        let mut path_count = 0;
        while path.read_record().unwrap().is_some() {
            path_count += 1;
        }
        assert_eq!(path_count, 2_000);
    }

    #[test]
    fn bgzf_crc_corruption_fails_generic_and_path_readers() {
        let mut encoded = bgzf_stream(b"@one\nACGT\n+\nIIII\n");
        let first_block_size = usize::from(u16::from_le_bytes([encoded[16], encoded[17]])) + 1;
        encoded[first_block_size - 8] ^= 0xff;

        assert!(generic_reader_fails(&encoded));
        assert!(path_reader_fails(&encoded));
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
    fn every_missing_gzip_trailer_byte_fails_generic_and_path_readers() {
        let encoded = gzip_member(b"@one\nACGT\n+\nIIII\n");
        for missing in 1..=8 {
            let truncated = &encoded[..encoded.len() - missing];
            assert!(
                generic_reader_fails(truncated),
                "generic reader accepted gzip missing {missing} trailer byte(s)"
            );
            assert!(
                path_reader_fails(truncated),
                "path reader accepted gzip missing {missing} trailer byte(s)"
            );
        }
    }
}
