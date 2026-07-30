use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use flate2::Compression as FlateCompression;
use flate2::write::GzEncoder;
use rsomics_common::{Result, RsomicsError};

use crate::record::{are_valid_printable_bytes, is_valid_header_byte};
use crate::{Compression, Format, OwnedRecord, Record};

const BUFFER_CAPACITY: usize = 256 * 1024;

enum Backend<W: Write> {
    Plain(BufWriter<W>),
    Gzip(Box<GzEncoder<BufWriter<W>>>),
}

pub struct Writer<W: Write> {
    backend: Backend<W>,
    format: Format,
}

impl<W: Write> Writer<W> {
    #[must_use]
    pub fn new(inner: W, format: Format) -> Self {
        Self {
            backend: Backend::Plain(BufWriter::with_capacity(BUFFER_CAPACITY, inner)),
            format,
        }
    }

    pub fn gzip(inner: W, format: Format, level: u32) -> Result<Self> {
        validate_gzip_level(level)?;
        let buffered = BufWriter::with_capacity(BUFFER_CAPACITY, inner);
        Ok(Self {
            backend: Backend::Gzip(Box::new(GzEncoder::new(
                buffered,
                FlateCompression::new(level),
            ))),
            format,
        })
    }

    #[must_use]
    pub fn format(&self) -> Format {
        self.format
    }

    pub fn write_record(&mut self, record: Record<'_>) -> Result<()> {
        validate_record(record, self.format)?;
        match self.format {
            Format::Fasta => {
                self.write_all(b">")?;
                self.write_all(record.id)?;
                self.write_all(b"\n")?;
                self.write_all(record.seq)?;
                self.write_all(b"\n")
            }
            Format::Fastq => {
                self.write_all(b"@")?;
                self.write_all(record.id)?;
                self.write_all(b"\n")?;
                self.write_all(record.seq)?;
                self.write_all(b"\n+\n")?;
                self.write_all(record.qual.expect("validated FASTQ quality"))?;
                self.write_all(b"\n")
            }
        }
    }

    pub fn write_owned(&mut self, record: &OwnedRecord) -> Result<()> {
        self.write_record(record.as_record())
    }

    pub fn flush(&mut self) -> Result<()> {
        match &mut self.backend {
            Backend::Plain(writer) => writer.flush().map_err(RsomicsError::Io),
            Backend::Gzip(writer) => writer.flush().map_err(RsomicsError::Io),
        }
    }

    pub fn finish(self) -> Result<()> {
        self.finish_into_inner().map(drop)
    }

    pub fn finish_into_inner(self) -> Result<W> {
        let mut buffered = match self.backend {
            Backend::Plain(mut writer) => {
                writer.flush().map_err(RsomicsError::Io)?;
                writer
            }
            Backend::Gzip(writer) => writer.finish().map_err(RsomicsError::Io)?,
        };
        buffered.flush().map_err(RsomicsError::Io)?;
        buffered
            .into_inner()
            .map_err(|error| RsomicsError::Io(error.into_error()))
    }

    fn write_all(&mut self, bytes: &[u8]) -> Result<()> {
        match &mut self.backend {
            Backend::Plain(writer) => writer.write_all(bytes).map_err(RsomicsError::Io),
            Backend::Gzip(writer) => writer.write_all(bytes).map_err(RsomicsError::Io),
        }
    }
}

pub fn create_path(path: &Path, format: Format, compression: Compression) -> Result<Writer<File>> {
    if let Compression::Gzip { level } = compression {
        validate_gzip_level(level)?;
    }
    let file = File::create(path).map_err(|error| {
        RsomicsError::Io(std::io::Error::new(
            error.kind(),
            format!("creating {}: {error}", path.display()),
        ))
    })?;
    match compression {
        Compression::Plain => Ok(Writer::new(file, format)),
        Compression::Gzip { level } => Writer::gzip(file, format, level),
    }
}

fn validate_gzip_level(level: u32) -> Result<()> {
    if level > 9 {
        return Err(RsomicsError::InvalidInput(format!(
            "gzip level must be in 0..=9 (got {level})"
        )));
    }
    Ok(())
}

fn validate_record(record: Record<'_>, format: Format) -> Result<()> {
    if record.id.is_empty() {
        return Err(RsomicsError::InvalidInput(
            "sequence identifier must not be empty".into(),
        ));
    }
    if record
        .id
        .iter()
        .copied()
        .any(|byte| !is_valid_header_byte(byte))
    {
        return Err(RsomicsError::InvalidInput(
            "sequence identifier contains an invalid byte".into(),
        ));
    }
    if !are_valid_printable_bytes(record.seq) {
        return Err(RsomicsError::InvalidInput(
            "sequence contains a byte outside printable non-space ASCII 33..=126".into(),
        ));
    }

    match (format, record.qual) {
        (Format::Fasta, None) => {
            if record.seq.is_empty() {
                return Err(RsomicsError::InvalidInput(
                    "FASTA sequence must not be empty".into(),
                ));
            }
        }
        (Format::Fasta, Some(_)) => {
            return Err(RsomicsError::InvalidInput(
                "cannot write quality scores in FASTA format".into(),
            ));
        }
        (Format::Fastq, None) => {
            return Err(RsomicsError::InvalidInput(
                "FASTQ record requires quality scores".into(),
            ));
        }
        (Format::Fastq, Some(quality)) => {
            if record.seq.len() != quality.len() {
                return Err(RsomicsError::InvalidInput(format!(
                    "FASTQ sequence/quality length mismatch: {} vs {}",
                    record.seq.len(),
                    quality.len()
                )));
            }
            if !are_valid_printable_bytes(quality) {
                return Err(RsomicsError::InvalidInput(
                    "FASTQ quality contains a byte outside ASCII 33..=126".into(),
                ));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::io::{Cursor, Read, Write as _};

    use flate2::read::MultiGzDecoder;

    use super::*;
    use crate::Reader;

    #[test]
    fn generic_cursor_fasta_round_trip() {
        let mut output = Vec::new();
        {
            let mut writer = Writer::new(&mut output, Format::Fasta);
            writer
                .write_record(Record {
                    id: b"one",
                    seq: b"ACGT",
                    qual: None,
                })
                .unwrap();
            writer.finish().unwrap();
        }

        let mut reader = Reader::detect(Cursor::new(&output)).unwrap();
        let record = reader.read_record().unwrap().unwrap();
        assert_eq!(record.id, b"one");
        assert_eq!(record.seq, b"ACGT");
        assert!(record.qual.is_none());
    }

    #[test]
    fn generic_cursor_fastq_round_trip() {
        let mut output = Vec::new();
        {
            let mut writer = Writer::new(&mut output, Format::Fastq);
            writer
                .write_record(Record {
                    id: b"one",
                    seq: b"ACGT",
                    qual: Some(b"IIII"),
                })
                .unwrap();
            writer.finish().unwrap();
        }

        let mut reader = Reader::detect(Cursor::new(&output)).unwrap();
        let record = reader.read_record().unwrap().unwrap();
        assert_eq!(record.qual, Some(b"IIII".as_slice()));
    }

    #[test]
    fn owned_writer_can_return_its_sink() {
        let mut writer = Writer::new(Vec::new(), Format::Fasta);
        writer
            .write_record(Record {
                id: b"one",
                seq: b"ACGT",
                qual: None,
            })
            .unwrap();
        assert_eq!(
            writer.finish_into_inner().unwrap(),
            b">one\nACGT\n".to_vec()
        );
    }

    #[test]
    fn gzip_writer_round_trip() {
        let mut compressed = Vec::new();
        {
            let mut writer = Writer::gzip(&mut compressed, Format::Fastq, 4).unwrap();
            writer
                .write_record(Record {
                    id: b"one",
                    seq: b"ACGT",
                    qual: Some(b"IIII"),
                })
                .unwrap();
            writer.finish().unwrap();
        }

        let mut decoded = Vec::new();
        MultiGzDecoder::new(compressed.as_slice())
            .read_to_end(&mut decoded)
            .unwrap();
        assert_eq!(decoded, b"@one\nACGT\n+\nIIII\n");
    }

    #[test]
    fn gzip_path_round_trip_uses_content_detection() {
        let file = tempfile::Builder::new().suffix(".data").tempfile().unwrap();
        let mut writer =
            create_path(file.path(), Format::Fastq, Compression::Gzip { level: 4 }).unwrap();
        writer
            .write_record(Record {
                id: b"one",
                seq: b"ACGT",
                qual: Some(b"IIII"),
            })
            .unwrap();
        writer.finish().unwrap();

        let mut reader = crate::open_path(file.path()).unwrap();
        assert_eq!(reader.format(), Format::Fastq);
        assert_eq!(reader.read_record().unwrap().unwrap().id, b"one");
    }

    #[test]
    fn writer_rejects_format_mismatch() {
        let mut writer = Writer::new(Vec::new(), Format::Fastq);
        assert!(matches!(
            writer.write_record(Record {
                id: b"one",
                seq: b"ACGT",
                qual: None,
            }),
            Err(RsomicsError::InvalidInput(_))
        ));
    }

    #[test]
    fn invalid_gzip_level_does_not_truncate_path() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(b"keep").unwrap();
        file.flush().unwrap();
        assert!(matches!(
            create_path(file.path(), Format::Fastq, Compression::Gzip { level: 10 }),
            Err(RsomicsError::InvalidInput(_))
        ));
        assert_eq!(std::fs::read(file.path()).unwrap(), b"keep");
    }

    #[test]
    fn writer_accepts_biological_symbols_and_rejects_control_bytes() {
        let mut writer = Writer::new(Vec::new(), Format::Fasta);
        writer
            .write_record(Record {
                id: b"one",
                seq: b"ACGTN-.*?",
                qual: None,
            })
            .unwrap();

        for invalid in [b"A\0C".as_slice(), b"A\x7fC"] {
            let mut writer = Writer::new(Vec::new(), Format::Fasta);
            assert!(matches!(
                writer.write_record(Record {
                    id: b"one",
                    seq: invalid,
                    qual: None,
                }),
                Err(RsomicsError::InvalidInput(_))
            ));
        }
    }

    #[test]
    fn writer_allows_empty_fastq_sequence_and_quality() {
        let mut writer = Writer::new(Vec::new(), Format::Fastq);
        writer
            .write_record(Record {
                id: b"one",
                seq: b"",
                qual: Some(b""),
            })
            .unwrap();
        assert_eq!(
            writer.finish_into_inner().unwrap(),
            b"@one\n\n+\n\n".to_vec()
        );
    }

    struct FailingSink;

    impl Write for FailingSink {
        fn write(&mut self, _bytes: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::other("intentional sink failure"))
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Err(std::io::Error::other("intentional sink failure"))
        }
    }

    #[test]
    fn writer_propagates_sink_failures() {
        let writer = Writer::new(FailingSink, Format::Fasta);
        assert!(matches!(writer.finish(), Err(RsomicsError::Io(_))));

        let mut writer = Writer::new(FailingSink, Format::Fasta);
        let long_sequence = vec![b'A'; BUFFER_CAPACITY + 1];
        assert!(matches!(
            writer.write_record(Record {
                id: b"one",
                seq: &long_sequence,
                qual: None,
            }),
            Err(RsomicsError::Io(_))
        ));
    }

    #[test]
    fn gzip_writer_propagates_finish_sink_failure() {
        let mut writer = Writer::gzip(FailingSink, Format::Fastq, 4).unwrap();
        writer
            .write_record(Record {
                id: b"one",
                seq: b"ACGT",
                qual: Some(b"IIII"),
            })
            .unwrap();
        assert!(matches!(
            writer.finish_into_inner(),
            Err(RsomicsError::Io(_))
        ));
    }
}
