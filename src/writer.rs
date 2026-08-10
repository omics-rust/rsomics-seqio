use std::io::{BufWriter, Write};

use rsomics_common::{Result, RsomicsError};

use crate::record::{are_valid_printable_bytes, is_valid_header_byte};
use crate::{Format, OwnedRecord, Record};

const BUFFER_CAPACITY: usize = 256 * 1024;

pub struct Writer<W: Write> {
    inner: BufWriter<W>,
    format: Format,
}

impl<W: Write> Writer<W> {
    #[must_use]
    pub fn new(inner: W, format: Format) -> Self {
        Self {
            inner: BufWriter::with_capacity(BUFFER_CAPACITY, inner),
            format,
        }
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
        self.inner.flush().map_err(RsomicsError::Io)
    }

    pub fn finish(self) -> Result<()> {
        self.finish_into_inner().map(drop)
    }

    pub fn finish_into_inner(self) -> Result<W> {
        let mut inner = self.inner;
        inner.flush().map_err(RsomicsError::Io)?;
        inner
            .into_inner()
            .map_err(|error| RsomicsError::Io(error.into_error()))
    }

    fn write_all(&mut self, bytes: &[u8]) -> Result<()> {
        self.inner.write_all(bytes).map_err(RsomicsError::Io)
    }
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
        (Format::Fasta, None) => {}
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
    use std::io::Cursor;

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
    fn empty_fasta_sequence_round_trips() {
        let mut output = Vec::new();
        let mut writer = Writer::new(&mut output, Format::Fasta);
        writer
            .write_record(Record {
                id: b"empty",
                seq: b"",
                qual: None,
            })
            .unwrap();
        writer.finish().unwrap();

        let mut reader = Reader::new(Cursor::new(&output), Format::Fasta);
        assert!(reader.read_record().unwrap().unwrap().seq.is_empty());
        assert!(reader.read_record().unwrap().is_none());
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
}
