use std::io::BufRead;

use rsomics_common::{Result, RsomicsError};

use crate::record::{is_valid_header_byte, is_valid_quality_byte, is_valid_sequence_byte};
use crate::{Format, Record};

const BUFFER_CAPACITY: usize = 64 * 1024;

pub struct Reader<R> {
    inner: R,
    format: Format,
    line: Vec<u8>,
    id: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    pending_header: Vec<u8>,
    has_pending_header: bool,
    line_number: u64,
    failed: bool,
}

impl<R: BufRead> Reader<R> {
    #[must_use]
    pub fn new(inner: R, format: Format) -> Self {
        Self {
            inner,
            format,
            line: Vec::with_capacity(BUFFER_CAPACITY),
            id: Vec::with_capacity(256),
            seq: Vec::with_capacity(BUFFER_CAPACITY),
            qual: Vec::with_capacity(BUFFER_CAPACITY),
            pending_header: Vec::with_capacity(256),
            has_pending_header: false,
            line_number: 0,
            failed: false,
        }
    }

    pub fn detect(mut inner: R) -> Result<Self> {
        let first = inner.fill_buf().map_err(RsomicsError::Io)?;
        let format = match first.first() {
            Some(b'>') => Format::Fasta,
            Some(b'@') => Format::Fastq,
            Some(&byte) => {
                return Err(RsomicsError::InvalidInput(format!(
                    "expected FASTA '>' or FASTQ '@' at byte 1, got {byte:?}"
                )));
            }
            None => {
                return Err(RsomicsError::InvalidInput(
                    "cannot detect sequence format from empty input".into(),
                ));
            }
        };
        Ok(Self::new(inner, format))
    }

    #[must_use]
    pub fn format(&self) -> Format {
        self.format
    }

    pub fn read_record(&mut self) -> Result<Option<Record<'_>>> {
        if self.failed {
            return Err(RsomicsError::InvalidInput(
                "reader cannot continue after a previous record error".into(),
            ));
        }

        let parsed = match self.format {
            Format::Fasta => self.read_fasta(),
            Format::Fastq => self.read_fastq(),
        };
        match parsed {
            Ok(false) => Ok(None),
            Ok(true) => Ok(Some(Record {
                id: &self.id,
                seq: &self.seq,
                qual: (self.format == Format::Fastq).then_some(self.qual.as_slice()),
            })),
            Err(error) => {
                self.failed = true;
                Err(error)
            }
        }
    }

    #[must_use]
    pub fn into_inner(self) -> R {
        self.inner
    }

    fn read_fasta(&mut self) -> Result<bool> {
        self.id.clear();
        self.seq.clear();
        self.qual.clear();

        if self.has_pending_header {
            std::mem::swap(&mut self.line, &mut self.pending_header);
            self.has_pending_header = false;
        } else if !read_line(&mut self.inner, &mut self.line, &mut self.line_number)? {
            return Ok(false);
        }

        if self.line.first() != Some(&b'>') {
            return Err(invalid_at(
                self.line_number,
                "expected FASTA header starting with '>'",
            ));
        }
        if self.line.len() == 1 {
            return Err(invalid_at(self.line_number, "empty FASTA identifier"));
        }
        validate_header(&self.line[1..], self.line_number)?;
        self.id.extend_from_slice(&self.line[1..]);

        let mut sequence_lines = 0u64;
        loop {
            if !read_line(&mut self.inner, &mut self.line, &mut self.line_number)? {
                break;
            }
            if self.line.first() == Some(&b'>') {
                std::mem::swap(&mut self.line, &mut self.pending_header);
                self.has_pending_header = true;
                break;
            }
            if self.line.is_empty() {
                return Err(invalid_at(
                    self.line_number,
                    "empty line inside FASTA sequence",
                ));
            }
            validate_sequence(&self.line, self.line_number)?;
            self.seq.extend_from_slice(&self.line);
            sequence_lines += 1;
        }

        if sequence_lines == 0 {
            return Err(invalid_at(self.line_number, "FASTA record has no sequence"));
        }

        Ok(true)
    }

    fn read_fastq(&mut self) -> Result<bool> {
        self.id.clear();
        self.seq.clear();
        self.qual.clear();

        if !read_line(&mut self.inner, &mut self.line, &mut self.line_number)? {
            return Ok(false);
        }
        if self.line.first() != Some(&b'@') {
            return Err(invalid_at(
                self.line_number,
                "expected FASTQ header starting with '@'",
            ));
        }
        if self.line.len() == 1 {
            return Err(invalid_at(self.line_number, "empty FASTQ identifier"));
        }
        validate_header(&self.line[1..], self.line_number)?;
        self.id.extend_from_slice(&self.line[1..]);

        if !read_line(&mut self.inner, &mut self.line, &mut self.line_number)? {
            return Err(invalid_at(
                self.line_number,
                "truncated FASTQ: missing sequence line",
            ));
        }
        validate_sequence(&self.line, self.line_number)?;
        self.seq.extend_from_slice(&self.line);

        if !read_line(&mut self.inner, &mut self.line, &mut self.line_number)? {
            return Err(invalid_at(
                self.line_number,
                "truncated FASTQ: missing '+' separator",
            ));
        }
        if self.line.first() != Some(&b'+') {
            return Err(invalid_at(self.line_number, "expected FASTQ '+' separator"));
        }
        if self.line.len() > 1 && self.line[1..] != self.id {
            return Err(invalid_at(
                self.line_number,
                "FASTQ '+' identifier does not match the record identifier",
            ));
        }

        if !read_line(&mut self.inner, &mut self.line, &mut self.line_number)? {
            return Err(invalid_at(
                self.line_number,
                "truncated FASTQ: missing quality line",
            ));
        }
        validate_quality(&self.line, self.line_number)?;
        self.qual.extend_from_slice(&self.line);

        if self.seq.len() != self.qual.len() {
            return Err(invalid_at(
                self.line_number,
                &format!(
                    "FASTQ sequence/quality length mismatch: {} vs {}",
                    self.seq.len(),
                    self.qual.len()
                ),
            ));
        }

        Ok(true)
    }
}

fn read_line<R: BufRead>(
    inner: &mut R,
    buffer: &mut Vec<u8>,
    line_number: &mut u64,
) -> Result<bool> {
    buffer.clear();
    let read = inner.read_until(b'\n', buffer).map_err(RsomicsError::Io)?;
    if read == 0 {
        return Ok(false);
    }
    *line_number += 1;
    if buffer.last() == Some(&b'\n') {
        buffer.pop();
    }
    if buffer.last() == Some(&b'\r') {
        buffer.pop();
    }
    Ok(true)
}

fn validate_header(bytes: &[u8], line: u64) -> Result<()> {
    if bytes
        .iter()
        .copied()
        .any(|byte| !is_valid_header_byte(byte))
    {
        return Err(invalid_at(line, "invalid byte in sequence identifier"));
    }
    Ok(())
}

fn validate_sequence(bytes: &[u8], line: u64) -> Result<()> {
    if bytes
        .iter()
        .copied()
        .any(|byte| !is_valid_sequence_byte(byte))
    {
        return Err(invalid_at(
            line,
            "sequence contains a byte outside printable non-space ASCII 33..=126",
        ));
    }
    Ok(())
}

fn validate_quality(bytes: &[u8], line: u64) -> Result<()> {
    if bytes
        .iter()
        .copied()
        .any(|byte| !is_valid_quality_byte(byte))
    {
        return Err(invalid_at(
            line,
            "FASTQ quality contains a byte outside ASCII 33..=126",
        ));
    }
    Ok(())
}

fn invalid_at(line: u64, message: &str) -> RsomicsError {
    RsomicsError::InvalidInput(format!("line {line}: {message}"))
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;

    #[test]
    fn multiline_fasta_borrows_reused_buffers() {
        let input = b">one description\nACGT\nTGCA\n>two\nNNNN\n";
        let mut reader = Reader::detect(Cursor::new(input)).unwrap();
        assert_eq!(reader.format(), Format::Fasta);

        let first = reader.read_record().unwrap().unwrap().to_owned();
        assert_eq!(first.id, b"one description");
        assert_eq!(first.seq, b"ACGTTGCA");
        assert!(first.qual.is_none());

        let second = reader.read_record().unwrap().unwrap();
        assert_eq!(second.id, b"two");
        assert_eq!(second.seq, b"NNNN");
        assert!(reader.read_record().unwrap().is_none());
    }

    #[test]
    fn fasta_crlf_is_normalized() {
        let mut reader = Reader::new(Cursor::new(b">one\r\nAC\r\nGT\r\n"), Format::Fasta);
        let record = reader.read_record().unwrap().unwrap();
        assert_eq!(record.id, b"one");
        assert_eq!(record.seq, b"ACGT");
    }

    #[test]
    fn fasta_header_without_sequence_errors() {
        let mut reader = Reader::new(Cursor::new(b">one\n"), Format::Fasta);
        assert!(matches!(
            reader.read_record(),
            Err(RsomicsError::InvalidInput(_))
        ));
    }

    #[test]
    fn strict_fastq_reads_four_lines() {
        let input = b"@r1\nACGT\n+r1\nIIII\n@r2\nTGCA\n+\nFFFF\n";
        let mut reader = Reader::detect(Cursor::new(input)).unwrap();
        assert_eq!(reader.format(), Format::Fastq);
        assert_eq!(
            reader.read_record().unwrap().unwrap().to_owned(),
            crate::OwnedRecord {
                id: b"r1".to_vec(),
                seq: b"ACGT".to_vec(),
                qual: Some(b"IIII".to_vec()),
            }
        );
        assert_eq!(reader.read_record().unwrap().unwrap().id, b"r2");
        assert!(reader.read_record().unwrap().is_none());
    }

    #[test]
    fn fastq_crlf_is_normalized() {
        let mut reader = Reader::new(Cursor::new(b"@r1\r\nACGT\r\n+\r\nIIII\r\n"), Format::Fastq);
        let record = reader.read_record().unwrap().unwrap();
        assert_eq!(record.id, b"r1");
        assert_eq!(record.seq, b"ACGT");
        assert_eq!(record.qual, Some(b"IIII".as_slice()));
    }

    #[test]
    fn wrapped_fastq_is_rejected() {
        let mut reader = Reader::new(Cursor::new(b"@r1\nAC\nGT\n+\nIIII\n"), Format::Fastq);
        assert!(matches!(
            reader.read_record(),
            Err(RsomicsError::InvalidInput(_))
        ));
    }

    #[test]
    fn truncated_fastq_is_rejected() {
        let mut reader = Reader::new(Cursor::new(b"@r1\nACGT\n+\n"), Format::Fastq);
        assert!(matches!(
            reader.read_record(),
            Err(RsomicsError::InvalidInput(_))
        ));
    }

    #[test]
    fn fastq_length_mismatch_is_rejected() {
        let mut reader = Reader::new(Cursor::new(b"@r1\nACGT\n+\nIII\n"), Format::Fastq);
        assert!(matches!(
            reader.read_record(),
            Err(RsomicsError::InvalidInput(_))
        ));
    }

    #[test]
    fn fastq_repeated_identifier_must_match() {
        let mut reader = Reader::new(
            Cursor::new(b"@r1 description\nACGT\n+r2\nIIII\n"),
            Format::Fastq,
        );
        assert!(matches!(
            reader.read_record(),
            Err(RsomicsError::InvalidInput(_))
        ));
    }

    #[test]
    fn empty_explicit_reader_is_clean_eof() {
        let mut reader = Reader::new(Cursor::new(b""), Format::Fastq);
        assert!(reader.read_record().unwrap().is_none());
    }

    #[test]
    fn empty_detected_reader_is_an_error() {
        assert!(matches!(
            Reader::detect(Cursor::new(b"")),
            Err(RsomicsError::InvalidInput(_))
        ));
    }

    #[test]
    fn complete_records_need_no_final_newline() {
        let mut fasta = Reader::new(Cursor::new(b">one\nACGT"), Format::Fasta);
        assert_eq!(fasta.read_record().unwrap().unwrap().seq, b"ACGT");
        assert!(fasta.read_record().unwrap().is_none());

        let mut fastq = Reader::new(Cursor::new(b"@one\nACGT\n+\nIIII"), Format::Fastq);
        assert_eq!(
            fastq.read_record().unwrap().unwrap().qual,
            Some(b"IIII".as_slice())
        );
        assert!(fastq.read_record().unwrap().is_none());
    }

    #[test]
    fn empty_fastq_sequence_and_quality_are_valid() {
        let mut reader = Reader::new(Cursor::new(b"@one\n\n+\n\n"), Format::Fastq);
        let record = reader.read_record().unwrap().unwrap();
        assert!(record.seq.is_empty());
        assert_eq!(record.qual, Some(b"".as_slice()));
    }

    #[test]
    fn sequence_accepts_biological_symbols_and_rejects_controls() {
        let mut valid = Reader::new(Cursor::new(b">one\nACGTN-.*?\n"), Format::Fasta);
        assert_eq!(valid.read_record().unwrap().unwrap().seq, b"ACGTN-.*?");

        for invalid in [b">one\nA\0C\n".as_slice(), b">one\nA\x7fC\n"] {
            let mut reader = Reader::new(Cursor::new(invalid), Format::Fasta);
            assert!(matches!(
                reader.read_record(),
                Err(RsomicsError::InvalidInput(_))
            ));
        }
    }

    #[test]
    fn reader_is_fail_closed_after_record_error() {
        let input = b"@bad\nACGT\n+\nIII\n@good\nACGT\n+\nIIII\n";
        let mut reader = Reader::new(Cursor::new(input), Format::Fastq);
        assert!(matches!(
            reader.read_record(),
            Err(RsomicsError::InvalidInput(_))
        ));
        let second = reader.read_record().unwrap_err();
        assert!(
            second
                .to_string()
                .contains("cannot continue after a previous record error")
        );
    }
}
