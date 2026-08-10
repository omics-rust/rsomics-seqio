use std::io::{self, Cursor, Write};

use rsomics_seqio::{Compression, Format, OutputWriter, Record, open_reader};

const BGZF_EOF: &[u8] = &[
    0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0, 0x1b, 0, 3, 0, 0, 0, 0, 0,
    0, 0, 0, 0,
];

#[test]
fn plain_output_returns_the_caller_stream() {
    let mut writer = OutputWriter::new(Vec::new(), Format::Fasta, Compression::Plain).unwrap();
    writer
        .write_record(Record {
            id: b"read-1",
            seq: b"ACGT",
            qual: None,
        })
        .unwrap();

    assert_eq!(writer.finish().unwrap(), b">read-1\nACGT\n");
}

#[test]
fn bgzf_output_round_trips_and_has_one_eof_marker() {
    let mut writer =
        OutputWriter::new(Vec::new(), Format::Fastq, Compression::Bgzf { level: 6 }).unwrap();
    writer
        .write_record(Record {
            id: b"read-1",
            seq: b"ACGT",
            qual: Some(b"IJKL"),
        })
        .unwrap();

    let encoded = writer.finish().unwrap();
    assert!(encoded.ends_with(BGZF_EOF));
    assert_eq!(
        encoded
            .windows(BGZF_EOF.len())
            .filter(|w| *w == BGZF_EOF)
            .count(),
        1
    );

    let mut reader = open_reader(Cursor::new(encoded)).unwrap();
    let record = reader.read_record().unwrap().unwrap();
    assert_eq!(record.id, b"read-1");
    assert_eq!(record.seq, b"ACGT");
    assert_eq!(record.qual, Some(b"IJKL".as_slice()));
    assert!(reader.read_record().unwrap().is_none());
}

#[test]
fn invalid_bgzf_level_is_rejected() {
    let error = OutputWriter::new(Vec::new(), Format::Fasta, Compression::Bgzf { level: 10 })
        .err()
        .unwrap();
    assert!(error.to_string().contains("compression level"));
}

#[derive(Debug)]
struct FailingSink;

impl Write for FailingSink {
    fn write(&mut self, _: &[u8]) -> io::Result<usize> {
        Err(io::Error::new(io::ErrorKind::BrokenPipe, "closed"))
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

#[test]
fn bgzf_finalization_failure_is_returned() {
    let mut writer =
        OutputWriter::new(FailingSink, Format::Fasta, Compression::Bgzf { level: 6 }).unwrap();
    writer
        .write_record(Record {
            id: b"read-1",
            seq: b"ACGT",
            qual: None,
        })
        .unwrap();

    let error = writer.finish().unwrap_err();
    assert!(matches!(
        error,
        rsomics_common::RsomicsError::Io(ref source)
            if source.kind() == io::ErrorKind::BrokenPipe
    ));
}
