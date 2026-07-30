use std::io::{Cursor, Write};

use rsomics_seqio::{Format, Reader, Record, Writer, open_path, open_reader};

#[test]
fn tabbed_fasta_and_fastq_headers_round_trip() {
    for (format, input, expected_quality) in [
        (Format::Fasta, b">seq1\talpha beta\nACGT\n".as_slice(), None),
        (
            Format::Fastq,
            b"@seq1\talpha beta\nACGT\n+seq1\talpha beta\nIIII\n".as_slice(),
            Some(b"IIII".as_slice()),
        ),
    ] {
        let mut reader = Reader::new(Cursor::new(input), format);
        let record = reader.read_record().unwrap().unwrap().to_owned();
        assert_eq!(record.id, b"seq1\talpha beta");
        assert_eq!(record.qual.as_deref(), expected_quality);

        let mut output = Vec::new();
        let mut writer = Writer::new(&mut output, format);
        writer.write_owned(&record).unwrap();
        writer.finish().unwrap();

        let mut round_trip = Reader::new(Cursor::new(output), format);
        assert_eq!(
            round_trip.read_record().unwrap().unwrap().to_owned(),
            record
        );
    }
}

#[test]
fn wrapped_fastq_preserves_borrowed_record_shape() {
    let input = b"@r1\twrapped\nAC\nGT\n+r1\twrapped\nII\nII\n@r2\nTGCA\n+\nFFFF\n";
    let mut reader = open_reader(Cursor::new(input)).unwrap();

    assert_eq!(
        reader.read_record().unwrap().unwrap().to_owned(),
        rsomics_seqio::OwnedRecord {
            id: b"r1\twrapped".to_vec(),
            seq: b"ACGT".to_vec(),
            qual: Some(b"IIII".to_vec()),
        }
    );
    let second = reader.read_record().unwrap().unwrap();
    assert_eq!(second.id, b"r2");
    assert_eq!(second.seq, b"TGCA");
    assert_eq!(second.qual, Some(b"FFFF".as_slice()));
    assert!(reader.read_record().unwrap().is_none());
}

#[test]
fn writer_output_with_plus_prefixed_sequence_round_trips() {
    let record = Record {
        id: b"plus",
        seq: b"+ACG",
        qual: Some(b"IIII"),
    };
    let mut output = Vec::new();
    let mut writer = Writer::new(&mut output, Format::Fastq);
    writer.write_record(record).unwrap();
    writer.finish().unwrap();

    let mut reader = Reader::new(Cursor::new(output), Format::Fastq);
    assert_eq!(reader.read_record().unwrap().unwrap(), record);
    assert!(reader.read_record().unwrap().is_none());
}

#[test]
fn valid_prefix_does_not_hide_trailing_fastq_damage() {
    let mut input = Vec::new();
    for index in 0..64 {
        writeln!(input, "@read{index}\nAC\nGT\n+\nII\nII").unwrap();
    }
    input.extend_from_slice(b"@broken\nACGT\n+\nIII\n");

    let mut generic = open_reader(Cursor::new(&input)).unwrap();
    for _ in 0..64 {
        assert!(generic.read_record().unwrap().is_some());
    }
    assert!(generic.read_record().is_err());

    let mut file = tempfile::NamedTempFile::new().unwrap();
    file.write_all(&input).unwrap();
    file.flush().unwrap();
    let mut path = open_path(file.path()).unwrap();
    for _ in 0..64 {
        assert!(path.read_record().unwrap().is_some());
    }
    assert!(path.read_record().is_err());
}

#[test]
fn public_writer_rejects_non_ascii_and_control_headers() {
    for id in [
        b"bad\x0bheader".as_slice(),
        b"bad\x7fheader",
        b"bad\x80header",
    ] {
        let mut writer = Writer::new(Vec::new(), Format::Fasta);
        assert!(
            writer
                .write_record(Record {
                    id,
                    seq: b"ACGT",
                    qual: None,
                })
                .is_err()
        );
    }
}
