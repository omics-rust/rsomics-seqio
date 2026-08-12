use std::fs;
use std::io::Write;

use noodles_bgzf as bgzf;
use rsomics_seqio::IndexedFasta;

fn write_reference(path: &std::path::Path) {
    fs::write(path, b">chr1\nACGTAC\nGTAA\n>chr2\nTTGG\n").unwrap();
    fs::write(
        path.with_extension("fa.fai"),
        b"chr1\t10\t6\t6\t7\nchr2\t4\t24\t4\t5\n",
    )
    .unwrap();
}

#[test]
fn zero_based_ranges_cross_fasta_lines_and_reuse_the_cache() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("reference.fa");
    write_reference(&path);

    let mut reference = IndexedFasta::with_cache_capacity(&path, 4).unwrap();
    assert_eq!(reference.len(b"chr1").unwrap(), 10);
    assert_eq!(reference.fetch(b"chr1", 4..9).unwrap(), b"ACGTA");
    assert_eq!(reference.fetch(b"chr1", 5..7).unwrap(), b"CG");
    assert_eq!(reference.fetch(b"chr2", 1..4).unwrap(), b"TGG");
}

#[test]
fn invalid_name_and_ranges_fail_with_reference_context() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("reference.fa");
    write_reference(&path);

    let mut reference = IndexedFasta::open(&path).unwrap();
    let missing = reference.fetch(b"absent", 0..1).unwrap_err().to_string();
    assert!(missing.contains("absent"), "{missing}");
    assert!(missing.contains("reference.fa"), "{missing}");

    let reversed_range = std::ops::Range { start: 4, end: 3 };
    let reversed = reference
        .fetch(b"chr1", reversed_range)
        .unwrap_err()
        .to_string();
    assert!(reversed.contains("4..3"), "{reversed}");
    let past_end = reference.fetch(b"chr1", 9..11).unwrap_err().to_string();
    assert!(past_end.contains("length 10"), "{past_end}");
}

#[test]
fn bgzf_reference_uses_fai_and_gzi_indexes() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("reference.fa.gz");
    let mut writer = bgzf::io::Writer::new(fs::File::create(&path).unwrap());
    writer
        .write_all(b">chr1\nACGTAC\nGTAA\n>chr2\nTTGG\n")
        .unwrap();
    writer.try_finish().unwrap();
    fs::write(
        path.with_extension("gz.fai"),
        b"chr1\t10\t6\t6\t7\nchr2\t4\t24\t4\t5\n",
    )
    .unwrap();
    fs::write(path.with_extension("gz.gzi"), 0_u64.to_le_bytes()).unwrap();

    let mut reference = IndexedFasta::open(&path).unwrap();
    assert_eq!(reference.fetch(b"chr1", 4..9).unwrap(), b"ACGTA");
    assert_eq!(reference.fetch(b"chr2", 0..4).unwrap(), b"TTGG");
}

#[test]
fn duplicate_reference_names_are_rejected() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("reference.fa");
    fs::write(&path, b">chr1\nAC\n>chr1\nGT\n").unwrap();
    fs::write(
        path.with_extension("fa.fai"),
        b"chr1\t2\t6\t2\t3\nchr1\t2\t15\t2\t3\n",
    )
    .unwrap();

    let error = IndexedFasta::open(&path).err().unwrap().to_string();
    assert!(error.contains("duplicate reference chr1"), "{error}");
}
