use std::fs;
use std::hint::black_box;
use std::io::Write;

use criterion::{Criterion, criterion_group, criterion_main};
use noodles_core::{Position, Region};
use noodles_fasta as fasta;
use rsomics_seqio::IndexedFasta;

fn indexed_ranges(criterion: &mut Criterion) {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("reference.fa");
    let sequence = vec![b'A'; 4 * 1024 * 1024];
    let mut file = fs::File::create(&path).unwrap();
    file.write_all(b">chr1\n").unwrap();
    for line in sequence.chunks(80) {
        file.write_all(line).unwrap();
        file.write_all(b"\n").unwrap();
    }
    fs::write(
        path.with_extension("fa.fai"),
        format!("chr1\t{}\t6\t80\t81\n", sequence.len()),
    )
    .unwrap();

    let mut cached = IndexedFasta::open(&path).unwrap();
    criterion.bench_function("indexed_fasta_cached_64b", |bencher| {
        bencher.iter(|| {
            black_box(cached.fetch(b"chr1", 2_000_000..2_000_064).unwrap());
        });
    });

    let mut direct = fasta::io::indexed_reader::Builder::default()
        .build_from_path(&path)
        .unwrap();
    let start = Position::try_from(2_000_001).unwrap();
    let end = Position::try_from(2_000_064).unwrap();
    let region = Region::new("chr1", start..=end);
    criterion.bench_function("noodles_indexed_fasta_direct_64b", |bencher| {
        bencher.iter(|| black_box(direct.query(&region).unwrap()));
    });
}

criterion_group!(benches, indexed_ranges);
criterion_main!(benches);
