# rsomics-seqio

Strict, streaming FASTA and FASTQ I/O for the rsomics product family.

## Public API

`Reader<R: BufRead>` is the primary boundary. `read_record` lends a
`Record<'_>` from reusable internal buffers, so a streaming consumer does not
allocate an identifier, sequence, and quality vector for every record.

```rust
use std::io::{BufReader, stdin};
use rsomics_seqio::Reader;

let input = stdin().lock();
let mut reader = Reader::detect(BufReader::new(input))?;
while let Some(record) = reader.read_record()? {
    consume(record.id, record.seq, record.qual);
}
```

`open_path` detects plain, gzip, and BGZF input from content, then detects
FASTA or FASTQ after decompression. `open_reader` provides the same transparent
compression and format detection for a generic `Read` source. Compression
probing consumes only the two-byte gzip magic and replays those bytes, so it is
safe for sources that return short reads. Consumers that already know the
format can construct `Reader` directly.

`Writer<W: Write>` writes canonical single-line FASTA or strict four-line
FASTQ to any writer. Compression and transactional file policy stay with the
consuming product. Call `finish` or `finish_into_inner` to flush the sink.

`Record::to_owned` produces the canonical `OwnedRecord`, whose optional
quality field represents both FASTA and FASTQ without a second ambiguous owned
record type.

## Validation

- FASTA supports multiline sequences and requires a non-empty identifier and
  at least one non-empty sequence line.
- FASTQ accepts single-line or wrapped sequence and quality bodies. The `+`
  separator ends the sequence body; quality lines are accumulated until their
  byte count exactly equals the sequence length. Repeated identifiers after
  `+` must match, and short, overlong, or truncated qualities are errors.
- Headers accept printable ASCII plus TAB. CR and LF remain line boundaries;
  other controls and non-ASCII bytes are rejected.
- Sequence bytes must be printable, non-space ASCII (`!` through `~`). This
  includes common biological symbols such as letters, `-`, `.`, `*`, and `?`,
  while rejecting spaces, NUL, DEL, other ASCII controls, and non-ASCII bytes.
- Empty FASTQ sequence and quality lines are accepted when both lengths are
  zero. FASTA still requires at least one non-empty sequence line.
- CRLF and LF inputs are normalized.
- Empty input is clean EOF when the caller supplied a `Format`; automatic
  format detection rejects it because no format can be inferred.
- Malformed records, truncated gzip streams, writer failures, and
  decompressor-thread failures are errors rather than short successful
  operations. A reader is fail-closed after its first record error; callers
  must discard it rather than attempt to resume at an uncertain boundary.

## Compression

Path gzip and BGZF input uses a private flate2/zlib-rs decoder producer thread.
The file is opened once, probed through that handle, and replayed into the
decoder, avoiding a path re-open race. Generic readers and gzip writers use
the same decoder family. BGZF is recognized as gzip framing and intentionally
uses the same decoder route rather than a separate runtime classification.
Backend types are not public API.

Third-party license metadata:

- flate2: MIT OR Apache-2.0
- zlib-rs: Zlib
- crossbeam-channel: MIT OR Apache-2.0

## Compatibility

Version 0.4 uses the current `rsomics-common` error contract. Its record and
I/O APIs retain the version 0.3 behavior described below.

Version 0.2 is an intentional source-breaking redesign and does not claim 0.1
API compatibility:

- `OwnedFastxRecord` is replaced by the canonical
  `OwnedRecord { id, seq, qual: Option<Vec<u8>> }`.
- Sequence validation is stricter and consistent between reader and writer:
  bytes outside printable non-space ASCII are rejected.
- gzip/BGZF decoding no longer selects `rsomics-igzip` on Linux.

New product code should use `open_path`, `open_reader`, or `Reader` directly.

## Inspected team-owned source assets

The redesign reviewed these repository revisions:

- rsomics-seqio: `979b609cb87dbe468c122f055148822250521746`
- rsomics-fqgz: `c5e7de12d21e72cdcc3b62f84302653c34b5dc54`
- rsomics-fasta-utils: `93e81c2ab88524f97bdaaab6f34105743d798b96`
- rsomics-fastq-utils: `7f0551d463977e6343af3c9477879799d1bcbc81`
- rsomics-fasta-validate: `93da2fdc4f8596899e76822271ab6a718df519e6`
- rsomics-fastq-validate: `69f6af79d6af37a6384ac3af3938062d11c60ae7`
- rsomics-fastx-convert: `6ac272a61a7017ce0bf26b2209fbea113fb561a1`

All listed rsomics sources are team-owned. No third-party implementation was
copied.

## Origin

The path decoder thread and large-block buffering follow the architecture
described by fastp (Chen et al. 2018, doi:10.1093/bioinformatics/bty560;
upstream MIT). BGZF framing follows SAM/BAM specification section 4.1.

License: MIT OR Apache-2.0.
