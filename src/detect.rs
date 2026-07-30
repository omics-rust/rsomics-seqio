use std::io::{Cursor, Read};

use rsomics_common::{Result, RsomicsError};

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CompressionKind {
    Plain,
    Gzip,
}

pub(crate) struct ReplayReader<R> {
    prefix: Cursor<Vec<u8>>,
    inner: R,
}

impl<R: Read> Read for ReplayReader<R> {
    fn read(&mut self, output: &mut [u8]) -> std::io::Result<usize> {
        let replayed = self.prefix.read(output)?;
        if replayed != 0 {
            return Ok(replayed);
        }
        self.inner.read(output)
    }
}

pub(crate) fn probe<R: Read>(mut source: R) -> Result<(CompressionKind, ReplayReader<R>)> {
    let mut prefix = Vec::with_capacity(GZIP_MAGIC.len());
    while prefix.len() < GZIP_MAGIC.len() {
        let mut byte = [0u8; 1];
        match source.read(&mut byte) {
            Ok(0) => break,
            Ok(_) => prefix.push(byte[0]),
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => {}
            Err(error) => return Err(RsomicsError::Io(error)),
        }
    }

    let kind = if prefix == GZIP_MAGIC {
        CompressionKind::Gzip
    } else {
        CompressionKind::Plain
    };
    Ok((
        kind,
        ReplayReader {
            prefix: Cursor::new(prefix),
            inner: source,
        },
    ))
}

// BGZF uses the gzip decoder too. This parser is retained only for tests and
// specification checks; runtime routing needs no more than gzip's two-byte magic.
#[cfg(test)]
fn is_bgzf_header(probe: &[u8]) -> bool {
    if probe.len() < 12 || probe[..2] != GZIP_MAGIC || probe[2] != 0x08 || probe[3] & 0x04 == 0 {
        return false;
    }
    let extra_len = usize::from(u16::from_le_bytes([probe[10], probe[11]]));
    let Some(extra_end) = 12usize.checked_add(extra_len) else {
        return false;
    };
    if probe.len() < extra_end {
        return false;
    }

    let mut offset = 12;
    while offset + 4 <= extra_end {
        let field_len = usize::from(u16::from_le_bytes([probe[offset + 2], probe[offset + 3]]));
        let Some(field_end) = offset.checked_add(4 + field_len) else {
            return false;
        };
        if field_end > extra_end {
            return false;
        }
        if probe[offset] == b'B' && probe[offset + 1] == b'C' && field_len == 2 {
            return true;
        }
        offset = field_end;
    }
    false
}

#[cfg(test)]
mod tests {
    use std::io::{Cursor, Read};

    use super::*;

    const BGZF_HEADER: &[u8] = &[
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, b'B', b'C', 0x02,
        0x00, 0x00, 0x00,
    ];

    const BGZF_HEADER_WITH_PRIOR_EXTRA: &[u8] = &[
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x0b, 0x00, b'X', b'Y', 0x01,
        0x00, 0xff, b'B', b'C', 0x02, 0x00, 0x00, 0x00,
    ];

    struct OneByte<R>(R);

    impl<R: Read> Read for OneByte<R> {
        fn read(&mut self, output: &mut [u8]) -> std::io::Result<usize> {
            let limit = output.len().min(1);
            self.0.read(&mut output[..limit])
        }
    }

    #[test]
    fn probe_accumulates_and_replays_short_reads() {
        for (bytes, expected) in [
            (b">record".as_slice(), CompressionKind::Plain),
            (&[0x1f, 0x8b, 0x08, 0x00], CompressionKind::Gzip),
        ] {
            let (kind, mut replayed) = probe(OneByte(Cursor::new(bytes))).unwrap();
            assert_eq!(kind, expected);
            let mut recovered = Vec::new();
            replayed.read_to_end(&mut recovered).unwrap();
            assert_eq!(recovered, bytes);
        }
    }

    #[test]
    fn one_byte_plain_input_is_replayed() {
        let (kind, mut replayed) = probe(OneByte(Cursor::new(b">"))).unwrap();
        assert_eq!(kind, CompressionKind::Plain);
        let mut recovered = Vec::new();
        replayed.read_to_end(&mut recovered).unwrap();
        assert_eq!(recovered, b">");
    }

    #[test]
    fn bgzf_headers_remain_recognizable_as_gzip_extensions() {
        assert!(is_bgzf_header(BGZF_HEADER));
        assert!(is_bgzf_header(BGZF_HEADER_WITH_PRIOR_EXTRA));
        assert!(!is_bgzf_header(&[0x1f, 0x8b, 0x08, 0x00]));
    }
}
