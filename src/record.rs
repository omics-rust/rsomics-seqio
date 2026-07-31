#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Record<'a> {
    pub id: &'a [u8],
    pub seq: &'a [u8],
    pub qual: Option<&'a [u8]>,
}

impl Record<'_> {
    #[must_use]
    pub fn to_owned(self) -> OwnedRecord {
        OwnedRecord {
            id: self.id.to_vec(),
            seq: self.seq.to_vec(),
            qual: self.qual.map(<[u8]>::to_vec),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OwnedRecord {
    pub id: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>,
}

impl OwnedRecord {
    #[must_use]
    pub fn as_record(&self) -> Record<'_> {
        Record {
            id: &self.id,
            seq: &self.seq,
            qual: self.qual.as_deref(),
        }
    }
}

pub(crate) fn is_valid_header_byte(byte: u8) -> bool {
    byte == b'\t' || (b' '..=b'~').contains(&byte)
}

pub(crate) fn is_valid_printable_byte(byte: u8) -> bool {
    (b'!'..=b'~').contains(&byte)
}

pub(crate) fn are_valid_printable_bytes(bytes: &[u8]) -> bool {
    const BYTE_ONES: u64 = 0x0101_0101_0101_0101;
    const BYTE_HIGHS: u64 = 0x8080_8080_8080_8080;
    const BYTE_DELS: u64 = 0x7f7f_7f7f_7f7f_7f7f;

    let mut chunks = bytes.chunks_exact(8);
    for chunk in &mut chunks {
        let word = u64::from_ne_bytes(chunk.try_into().expect("chunk length is exactly eight"));
        let below_exclamation = word.wrapping_sub(BYTE_ONES * u64::from(b'!')) & !word & BYTE_HIGHS;
        let del = word ^ BYTE_DELS;
        let contains_del = del.wrapping_sub(BYTE_ONES) & !del & BYTE_HIGHS;

        // A lane is printable non-space ASCII when its high bit is clear, it
        // is not DEL, and the lane-wise subtraction finds nothing below '!'.
        // Cross-byte borrows can smear the mask location, so these masks are
        // used only as an any-invalid predicate, never to locate a bad byte.
        if word & BYTE_HIGHS != 0 || below_exclamation != 0 || contains_del != 0 {
            return false;
        }
    }
    chunks
        .remainder()
        .iter()
        .copied()
        .all(is_valid_printable_byte)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn word_scan_matches_scalar_boundaries_across_chunks() {
        for position in 0..16 {
            for first in u8::MIN..=u8::MAX {
                for second in u8::MIN..=u8::MAX {
                    let mut bytes = [b'A'; 17];
                    bytes[position] = first;
                    bytes[position + 1] = second;
                    assert_eq!(
                        are_valid_printable_bytes(&bytes),
                        bytes.iter().copied().all(is_valid_printable_byte),
                        "position={position}, first={first}, second={second}"
                    );
                }
            }
        }
    }
}
