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

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LegacyFastqRecord {
    pub id: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

pub(crate) fn is_valid_header_byte(byte: u8) -> bool {
    (b' '..=b'~').contains(&byte)
}

pub(crate) fn is_valid_sequence_byte(byte: u8) -> bool {
    (b'!'..=b'~').contains(&byte)
}

pub(crate) fn is_valid_quality_byte(byte: u8) -> bool {
    (b'!'..=b'~').contains(&byte)
}
