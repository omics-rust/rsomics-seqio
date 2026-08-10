use std::io::{self, Write};

use noodles_bgzf::io::{Writer as BgzfWriter, writer};
use rsomics_common::{Result, RsomicsError};

use crate::{Format, Record, Writer};

/// Encoding applied after FASTA or FASTQ serialization.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Compression {
    Plain,
    Bgzf { level: u8 },
}

enum Encoder<W: Write> {
    Plain(W),
    Bgzf(BgzfWriter<W>),
}

impl<W: Write> Encoder<W> {
    fn new(inner: W, compression: Compression) -> Result<Self> {
        match compression {
            Compression::Plain => Ok(Self::Plain(inner)),
            Compression::Bgzf { level } => {
                if level > 9 {
                    return Err(RsomicsError::InvalidInput(format!(
                        "BGZF compression level must be between 0 and 9, got {level}"
                    )));
                }
                let level = writer::CompressionLevel::new(level).unwrap();
                let writer = writer::Builder::default()
                    .set_compression_level(level)
                    .build_from_writer(inner);
                Ok(Self::Bgzf(writer))
            }
        }
    }

    fn finish(self) -> Result<W> {
        match self {
            Self::Plain(mut inner) => {
                inner.flush().map_err(RsomicsError::Io)?;
                Ok(inner)
            }
            Self::Bgzf(inner) => inner.finish().map_err(RsomicsError::Io),
        }
    }
}

impl<W: Write> Write for Encoder<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Self::Plain(inner) => inner.write(buf),
            Self::Bgzf(inner) => inner.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::Plain(inner) => inner.flush(),
            Self::Bgzf(inner) => inner.flush(),
        }
    }
}

/// A strict FASTA/FASTQ writer with explicit stream finalization.
pub struct OutputWriter<W: Write> {
    inner: Writer<Encoder<W>>,
}

impl<W: Write> OutputWriter<W> {
    /// Creates a writer over a caller-owned destination.
    pub fn new(inner: W, format: Format, compression: Compression) -> Result<Self> {
        Ok(Self {
            inner: Writer::new(Encoder::new(inner, compression)?, format),
        })
    }

    #[must_use]
    pub fn format(&self) -> Format {
        self.inner.format()
    }

    pub fn write_record(&mut self, record: Record<'_>) -> Result<()> {
        self.inner.write_record(record)
    }

    /// Flushes records, completes compression, and returns the destination.
    pub fn finish(self) -> Result<W> {
        self.inner.finish_into_inner()?.finish()
    }
}
