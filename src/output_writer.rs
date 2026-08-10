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

/// A plain or BGZF byte-stream encoder over a caller-owned destination.
pub struct OutputEncoder<W: Write> {
    inner: Encoder<W>,
}

enum Encoder<W: Write> {
    Plain(W),
    Bgzf(BgzfWriter<W>),
}

impl<W: Write> OutputEncoder<W> {
    pub fn new(inner: W, compression: Compression) -> Result<Self> {
        let inner = match compression {
            Compression::Plain => Encoder::Plain(inner),
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
                Encoder::Bgzf(writer)
            }
        };
        Ok(Self { inner })
    }

    pub fn finish(self) -> Result<W> {
        match self.inner {
            Encoder::Plain(mut inner) => {
                inner.flush().map_err(RsomicsError::Io)?;
                Ok(inner)
            }
            Encoder::Bgzf(inner) => inner.finish().map_err(RsomicsError::Io),
        }
    }
}

impl<W: Write> Write for OutputEncoder<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match &mut self.inner {
            Encoder::Plain(inner) => inner.write(buf),
            Encoder::Bgzf(inner) => inner.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match &mut self.inner {
            Encoder::Plain(inner) => inner.flush(),
            Encoder::Bgzf(inner) => inner.flush(),
        }
    }
}

/// A strict FASTA/FASTQ writer with explicit stream finalization.
pub struct OutputWriter<W: Write> {
    inner: Writer<OutputEncoder<W>>,
}

impl<W: Write> OutputWriter<W> {
    /// Creates a writer over a caller-owned destination.
    pub fn new(inner: W, format: Format, compression: Compression) -> Result<Self> {
        Ok(Self {
            inner: Writer::new(OutputEncoder::new(inner, compression)?, format),
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
