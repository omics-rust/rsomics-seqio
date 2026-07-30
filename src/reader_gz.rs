use std::collections::VecDeque;
use std::io::{BufReader, Cursor, Read};
use std::thread;

use crossbeam_channel::{Receiver, bounded};
use flate2::read::MultiGzDecoder;
use rsomics_common::{Result, RsomicsError};

const DECODE_BUFFER: usize = 8 * 1024 * 1024;
const CHANNEL_DEPTH: usize = 4;

pub(crate) struct GzipStream {
    receiver: Receiver<std::io::Result<Vec<u8>>>,
    chunks: VecDeque<Cursor<Vec<u8>>>,
    finished: bool,
}

impl GzipStream {
    pub(crate) fn new<R>(source: R) -> Result<Self>
    where
        R: Read + Send + 'static,
    {
        let (sender, receiver) = bounded(CHANNEL_DEPTH);
        thread::Builder::new()
            .name("seqio-gzip".into())
            .spawn(move || produce(source, &sender))
            .map_err(RsomicsError::Io)?;
        Ok(Self {
            receiver,
            chunks: VecDeque::new(),
            finished: false,
        })
    }

    fn receive_chunk(&mut self) -> std::io::Result<bool> {
        if self.finished {
            return Ok(false);
        }
        match self.receiver.recv() {
            Ok(Ok(chunk)) if chunk.is_empty() => {
                self.finished = true;
                Ok(false)
            }
            Ok(Ok(chunk)) => {
                self.chunks.push_back(Cursor::new(chunk));
                Ok(true)
            }
            Ok(Err(error)) => {
                self.finished = true;
                Err(error)
            }
            Err(_) => {
                self.finished = true;
                Err(std::io::Error::new(
                    std::io::ErrorKind::UnexpectedEof,
                    "gzip decode thread ended without an EOF marker",
                ))
            }
        }
    }
}

impl Read for GzipStream {
    fn read(&mut self, output: &mut [u8]) -> std::io::Result<usize> {
        if output.is_empty() {
            return Ok(0);
        }
        loop {
            if let Some(chunk) = self.chunks.front_mut() {
                let read = chunk.read(output)?;
                if read > 0 {
                    return Ok(read);
                }
                self.chunks.pop_front();
                continue;
            }
            if !self.receive_chunk()? {
                return Ok(0);
            }
        }
    }
}

fn produce<R>(source: R, sender: &crossbeam_channel::Sender<std::io::Result<Vec<u8>>>)
where
    R: Read,
{
    let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        produce_inner(source, sender)
    }));
    let error = match outcome {
        Ok(Ok(())) => return,
        Ok(Err(error)) => error,
        Err(_) => std::io::Error::other("gzip decode thread panicked"),
    };
    let _ = sender.send(Err(error));
}

fn produce_inner<R>(
    source: R,
    sender: &crossbeam_channel::Sender<std::io::Result<Vec<u8>>>,
) -> std::io::Result<()>
where
    R: Read,
{
    let decoder = MultiGzDecoder::new(source);
    let mut decoder = BufReader::with_capacity(DECODE_BUFFER, decoder);
    loop {
        let mut block = vec![0u8; DECODE_BUFFER];
        let read = decoder.read(&mut block)?;
        if read == 0 {
            let _ = sender.send(Ok(Vec::new()));
            return Ok(());
        }
        block.truncate(read);
        if sender.send(Ok(block)).is_err() {
            return Ok(());
        }
    }
}

#[cfg(test)]
mod tests {
    use std::io::{Cursor, Write};

    use super::*;

    fn gzip_member(bytes: &[u8]) -> Vec<u8> {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(bytes).unwrap();
        encoder.finish().unwrap()
    }

    #[test]
    fn stream_decodes_multiple_members() {
        let mut encoded = gzip_member(b">one\nAC");
        encoded.extend_from_slice(&gzip_member(b"GT\n>two\nTGCA\n"));

        let mut decoded = Vec::new();
        GzipStream::new(Cursor::new(encoded))
            .unwrap()
            .read_to_end(&mut decoded)
            .unwrap();
        assert_eq!(decoded, b">one\nACGT\n>two\nTGCA\n");
    }

    #[test]
    fn truncated_gzip_errors_loudly() {
        let mut encoded = gzip_member(b"@r1\nACGT\n+\nIIII\n");
        encoded.truncate(encoded.len() - 6);

        let mut decoded = Vec::new();
        assert!(
            GzipStream::new(Cursor::new(encoded))
                .unwrap()
                .read_to_end(&mut decoded)
                .is_err()
        );
    }
}
