use bio::io::fasta::IndexedReader;
use std::path::Path;

use thiserror;

#[derive(thiserror::Error, Debug)]
pub enum GenomeError {
    #[error("Cannot read indexed fasta. Message: {0}")]
    IndexedReaderError(String),
}

pub struct Genome {
    reader: IndexedReader<std::fs::File>,
}

impl Genome {
    pub fn new_from_file_path(file_path: &Path) -> Result<Self, GenomeError> {
        Ok(Self {
            reader: IndexedReader::from_file(&file_path)
                .map_err(|e| GenomeError::IndexedReaderError(e.to_string()))?,
        })
    }

    pub fn get(&mut self, seq_name: &str, start: u64, stop: u64) -> Option<Vec<u8>> {
        self.reader.fetch(seq_name, start, stop).ok()?;
        let mut seq = Vec::with_capacity((stop - start) as usize);
        self.reader.read(&mut seq).ok()?;
        Some(seq)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_genome_get() {
        let mut genome = Genome::new_from_file_path(&Path::new("tests/data/test.fasta")).unwrap();
        assert_eq!(
            genome.get("chr5", 0, 180).unwrap(),
            concat!(
                "CAAGTTTGCTGGGCTTTCGTCATCCTGTAGACAAGCTTCTTTCTCGGTCA",
                "GGGTAATAACGTGGTGCGTGAACTGTACTTTTACTCACGTATGAAGCGCG",
                "GGAGTCAGGGAAAGTGAAGGAGCGCAAAGCATCTGCCGCCAGAGCACAGC",
                "ATCCGTACAGTAGGTCGCTACGACAGCAAG"
            )
            .as_bytes()
        );
    }
}
