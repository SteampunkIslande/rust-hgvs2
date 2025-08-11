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
pub mod mock_genome {

    use std::collections::HashMap;
    use std::fs::File;
    use std::io::{BufRead, BufReader, Cursor, Read};
    use std::path::Path;

    use rle_vec::RleVec;
    /// Génome mock utilisant une structure d'intervalles pour les requêtes efficaces
    pub struct MockGenome {
        genome: HashMap<String, RleVec<u8>>,
    }

    impl MockGenome {
        pub fn new_from_file_path(file_path: &Path) -> Self {
            let genome = File::open(file_path).map(MockGenome::from_reader).unwrap();
            Self { genome }
        }

        pub fn new_from_static_str(content: &'static str) -> Self {
            let genome = MockGenome::from_reader(Cursor::new(content));
            Self { genome }
        }

        pub fn from_reader(reader: impl Read) -> HashMap<String, RleVec<u8>> {
            let mut genome: HashMap<String, RleVec<u8>> = HashMap::new();
            let reader = BufReader::new(reader);

            for line in reader.lines() {
                if let Ok(line) = line {
                    let parts: Vec<&str> = line.trim().split('\t').collect();
                    let chrom = parts[0];
                    let start: usize = parts[1].parse().unwrap();
                    let end: usize = parts[2].parse().unwrap();
                    let seq = parts[3];
                    let mock_chr = genome.entry(chrom.to_string()).or_insert(RleVec::new());

                    mock_chr.push_n(end, b'N');
                    for (i, c) in seq.chars().map(|e| e as u8).enumerate() {
                        mock_chr.insert(start + i, c);
                    }
                }
            }
            genome
        }

        pub fn get(&mut self, seq_name: &str, start: u64, stop: u64) -> Option<Vec<u8>> {
            let v = self.genome.get(seq_name)?;
            let mut res: Vec<u8> = Vec::with_capacity((stop - start) as usize);
            for i in start..stop {
                res.push(v[i as usize]);
            }

            Some(res)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::mock_genome::MockGenome;
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

    #[test]
    fn test_genome_get_mock() {
        let mut genome = MockGenome::new_from_file_path(&Path::new("tests/data/test_hgvs.genome"));
        assert_eq!(
            &genome.get("chr1", 5933334, 5933375).unwrap(),
            b"CGCTGGACTTCCAAGGTGACACGGCGTCCATGCCCTTCTCG"
        );
        assert_eq!(
            &genome.get("chr1", 5933334, 5933374).unwrap(),
            b"CGCTGGACTTCCAAGGTGACACGGCGTCCATGCCCTTCTC"
        );
        assert_eq!(&genome.get("chr1", 2337999, 2338000).unwrap(), b"C");
    }

    #[test]
    fn test_genome_get_mock_from_static_str() {
        let mut genome = MockGenome::new_from_static_str(
            "chr1\t5927848\t5927889\tGAGCTCCGGGTGATAGAAGCGGAAGACCTGGTCCACCACGT",
        );
        assert_eq!(
            &genome.get("chr1", 5927848, 5927889).unwrap(),
            b"GAGCTCCGGGTGATAGAAGCGGAAGACCTGGTCCACCACGT"
        );
    }
}
