use std::collections::HashMap;
use std::ops::Index;
use std::ops::Range;

#[derive(Debug, Default)]
struct Genome {}

#[derive(Debug, Default)]
struct ChromosomeSubset {
    name: String,
    genome: GenomeSubset,
}

#[derive(Debug, Default)]
struct GenomeSubset {
    genome: Genome,
    chrom: String,
    start: i64,
    end: i64,
    seqid: String,
    _chroms: HashMap<String, ChromosomeSubset>,
}

impl Genome {}

impl Index<&str> for Genome {
    type Output = Vec<u8>;

    fn index(&self, _index: &str) -> &Self::Output {
        todo!("Implement retrieval of chromosome sequence by name")
    }
}

impl ChromosomeSubset {
    pub fn new(name: String, genome: GenomeSubset) -> Self {
        Self { name, genome }
    }
}

impl Index<Range<i64>> for ChromosomeSubset {
    type Output = [u8];

    fn index(&self, index: Range<i64>) -> &Self::Output {
        let (mut start, mut end) = (index.start, index.end);
        start -= self.genome.start;
        end -= self.genome.start;
        &self.genome.genome[&self.genome.seqid][start as usize..end as usize]
    }
}

impl GenomeSubset {
    pub fn new(genome: Genome, chrom: String, start: i64, end: i64, seqid: String) -> Self {
        Self {
            genome,
            chrom,
            start,
            end,
            seqid,
            ..Default::default()
        }
    }
}

impl Index<String> for GenomeSubset {
    type Output = ChromosomeSubset;

    fn index(&self, chrom: String) -> &Self::Output {
        todo!()
    }
}
