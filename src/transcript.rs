use std::ops::Not;

use crate::{cdna::CDNACoord, variants::Position};

pub struct CDNAMatch {
    pub cdna_start: u64,
    pub cdna_end: u64,
    pub transcript: Transcript,
    pub tx_position: Position,
}

/// A gene may have multiple transcripts with different combinations of exons.
///     We need both exons and cdna_match as need to know exact exon boundaries to work out flanking
pub struct Transcript {
    pub name: String,
    pub version: Option<i64>,
    pub gene_name: Option<String>,
    pub tx_position: Position,
    pub cds_position: Position,
    pub is_default: bool,
    cdna_match: Vec<CDNAMatch>,
    start_codon_transcript_pos: Option<u64>,
    stop_codon_transcript_pos: Option<u64>,
}

impl Transcript {
    pub fn new(
        name: String,
        version: Option<i64>,
        gene_name: String,
        tx_position: Position,
        cds_position: Position,
        is_default: bool,
        cdna_match: Option<Vec<CDNAMatch>>,
        start_codon_transcript_pos: Option<u64>,
        stop_codon_transcript_pos: Option<u64>,
    ) -> Self {
        // Ordered cdna_match
        let cdna_match = match cdna_match {
            Some(mut cdna_match) => {
                cdna_match.sort_by_key(|k| k.tx_position.chrom_start);
                if tx_position.is_forward_strand.not() {
                    cdna_match.reverse();
                }
                cdna_match
            }
            None => Vec::new(),
        };
        Self {
            name,
            version,
            gene_name: Some(gene_name),
            tx_position,
            cds_position,
            is_default,
            cdna_match: cdna_match,
            start_codon_transcript_pos,
            stop_codon_transcript_pos,
        }
    }

    pub fn full_name(&self) -> String {
        match self.version {
            Some(v) => {
                format!("{}.{}", self.name, v)
            }
            None => self.name.to_string(),
        }
    }

    pub fn is_coding(&self) -> bool {
        self.cds_position.chrom_stop - self.cds_position.chrom_start > 0
    }

    pub fn strand(&self) -> &'static str {
        if self.tx_position.is_forward_strand {
            "+"
        } else {
            "-"
        }
    }

    /// Returns the unsorted cdna_match vector
    pub fn cdna_match(&self) -> &Vec<CDNAMatch> {
        todo!()
    }

    pub fn cdna_to_genomic_coord(&self, coord: CDNACoord) {
        todo!()
    }
}
