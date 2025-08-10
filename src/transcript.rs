use crate::{cdna::CDNACoord, variants::Position};

/// A gene may have multiple transcripts with different combinations of exons.
pub struct Transcript {
    pub name: String,
    pub version: Option<i64>,
    pub gene_name: Option<String>,
    pub tx_position: Position,
    pub cds_position: Position,
    pub is_default: bool,
}

impl Transcript {
    pub fn new(
        name: String,
        version: Option<i64>,
        gene_name: String,
        tx_position: Position,
        cds_position: Position,
        is_default: bool,
    ) -> Self {
        Self {
            name,
            version,
            gene_name: Some(gene_name),
            tx_position,
            cds_position,
            is_default,
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

    pub fn cdna_to_genomic_coord(&self, coord: &CDNACoord) -> u64 {
        todo!()
    }
}
