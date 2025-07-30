/*
/// HGVS language currently implemented.
///
/// let HGVS: &'static str = ALLELE
///      | PREFIX_NAME : ALLELE
///
/// let PREFIX_NAME: &'static str = TRANSCRIPT
///             | TRANSCRIPT '(' GENE ')'
///
/// let TRANSCRIPT: &'static str = TRANSCRIPT_NAME
///            | TRANSCRIPT_NAME '.' TRANSCRIPT_VERSION
///
/// let TRANSCRIPT_VERSION: &'static str = NUMBER
///
/// let ALLELE: &'static str = 'c.' CDNA_ALLELE    # cDNA
///        | 'g.' GENOMIC_ALLELE # genomic
/// NC_ALLELE =
/// RNA_ALLELE =
/// let CDNA_ALLELE: &'static str = CDNA_COORD SINGLE_BASE_CHANGE
///             | CDNA_COORD_RANGE MULTI_BASE_CHANGE
///
/// GENOMIC_ALLELE =
/// let MIT_ALLELE: &'static str = COORD SINGLE_BASE_CHANGE
///            | COORD_RANGE MULTI_BASE_CHANGE
///
/// let SINGLE_BASE_CHANGE: &'static str = CDNA_ALLELE = CDNA_COORD BASE '='        # no change
///                    | CDNA_COORD BASE '>' BASE                 # substitution
///                    | CDNA_COORD 'ins' BASE                    # 1bp insertion
///                    | CDNA_COORD 'del' BASE                    # 1bp deletion
///                    | CDNA_COORD 'dup' BASE                    # 1bp duplication
///                    | CDNA_COORD 'ins'                         # 1bp insertion
///                    | CDNA_COORD 'del'                         # 1bp deletion
///                    | CDNA_COORD 'dup'                         # 1bp duplication
///                    | CDNA_COORD 'del' BASE 'ins' BASE         # 1bp indel
///                    | CDNA_COORD 'delins'
///        | 'm.' MIT_ALLELE     # mitochondrial sequence
///        | 'n.' NC_ALLELE      # non-coding RNA reference sequence
///        | 'r.' RNA_ALLELE     # RNA sequence (like r.76a>u)
///        | 'p.' PROTEIN_ALLELE # protein sequence (like  p.Lys76Asn)
///
/// NC_ALLELE =
/// RNA_ALLELE =
/// let CDNA_ALLELE: &'static str = CDNA_COORD SINGLE_BASE_CHANGE
///             | CDNA_COORD_RANGE MULTI_BASE_CHANGE
///
/// GENOMIC_ALLELE =
/// let MIT_ALLELE: &'static str = COORD SINGLE_BASE_CHANGE
///            | COORD_RANGE MULTI_BASE_CHANGE
///
/// let SINGLE_BASE_CHANGE: &'static str = CDNA_ALLELE = CDNA_COORD BASE '='        # no change
///                    | CDNA_COORD BASE '>' BASE                 # substitution
///                    | CDNA_COORD 'ins' BASE                    # 1bp insertion
///                    | CDNA_COORD 'del' BASE                    # 1bp deletion
///                    | CDNA_COORD 'dup' BASE                    # 1bp duplication
///                    | CDNA_COORD 'ins'                         # 1bp insertion
///                    | CDNA_COORD 'del'                         # 1bp deletion
///                    | CDNA_COORD 'dup'                         # 1bp duplication
///                    | CDNA_COORD 'del' BASE 'ins' BASE         # 1bp indel
///                    | CDNA_COORD 'delins' BASE                 # 1bp indel
///
/// let MULTI_BASE_CHANGE: &'static str = COORD_RANGE 'del' BASES             # deletion
///                   | COORD_RANGE 'ins' BASES             # insertion
///                   | COORD_RANGE 'dup' BASES             # duplication
///                   | COORD_RANGE 'del'                   # deletion
///                   | COORD_RANGE 'dup'                   # duplication
///                   | COORD_RANGE 'del' BASES 'ins' BASES # indel
///                   | COORD_RANGE 'delins' BASES          # indel
///
///
/// let AMINO1: &'static str = [GAVLIMFWPSTCYNQDEKRH]
///
/// let AMINO3: &'static str = 'Gly' | 'Ala' | 'Val' | 'Leu' | 'Ile' | 'Met' | 'Phe' | 'Trp' | 'Pro'
///        | 'Ser' | 'Thr' | 'Cys' | 'Tyr' | 'Asn' | 'Gln' | 'Asp' | 'Glu' | 'Lys'
///        | 'Arg' | 'His'
///
/// let PROTEIN_ALLELE: &'static str = AMINO3 COORD '='               # no peptide change
///                | AMINO1 COORD '='               # no peptide change
///                | AMINO3 COORD AMINO3 PEP_EXTRA  # peptide change
///                | AMINO1 COORD AMINO1 PEP_EXTRA  # peptide change
///                | AMINO3 COORD '_' AMINO3 COORD PEP_EXTRA        # indel
///                | AMINO1 COORD '_' AMINO1 COORD PEP_EXTRA        # indel
///                | AMINO3 COORD '_' AMINO3 COORD PEP_EXTRA AMINO3 # indel
///                | AMINO1 COORD '_' AMINO1 COORD PEP_EXTRA AMINO1 # indel
///
/// # A genomic range:
/// let COORD_RANGE: &'static str = COORD '_' COORD
///
/// # A cDNA range:
/// let CDNA_COORD_RANGE: &'static str = CDNA_COORD '_' CDNA_COORD
///
/// # A cDNA coordinate:
/// let CDNA_COORD: &'static str = COORD_PREFIX COORD
///            | COORD_PREFIX COORD OFFSET_PREFIX OFFSET
/// let COORD_PREFIX: &'static str = '' | '-' | '*'
/// let COORD: &'static str = NUMBER
/// let OFFSET_PREFIX: &'static str = '-' | '+'
/// let OFFSET: &'static str = NUMBER
///
/// # Primatives:
/// let NUMBER: &'static str = \d+
/// let BASE: &'static str = [ACGT]
/// let BASES: &'static str = BASE+
*/

use crate::{
    cdna::{CDNACoord, CDNAError, Landmark},
    hgvs_name::refseq_prefixes::get_refseq_type,
    transcript::Transcript,
};
use lazy_regex::lazy_regex;
use std::{
    fmt::{Debug, Display},
    ops::Not,
    str::FromStr,
};

#[derive(Debug, thiserror::Error)]
pub struct InvalidHGVSName {
    name: Option<String>,
    part: String,
    reason: String,
}

impl Display for InvalidHGVSName {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self.name {
            Some(n) => f.write_fmt(format_args!(
                "Invalid HGVS {} {}: {}",
                self.part, n, self.reason
            )),
            None => f.write_fmt(format_args!("Invalid HGVS {}: {}", self.part, self.reason)),
        }
    }
}

impl InvalidHGVSName {
    pub fn with_reason(reason: &str) -> Self {
        Self {
            name: None,
            part: "".to_string(),
            reason: reason.to_string(),
        }
    }
    pub fn new(name: Option<&str>, part: &str, reason: &str) -> Self {
        Self {
            name: name.map(String::from),
            part: part.to_string(),
            reason: reason.to_string(),
        }
    }
}

#[derive(thiserror::Error, Debug)]
pub enum HGVSNameError {
    #[error(transparent)]
    InvalidHGVSNameError(#[from] InvalidHGVSName),

    #[error(transparent)]
    CDNAError(#[from] CDNAError),
}

pub mod hgvs_regex {

    use lazy_regex::regex;
    use lazy_regex::regex::Regex;

    pub static CDNA_ALLELE_REGEXES: &'static [&lazy_regex::Lazy<Regex>] = &[
        regex!(
            r"^(?P<start>(?P<start_coord_prefix>|-|\*)(?P<start_coord>\d+)((?P<start_offset_prefix>-|\+)(?P<start_offset>\d+))?)(?P<mutation_type>=)$"
        ),
        regex!(
            r"^(?P<start>(?P<start_coord_prefix>|-|\*)(?P<start_coord>\d+)((?P<start_offset_prefix>-|\+)(?P<start_offset>\d+))?)(?P<ref>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+)(?P<mutation_type>=)$"
        ),
        regex!(
            r"^(?P<start>(?P<start_coord_prefix>|-|\*)(?P<start_coord>\d+)((?P<start_offset_prefix>-|\+)(?P<start_offset>\d+))?)(?P<ref>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+)(?P<mutation_type>>)(?P<alt>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+)$"
        ),
        regex!(
            r"^(?P<start>(?P<start_coord_prefix>|-|\*)(?P<start_coord>\d+)((?P<start_offset_prefix>-|\+)(?P<start_offset>\d+))?)(?P<mutation_type>ins)(?P<alt>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+)$"
        ),
        regex!(
            r"^(?P<start>(?P<start_coord_prefix>|-|\*)(?P<start_coord>\d+)((?P<start_offset_prefix>-|\+)(?P<start_offset>\d+))?)(?P<mutation_type>del)(?P<ref>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+)$"
        ),
        regex!(
            r"^(?P<start>(?P<start_coord_prefix>|-|\*)(?P<start_coord>\d+)((?P<start_offset_prefix>-|\+)(?P<start_offset>\d+))?)(?P<mutation_type>dup)(?P<ref>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+)$"
        ),
        regex!(
            r"^(?P<start>(?P<start_coord_prefix>|-|\*)(?P<start_coord>\d+)((?P<start_offset_prefix>-|\+)(?P<start_offset>\d+))?)(?P<mutation_type>del)$"
        ),
        regex!(
            r"^(?P<start>(?P<start_coord_prefix>|-|\*)(?P<start_coord>\d+)((?P<start_offset_prefix>-|\+)(?P<start_offset>\d+))?)(?P<mutation_type>dup)$"
        ),
        regex!(
            r"^(?P<start>(?P<start_coord_prefix>|-|\*)(?P<start_coord>\d+)((?P<start_offset_prefix>-|\+)(?P<start_offset>\d+))?)_(?P<end>(?P<end_coord_prefix>|-|\*)(?P<end_coord>\d+)((?P<end_offset_prefix>-|\+)(?P<end_offset>\d+))?)(?P<mutation_type>ins)(?P<alt>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+)$"
        ),
        regex!(
            r"^(?P<start>(?P<start_coord_prefix>|-|\*)(?P<start_coord>\d+)((?P<start_offset_prefix>-|\+)(?P<start_offset>\d+))?)_(?P<end>(?P<end_coord_prefix>|-|\*)(?P<end_coord>\d+)((?P<end_offset_prefix>-|\+)(?P<end_offset>\d+))?)(?P<mutation_type>del)(?P<ref>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+)$"
        ),
        regex!(
            r"^(?P<start>(?P<start_coord_prefix>|-|\*)(?P<start_coord>\d+)((?P<start_offset_prefix>-|\+)(?P<start_offset>\d+))?)_(?P<end>(?P<end_coord_prefix>|-|\*)(?P<end_coord>\d+)((?P<end_offset_prefix>-|\+)(?P<end_offset>\d+))?)(?P<mutation_type>dup)(?P<ref>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+)$"
        ),
        regex!(
            r"^(?P<start>(?P<start_coord_prefix>|-|\*)(?P<start_coord>\d+)((?P<start_offset_prefix>-|\+)(?P<start_offset>\d+))?)_(?P<end>(?P<end_coord_prefix>|-|\*)(?P<end_coord>\d+)((?P<end_offset_prefix>-|\+)(?P<end_offset>\d+))?)(?P<mutation_type>del)$"
        ),
        regex!(
            r"^(?P<start>(?P<start_coord_prefix>|-|\*)(?P<start_coord>\d+)((?P<start_offset_prefix>-|\+)(?P<start_offset>\d+))?)_(?P<end>(?P<end_coord_prefix>|-|\*)(?P<end_coord>\d+)((?P<end_offset_prefix>-|\+)(?P<end_offset>\d+))?)(?P<mutation_type>dup)$"
        ),
        regex!(
            r"^(?P<delins>(?P<start>(?P<start_coord_prefix>|-|\*)(?P<start_coord>\d+)((?P<start_offset_prefix>-|\+)(?P<start_offset>\d+))?)del(?P<ref>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+)ins(?P<alt>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+))$"
        ),
        regex!(
            r"^(?P<delins>(?P<start>(?P<start_coord_prefix>|-|\*)(?P<start_coord>\d+)((?P<start_offset_prefix>-|\+)(?P<start_offset>\d+))?)_(?P<end>(?P<end_coord_prefix>|-|\*)(?P<end_coord>\d+)((?P<end_offset_prefix>-|\+)(?P<end_offset>\d+))?)del(?P<ref>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+)ins(?P<alt>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+))$"
        ),
        regex!(
            r"^(?P<delins>(?P<start>(?P<start_coord_prefix>|-|\*)(?P<start_coord>\d+)((?P<start_offset_prefix>-|\+)(?P<start_offset>\d+))?)delins(?P<alt>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+))$"
        ),
        regex!(
            r"^(?P<delins>(?P<start>(?P<start_coord_prefix>|-|\*)(?P<start_coord>\d+)((?P<start_offset_prefix>-|\+)(?P<start_offset>\d+))?)_(?P<end>(?P<end_coord_prefix>|-|\*)(?P<end_coord>\d+)((?P<end_offset_prefix>-|\+)(?P<end_offset>\d+))?)delins(?P<alt>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+))$"
        ),
    ];

    pub static PEP_ALLELE_REGEXES: &'static [&lazy_regex::Lazy<Regex>] = &[
        regex!(r"^(?P<ref>([A-Z]([a-z]{2}))+)(?P<start>\d+)(?P<extra>(|=|\?)(|fs))$"),
        regex!(
            r"^(?P<ref>([A-Z]([a-z]{2}))+)(?P<start>\d+)(?P<alt>([A-Z]([a-z]{2}))+)(?P<extra>(|=|\?)(|fs))$"
        ),
        regex!(
            r"^(?P<delins>(?P<ref>([A-Z]([a-z]{2}))+)(?P<start>\d+)_(?P<ref2>([A-Z]([a-z]{2}))+)(?P<end>\d+)(?P<extra>(|=|\?)(|fs)))$"
        ),
        regex!(
            r"^(?P<delins>(?P<ref>([A-Z]([a-z]{2}))+)(?P<start>\d+)_(?P<ref2>([A-Z]([a-z]{2}))+)(?P<end>\d+)(?P<alt>([A-Z]([a-z]{2}))+)(?P<extra>(|=|\?)(|fs)))$"
        ),
    ];

    pub static GENOMIC_ALLELE_REGEXES: &'static [&lazy_regex::Lazy<Regex>] = &[
        regex!(r"^(?P<start>\d+)(?P<mutation_type>=)$"),
        regex!(
            r"^(?P<start>\d+)(?P<ref>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+)(?P<mutation_type>=)$"
        ),
        regex!(
            r"^(?P<start>\d+)(?P<ref>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+)(?P<mutation_type>>)(?P<alt>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+)$"
        ),
        regex!(
            r"^(?P<start>\d+)(?P<mutation_type>ins)(?P<alt>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+)$"
        ),
        regex!(
            r"^(?P<start>\d+)(?P<mutation_type>del)(?P<ref>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+)$"
        ),
        regex!(
            r"^(?P<start>\d+)(?P<mutation_type>dup)(?P<ref>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+)$"
        ),
        regex!(r"^(?P<start>\d+)(?P<mutation_type>del)$"),
        regex!(r"^(?P<start>\d+)(?P<mutation_type>dup)$"),
        regex!(r"^(?P<start>\d+)_(?P<end>\d+)(?P<mutation_type>=)$"),
        regex!(
            r"^(?P<start>\d+)_(?P<end>\d+)(?P<mutation_type>ins)(?P<alt>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+)$"
        ),
        regex!(
            r"^(?P<start>\d+)_(?P<end>\d+)(?P<mutation_type>del)(?P<ref>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+)$"
        ),
        regex!(
            r"^(?P<start>\d+)_(?P<end>\d+)(?P<mutation_type>dup)(?P<ref>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+)$"
        ),
        regex!(r"^(?P<start>\d+)_(?P<end>\d+)(?P<mutation_type>del)$"),
        regex!(r"^(?P<start>\d+)_(?P<end>\d+)(?P<mutation_type>dup)$"),
        regex!(
            r"^(?P<delins>(?P<start>\d+)del(?P<ref>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+)ins(?P<alt>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+))$"
        ),
        regex!(
            r"^(?P<delins>(?P<start>\d+)_(?P<end>\d+)del(?P<ref>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+)ins(?P<alt>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+))$"
        ),
        regex!(r"^(?P<delins>(?P<start>\d+)delins(?P<alt>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+))$"),
        regex!(
            r"^(?P<delins>(?P<start>\d+)_(?P<end>\d+)delins(?P<alt>[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+))$"
        ),
    ];
}

pub mod refseq_prefixes {
    use lazy_static::lazy_static;
    use std::collections::HashMap;

    lazy_static! {
        static ref REFSEQ_PREFIX_LOOKUP: HashMap<&'static str, (&'static str, &'static str)> = {
            let mut m = HashMap::new();
            m.insert(
                "AC_",
                (
                    "genomic",
                    "Complete genomic molecule, usually alternate assembly",
                ),
            );
            m.insert(
                "NC_",
                (
                    "genomic",
                    "Complete genomic molecule, usually reference assembly",
                ),
            );
            m.insert("NG_", ("genomic", "Incomplete genomic region"));
            m.insert("NT_", ("genomic", "Contig or scaffold, clone-based or WGS"));
            m.insert("NW_", ("genomic", "Contig or scaffold, primarily WGS"));
            m.insert("NS_", ("genomic", "Environmental sequence"));
            m.insert("NZ_", ("genomic", "Unfinished WGS"));
            m.insert("NM_", ("mRNA", ""));
            m.insert("NR_", ("RNA", ""));
            m.insert("XM_", ("mRNA", "Predicted model"));
            m.insert("XR_", ("RNA", "Predicted model"));
            m.insert("AP_", ("Protein", "Annotated on AC_ alternate assembly"));
            m.insert(
                "NP_",
                ("Protein", "Associated with an NM_ or NC_ accession"),
            );
            m.insert("YP_", ("Protein", ""));
            m.insert(
                "XP_",
                (
                    "Protein",
                    "Predicted model, associated with an XM_ accession",
                ),
            );
            m.insert(
                "ZP_",
                (
                    "Protein",
                    "Predicted model, annotated on NZ_ genomic records",
                ),
            );
            m
        };
    }

    /// Return the RefSeq type for a refseq name.
    pub fn get_refseq_type(name: &str) -> Option<&'static str> {
        REFSEQ_PREFIX_LOOKUP.get(&name[..3]).map(|o| o.0)
    }
}

#[derive(Default)]
pub struct HGVSName {
    name: String,
    prefix: String,
    chrom: String,
    transcript: String,
    gene: String,
    kind: String,
    mutation_type: Option<String>,
    start: i64,
    end: i64,
    ref_allele: String,
    ref2_allele: String,
    alt_allele: String,
    cdna_start: Option<CDNACoord>,
    cdna_end: Option<CDNACoord>,
    pep_extra: String,
}

impl Debug for HGVSName {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("HGVSName({})", self.format(None, None)))
    }
}

impl Display for HGVSName {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.format(None, None))
    }
}

impl HGVSName {
    pub fn new_from_name(name: String) -> Result<Self, HGVSNameError> {
        let mut res = Self {
            name,
            ..Default::default()
        };

        if res.name.len() > 0 {
            res = res.parse()?;
        }
        Ok(res)
    }

    fn parse(self) -> Result<Self, HGVSNameError> {
        // Extract the data we need into owned strings to avoid borrowing issues
        let (prefix, allele) = if let Some((p, a)) = self.name.split_once(':') {
            (p.to_string(), a.to_string())
        } else {
            (String::new(), self.name.clone())
        };

        self.parse_prefix(&prefix).parse_allele(&allele)?.validate()
    }

    /// Parse a HGVS prefix (gene/transcript/chromosome).
    /// Some examples of full hgvs names with transcript include:
    ///     NM_007294.3:c.2207A>C
    ///     NM_007294.3(BRCA1):c.2207A>C
    ///     BRCA1{NM_007294.3}:c.2207A>C
    fn parse_prefix(mut self, prefix: &str) -> Self {
        self.prefix = prefix.to_string();
        if self.prefix.is_empty() {
            return self;
        }

        // Transcript and gene given with parens.
        // example: NM_007294.3(BRCA1):c.2207A>C
        let reg_captures =
            lazy_regex!(r"^(?P<transcript>[^(]+)\((?P<gene>[^)]+)\)$").captures(&self.prefix);
        if let Some(reg_matches) = reg_captures {
            self.transcript = reg_matches
                .name("transcript")
                .map(|o| o.as_str().to_string())
                .unwrap_or_default();
            self.gene = reg_matches
                .name("gene")
                .map(|o| o.as_str().to_string())
                .unwrap_or_default();
            return self;
        }

        // Transcript and gene given with braces.
        // example: BRCA1{NM_007294.3}:c.2207A>C
        let reg_captures =
            lazy_regex!(r"^(?P<gene>[^{]+)\{(?P<transcript>[^}]+)\}$").captures(&self.prefix);
        if let Some(reg_matches) = reg_captures {
            self.transcript = reg_matches
                .name("transcript")
                .map(|o| o.as_str().to_string())
                .unwrap_or_default();
            self.gene = reg_matches
                .name("gene")
                .map(|o| o.as_str().to_string())
                .unwrap_or_default();
            return self;
        }
        // Determine using Ensembl type.
        if prefix.starts_with("ENST") {
            self.transcript = prefix.to_string();
            return self;
        }
        // Determine using LRG type.
        if prefix.starts_with("LRG_") {
            self.transcript = prefix.to_string();
            return self;
        }

        // Determine using refseq type
        if let Some(refseq_type) = get_refseq_type(prefix) {
            match refseq_type {
                "mRNA" | "RNA" => {
                    self.transcript = prefix.to_string();
                    return self;
                }
                "genomic" => {
                    self.chrom = prefix.to_string();
                    return self;
                }
                _ => {}
            }
        }

        // Assume chrom
        if prefix.starts_with("chr") {
            self.chrom = prefix.to_string();
            return self;
        }

        // Assume gene name
        self.gene = prefix.to_string();

        self
    }

    /// Parse a HGVS allele description.
    ///    Some examples include:
    ///      cDNA substitution: c.101A>C,
    ///      cDNA indel: c.3428delCinsTA, c.1000_1003delATG, c.1000_1001insATG
    ///      No protein change: p.Glu1161=
    ///      Protein change: p.Glu1161Ser
    ///      Protein frameshift: p.Glu1161_Ser1164?fs
    ///      Genomic substitution: g.1000100A>T
    ///      Genomic indel: g.1000100_1000102delATG
    fn parse_allele(mut self, allele: &str) -> Result<Self, HGVSNameError> {
        if allele.contains(".").not() {
            return Err(InvalidHGVSName::new(
                Some(allele),
                "allele",
                "expected kind \"c.\",\"p.\",\"g.\", etc",
            )
            .into());
        }
        let (kind, details) = if let Some((kind, details)) = allele.split_once(".") {
            self.kind = kind.to_string();
            self.mutation_type = None;
            (kind, details)
        } else {
            (allele, "")
        };
        match kind {
            "c" | "n" => {
                self = self.parse_cdna(details)?;
                if kind == "n" {
                    if self.cdna_start.as_ref().is_some_and(|o| o.coord < 0) {
                        return Err(InvalidHGVSName::new(
                            Some(allele),
                            "allele",
                            "Non-coding transcript cannot contain negative (5'UTR) coordinates",
                        )
                        .into());
                    }
                    if self
                        .cdna_start
                        .as_ref()
                        .is_some_and(|o| o.landmark == Landmark::CdnaStopCodon)
                        || self
                            .cdna_end
                            .as_ref()
                            .is_some_and(|o| o.landmark == Landmark::CdnaStopCodon)
                    {
                        return Err(InvalidHGVSName::new(
                            Some(allele),
                            "allele",
                            "Non-coding transcript cannot contain '*' (3'UTR) coordinate",
                        )
                        .into());
                    }
                }
            }
            "p" => {
                self = self.parse_protein(&details)?;
            }
            "g" | "m" => {
                self = self.parse_genome(&details)?;
            }
            _ => {
                return Err(
                    InvalidHGVSName::with_reason(&format!("unknown kind: {}", allele)).into(),
                );
            }
        }
        Ok(self)
    }

    ///Parse a HGVS cDNA name.
    ///    Some examples include:
    ///      Substitution: 101A>C,
    ///      Indel: 3428delCinsTA, 1000_1003delATG, 1000_1001insATG
    fn parse_cdna(mut self, details: &str) -> Result<Self, HGVSNameError> {
        for regex in hgvs_regex::CDNA_ALLELE_REGEXES {
            if let Some(groups) = regex.captures(details) {
                // Parse mutation type
                if groups.name("delins").is_some() {
                    self.mutation_type = Some("delins".to_string());
                } else {
                    self.mutation_type =
                        groups.name("mutation_type").map(|o| o.as_str().to_string());
                }

                // Parse coordinates
                let start_val = groups.name("start").map(|o| o.as_str()).ok_or(
                    InvalidHGVSName::with_reason(
                        "Ill-formed regular expression doesn't contain group `start`",
                    ),
                )?;
                self.cdna_start = CDNACoord::from_str(start_val).ok();
                self.cdna_end = if let Some(cdna_end) = groups.name("end").map(|o| o.as_str()) {
                    Some(CDNACoord::from_str(cdna_end)?)
                } else {
                    Some(CDNACoord::from_str(start_val)?)
                };

                // Parse alleles
                self.ref_allele = groups
                    .name("ref")
                    .map(|o| o.as_str().to_string())
                    .unwrap_or_default();
                self.alt_allele = groups
                    .name("alt")
                    .map(|o| o.as_str().to_string())
                    .unwrap_or_default();

                // Convert numerical alleles
                if let Some(rep_count) = self.ref_allele.parse::<i64>().ok() {
                    self.ref_allele = "N".repeat(rep_count as usize);
                }
                if let Some(rep_count) = self.alt_allele.parse::<i64>().ok() {
                    self.alt_allele = "N".repeat(rep_count as usize);
                }

                // Convert duplication alleles
                if self.mutation_type.as_ref().is_some_and(|x| x == "dup") {
                    self.alt_allele = self.ref_allele.repeat(2);
                }

                // Convert no match alleles
                if self.mutation_type.as_ref().is_some_and(|x| x == "=") {
                    // Copy the values (to avoid partial move)
                    self.alt_allele = self.ref_allele.to_string();
                }

                return Ok(self);
            }
        }
        return Err(InvalidHGVSName::new(Some(details), "cDNA allele", "").into());
    }

    /// Parse a HGVS protein name.
    ///    Some examples include:
    ///      No change: Glu1161=
    ///      Change: Glu1161Ser
    ///      Frameshift: Glu1161_Ser1164?fs
    fn parse_protein(mut self, details: &str) -> Result<Self, HGVSNameError> {
        for regex in hgvs_regex::PEP_ALLELE_REGEXES {
            if let Some(groups) = regex.captures(details) {
                //Parse mutation type
                self.mutation_type = Some(
                    if let Some(d) = groups.name("delins").map(|o| o.as_str()) {
                        d
                    } else {
                        ">"
                    }
                    .to_string(),
                );

                // Parse coordinates
                self.start = groups
                    .name("start")
                    .expect("Regex should have captured named group `start`")
                    .as_str()
                    .parse()
                    .expect(
                        "Regex captured named group `start` which should match valid int64 only",
                    );
                self.end = if let Some(end_val) = groups.name("end") {
                    end_val.as_str().parse().expect(
                        "Regex captured named group `end` which should match valid int64 only",
                    )
                } else {
                    self.start
                };

                // Parse alleles
                self.ref_allele = groups
                    .name("ref")
                    .map(|o| o.as_str())
                    .unwrap_or("")
                    .to_string();

                if let Some(ref2_allele) = groups.name("ref2").map(|o| o.as_str()) {
                    self.ref2_allele = ref2_allele.to_string();
                    self.alt_allele = groups
                        .name("alt")
                        .map(|o| o.as_str().to_string())
                        .unwrap_or_default();
                } else {
                    // If alt is not given, assume matching with ref
                    self.ref2_allele = self.ref_allele.to_string();
                    self.alt_allele = groups
                        .name("alt")
                        .map(|o| o.as_str())
                        .unwrap_or(&self.ref_allele)
                        .to_string();
                }
                self.pep_extra = groups
                    .name("extra")
                    .map(|o| o.as_str())
                    .map(String::from)
                    .unwrap_or_default();
                return Ok(self);
            }
        }
        return Err(InvalidHGVSName::new(Some(details), "protein allele", "").into());
    }

    ///Parse a HGVS genomic name.
    ///    Som examples include:
    ///      Substitution: 1000100A>T
    ///      Indel: 1000100_1000102delATG
    fn parse_genome(mut self, details: &str) -> Result<Self, HGVSNameError> {
        for regex in hgvs_regex::GENOMIC_ALLELE_REGEXES {
            if let Some(groups) = regex.captures(details) {
                //Parse mutation type
                self.mutation_type = Some(
                    if groups.name("delins").map(|o| o.as_str()).is_some() {
                        "delins"
                    } else {
                        groups
                            .name("mutation_type")
                            .expect("Invalid regex, did not capture named group `mutation_type`")
                            .as_str()
                    }
                    .to_string(),
                );

                //Parse coordinates
                self.start = groups
                    .name("start")
                    .expect("Invalid regex, did not capture named group `start`")
                    .as_str()
                    .parse()
                    .expect("Regex captured named group `start` which should be a valid i64");

                self.end = if let Some(end_val) = groups.name("end").map(|o| o.as_str()) {
                    end_val
                        .parse()
                        .expect("Regex captured named group `start` which should be a valid i64")
                } else {
                    self.start
                };

                // Parse alleles
                self.ref_allele = groups
                    .name("ref")
                    .map(|o| o.as_str().to_string())
                    .unwrap_or_default();
                self.alt_allele = groups
                    .name("alt")
                    .map(|o| o.as_str().to_string())
                    .unwrap_or_default();

                // Convert numerical alleles
                if let Some(rep_count) = self.ref_allele.parse::<i64>().ok() {
                    self.ref_allele = "N".repeat(rep_count as usize);
                }
                if let Some(rep_count) = self.alt_allele.parse::<i64>().ok() {
                    self.alt_allele = "N".repeat(rep_count as usize);
                }

                // Convert duplication alleles
                if self.mutation_type.as_ref().is_some_and(|x| x == "dup") {
                    self.alt_allele = self.ref_allele.repeat(2);
                }

                // Convert no match alleles
                if self.mutation_type.as_ref().is_some_and(|x| x == "=") {
                    // Copy the values (to avoid partial move)
                    self.alt_allele = self.ref_allele.to_string();
                }
            }
        }
        return Err(InvalidHGVSName::new(Some(details), "genomic allele", "").into());
    }

    fn validate(self) -> Result<Self, HGVSNameError> {
        if self.start > self.end {
            Err(InvalidHGVSName::with_reason("Coordinates are nonincreasing").into())
        } else {
            Ok(self)
        }
    }

    /// Generate a HGVS name as a string.
    fn format(&self, use_prefix: Option<bool>, use_gene: Option<bool>) -> String {
        let use_prefix = use_prefix.unwrap_or(false);
        let use_gene = use_gene.unwrap_or(false);

        let allele = match self.kind.as_str() {
            "c" | "n" => format!("{}.{}", self.kind, self.format_cdna()),
            "p" => format!("{}.{}", self.kind, self.format_protein()),
            "g" | "m" => format!("{}.{}", self.kind, self.format_genome()),
            _ => "".to_string(), //TODO: Maybe should return an error
        };
        let prefix = if use_prefix {
            self.format_prefix(use_gene)
        } else {
            "".to_string()
        };
        if prefix.is_empty().not() {
            format!("{}:{}", prefix, allele)
        } else {
            allele
        }
    }

    /// Generate HGVS trancript/gene prefix.
    ///     Some examples of full hgvs names with transcript include:
    ///       NM_007294.3:c.2207A>C
    ///       NM_007294.3(BRCA1):c.2207A>C
    fn format_prefix(&self, use_gene: bool) -> String {
        match self.kind.as_str() {
            "g" | "m" => {
                if self.chrom.is_empty().not() {
                    return self.chrom.to_string();
                }
            }
            _ => {}
        }
        if self.transcript.is_empty().not() {
            if use_gene && self.gene.is_empty().not() {
                return format!("{}({})", &self.transcript, &self.gene);
            } else {
                return self.transcript.to_string();
            }
        } else {
            if use_gene {
                return self.gene.to_string();
            } else {
                return "".to_string();
            }
        }
    }

    /// Generate HGVS cDNA coordinates string.
    fn format_cdna_coords(&self) -> String {
        match (&self.cdna_start, &self.cdna_end) {
            (Some(cdna_start), Some(cdna_end)) => {
                if cdna_start == cdna_end {
                    cdna_start.to_string()
                } else {
                    format!("{}_{}", cdna_start, cdna_end)
                }
            }
            _ => "".to_string(),
        }
    }

    /// Generate HGVS DNA allele.
    fn format_dna_allele(&self) -> String {
        if let Some(ref mut_type) = self.mutation_type {
            match mut_type.as_str() {
                "=" => {
                    format!("{}=", self.ref_allele)
                }
                ">" => {
                    format!("{}>{}", self.ref_allele, self.alt_allele)
                }
                "delins" => {
                    format!("del{}ins{}", self.ref_allele, self.alt_allele)
                }
                "del" | "dup" => {
                    format!("{}{}", mut_type, self.ref_allele)
                }
                "ins" => {
                    format!("{}{}", mut_type, self.alt_allele)
                }
                "inv" => mut_type.to_string(),
                _ => "".to_string(), //TODO: Should raise!
            }
        } else {
            //TODO: Should raise too!
            "".to_string()
        }
    }

    /// Generate HGVS cDNA allele.
    ///     Some examples include:
    ///       Substitution: 101A>C,
    ///       Indel: 3428delCinsTA, 1000_1003delATG, 1000_1001insATG
    fn format_cdna(&self) -> String {
        format!("{}{}", self.format_cdna_coords(), self.format_dna_allele())
    }

    /// Generate HGVS protein name.
    ///     Some examples include:
    ///       No change: Glu1161=
    ///       Change: Glu1161Ser
    ///       Frameshift: Glu1161_Ser1164?fs
    fn format_protein(&self) -> String {
        if self.start == self.end {
            if self.ref_allele == self.ref2_allele && self.ref_allele == self.alt_allele {
                return format!(
                    "{}{}{}",
                    self.ref_allele,
                    self.start,
                    if self.pep_extra.is_empty().not() {
                        &self.pep_extra
                    } else {
                        "="
                    }
                );
            }
            if self.ref_allele == self.ref2_allele && self.ref_allele != self.alt_allele {
                return format!(
                    "{}{}{}{}",
                    self.ref_allele, self.start, self.alt_allele, self.pep_extra
                );
            }
        } else {
            return format!(
                "{}{}_{}{}{}",
                self.ref_allele, self.start, self.ref2_allele, self.end, self.pep_extra
            );
        }
        "".to_string()
    }

    /// Generate HGVS cDNA coordinates string.
    fn format_coords(&self) -> String {
        if self.start == self.end {
            self.start.to_string()
        } else {
            format!("{}_{}", self.start, self.end)
        }
    }

    /// Generate HGVS genomic allele.
    ///     Some examples include:
    ///       Substitution: 1000100A>T
    ///       Indel: 1000100_1000102delATG
    fn format_genome(&self) -> String {
        format!("{}{}", self.format_coords(), self.format_dna_allele())
    }

    /// Return genomic coordinates
    fn get_raw_coords(&self, transcript: Option<&Transcript>) -> Option<(String, i64, i64)> {
        let transcript = transcript?;
        match self.kind.as_ref() {
            "c" | "n" => {
                let chrom = transcript.tx_position.chrom.to_string();
                // let start = t
                None
            }
            _ => None,
        }
    }

    fn get_ref_coords(&self, transcript: Option<&Transcript>) -> Option<(String, i64, i64)> {
        todo!()
    }

    fn get_vcf_coords(&self, transcript: Option<&Transcript>) -> Option<(String, i64, i64)> {
        todo!()
    }
}
