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

use std::fmt::Display;
use thiserror::Error;

use crate::cdna::CDNACoord;

#[derive(Debug)]
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

#[derive(Error, Debug)]
pub enum HGVSNameError {
    #[error("{0}")]
    InvalidHGVSName(InvalidHGVSName),
}

pub mod hgvs_regex {

    use lazy_regex::regex;
    use lazy_regex::regex::Regex;

    static CDNA_ALLELE_REGEXES: &'static [&lazy_regex::Lazy<Regex>] = &[
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

    static PEP_ALLELE_REGEXES: &'static [&lazy_regex::Lazy<Regex>] = &[
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

    static GENOMIC_ALLELE_REGEXES: &'static [&lazy_regex::Lazy<Regex>] = &[
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

#[derive(Debug, Default)]
pub struct HGVSName {
    name: String,
    prefix: String,
    chrom: String,
    transcript: String,
    gene: String,
    kind: String,
    mutation_type: String,
    start: i64,
    end: i64,
    ref_allele: String,
    ref2_allele: String,
    alt_allele: String,
    cdna_start: CDNACoord,
    cdna_end: CDNACoord,
    pep_extra: String,
}

impl HGVSName {
    pub fn new(name: String) -> Self {
        let mut res = Self {
            name,
            ..Default::default()
        };

        if res.name.len() > 0 {
            res.parse();
        }
        res
    }

    fn parse(&mut self) {}
}
