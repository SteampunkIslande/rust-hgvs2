use derive_more::Display;
use lazy_regex::regex;
use std::fmt::{Debug, Display};
use std::num::ParseIntError;
use std::str::FromStr;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum CDNAError {
    #[error("Unknown coordinates format: `{0}`")]
    CoordinatesFormatError(String),

    #[error("Invalid number: {0}")]
    InvalidNumber(#[from] ParseIntError),

    #[error("Unknown offset prefix: `{0}`")]
    UnknownOffsetPrefix(String),

    #[error("Unknown coords prefix: `{0}`")]
    UnknownCoordPrefix(String),
}

/// A HGVS cDNA-based coordinate.
///
/// A cDNA coordinate can take one of these forms:
///
/// N = nucleotide N in protein coding sequence (e.g. 11A>G)
///
/// -N = nucleotide N 5' of the ATG translation initiation codon (e.g. -4A>G)
///         NOTE: so located in the 5'UTR or 5' of the transcription initiation
///         site (upstream of the gene, incl. promoter)
///
/// *N = nucleotide N 3' of the translation stop codon (e.g. *6A>G)
///         NOTE: so located in the 3'UTR or 3' of the polyA-addition site
///         (including downstream of the gene)
///
/// N+M = nucleotide M in the intron after (3' of) position N in the coding DNA
///         reference sequence (e.g. 30+4A>G)
///
/// N-M = nucleotide M in the intron before (5' of) position N in the coding
///         DNA reference sequence (e.g. 301-2A>G)
///
/// -N+M / -N-M = nucleotide in an intron in the 5'UTR (e.g. -45+4A>G)
///
/// *N+M / *N-M = nucleotide in an intron in the 3'UTR (e.g. *212-2A>G)
#[derive(Default, PartialEq)]
pub struct CDNACoord {
    pub coord: i64,
    pub offset: i64,
    pub landmark: Landmark,
}

#[derive(Default, Debug, PartialEq, Display)]
pub enum Landmark {
    #[display("cdna_start_codon")]
    #[default]
    CdnaStartCodon,
    #[display("cdna_stop_codon")]
    CdnaStopCodon,
}

impl CDNACoord {
    pub fn new(coord: Option<i64>, offset: Option<i64>, landmark: Option<Landmark>) -> Self {
        Self {
            coord: coord.unwrap_or(0),
            offset: offset.unwrap_or(0),
            landmark: landmark.unwrap_or_default(),
        }
    }
}

impl FromStr for CDNACoord {
    type Err = CDNAError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut res = Self::default();
        if let Some(regex_match) = regex!(r#"(|-|\*)(\d+)((-|\+)(\d+))?"#).captures(s) {
            let (_, [coord_prefix, coord, _, offset_prefix, offset]) = regex_match.extract();
            res.coord = coord.parse()?;
            res.offset = offset.parse().unwrap_or(0);

            match offset_prefix {
                "-" => res.offset *= -1,
                "+" | "" => {}
                _ => return Err(CDNAError::UnknownOffsetPrefix(offset_prefix.into())),
            }

            match coord_prefix {
                "" => res.landmark = Landmark::CdnaStartCodon,
                "-" => {
                    res.coord *= -1;
                    res.landmark = Landmark::CdnaStartCodon
                }
                "*" => res.landmark = Landmark::CdnaStopCodon,
                _ => {
                    return Err(CDNAError::UnknownCoordPrefix(coord_prefix.into()));
                }
            }

            return Ok(res);
        }
        return Err(CDNAError::CoordinatesFormatError(s.into()));
    }
}

impl Display for CDNACoord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let coord_prefix = if self.landmark == Landmark::CdnaStopCodon {
            "*"
        } else {
            ""
        };
        let offset = if self.offset < 0 {
            &format!("{}", self.offset)
        } else if self.offset > 0 {
            &format!("+{}", self.offset)
        } else {
            ""
        };
        f.write_fmt(format_args!("{}{}{}", coord_prefix, self.coord, offset))
    }
}

impl Debug for CDNACoord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.landmark != Landmark::CdnaStartCodon {
            f.write_fmt(format_args!(
                "CDNACoord({},{},'{}')",
                self.coord, self.offset, self.landmark
            ))
        } else {
            f.write_fmt(format_args!("CDNACoord({},{})", self.coord, self.offset))
        }
    }
}
