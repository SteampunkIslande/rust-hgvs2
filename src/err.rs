use crate::cdna::CDNAError;
use crate::genome::GenomeError;
use crate::hgvs_name::HGVSNameError;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum HGVSError {
    #[error("CDNAError: {0}")]
    CDNAError(#[from] CDNAError),

    #[error("HGVSNameError: {0}")]
    HGVSNameError(#[from] HGVSNameError),

    #[error(transparent)]
    GenomeError(#[from] GenomeError),
}
