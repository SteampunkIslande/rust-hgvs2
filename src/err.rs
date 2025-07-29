use crate::cdna::CDNAError;
use crate::hgvs_name::HGVSNameError;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum DataStoreError {
    #[error("CDNAError: {0}")]
    CDNAError(#[from] CDNAError),

    #[error("HGVSNameError: {0}")]
    HGVSNameError(#[from] HGVSNameError),
}
