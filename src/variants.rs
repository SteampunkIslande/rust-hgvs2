use std::fmt::Debug;

#[derive(Default, PartialEq)]
pub struct Position {
    pub chrom: String,
    pub chrom_start: u64,
    pub chrom_stop: u64,
    pub is_forward_strand: bool,
}

impl Debug for Position {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&format!(
            "<Position {}[{}:{}]>",
            self.chrom, self.chrom_start, self.chrom_stop
        ))
    }
}

/// Shamelessly stolen from: https://github.com/natir/rust_template/blob/5d58df7b0375b07a1f09cbf8357026b7ffa7d7e9/src/lib.rs#L15C5-L21C2
/// Returns the reverse complement of the input sequence
///
/// assert_eq!(revcomp("CGGTAA".as_bytes()), Vec::from("TTACCG"));
/// assert_eq!(revcomp("cggtaa".as_bytes()), Vec::from("ttaccg"));
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    // Complement the sequence
    seq.iter()
        .rev()
        .map(|c| if c & 2 == 0 { c ^ 21 } else { c ^ 4 })
        .collect()
}

impl Position {
    pub fn new(chrom: String, chrom_start: u64, chrom_end: u64, is_forward_strand: bool) -> Self {
        Self {
            chrom,
            chrom_start,
            chrom_stop: chrom_end,
            is_forward_strand,
        }
    }
}

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_revcomp() {
        assert_eq!(revcomp("CGGTAA".as_bytes()), Vec::from("TTACCG"));
        assert_eq!(revcomp("cggtaa".as_bytes()), Vec::from("ttaccg"));
    }
}
