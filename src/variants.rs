use std::{fmt::Debug, ops::Not};

use crate::genome::Genome;

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

pub enum IndelJustify {
    Left,
    Right,
}

/// Shamelessly stolen from: https://github.com/natir/rust_template/blob/5d58df7b0375b07a1f09cbf8357026b7ffa7d7e9/src/lib.rs#L15C5-L21C2
/// Returns the reverse complement of the input sequence
///
/// ```rust
/// use pyhgvs2::variants::revcomp;
/// assert_eq!(revcomp("CGGTAA".as_bytes()), Vec::from("TTACCG"));
/// assert_eq!(revcomp("cggtaa".as_bytes()), Vec::from("ttaccg"));
/// ```
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

/// Return a sequence for the genomic region.
///    Coordinates are 0-based, end-exclusive.
pub fn get_sequence(
    genome: &mut Genome,
    chrom: &str,
    start: u64,
    end: u64,
    is_forward_strand: bool,
) -> Option<Vec<u8>> {
    if start > end {
        return None;
    } else {
        let seq = genome.get(&chrom, start, end)?;
        if is_forward_strand.not() {
            return Some(revcomp(seq.as_slice()));
        } else {
            return Some(seq);
        }
    }
}

///Return a sequence for the genomic region
/// Position is 0-based, end-exclusive.
pub fn get_sequence_from_position(genome: &mut Genome, position: &Position) -> Option<Vec<u8>> {
    return get_sequence(
        genome,
        &position.chrom,
        position.chrom_start,
        position.chrom_stop,
        position.is_forward_strand,
    );
}

/// Justify an indel to the left or right along a sequence 'seq'.
///     start, end: 0-based, end-exclusive coordinates of 'indel' within the
///         sequence 'seq'. Inserts denote the insertion point using start=end
///         and deletions indicate the deleted region with (start,end).
///     indel: indel sequence, can be insertion or deletion.
///     seq: a larger sequence containing the indel. Can be a fragment from the
///         genome.
///     justify: Which direction to justify the indel ('left', 'right').
pub fn justify_indel(
    mut start: i64,
    mut end: i64,
    indel: &[u8],
    seq: &[u8],
    justify: IndelJustify,
) -> Option<(i64, i64, Vec<u8>)> {
    let mut indel = indel.to_vec();
    if indel.len() == 0 {
        Some((start, end, indel.to_vec()))
    } else {
        let (start, end, indel) = match justify {
            IndelJustify::Left => {
                while start > 0 && *seq.get((start - 1) as usize)? == indel[indel.len() - 1] {
                    let idx = indel.len() - 1 as usize;
                    indel[idx] = *seq.get((start - 1) as usize)?;
                    start -= 1;
                    end -= 1;
                }
                (start, end, indel)
            }
            IndelJustify::Right => {
                while end < seq.len() as i64 && *seq.get(end as usize)? == indel[0] {
                    indel.remove(0);
                    indel.push(seq[end as usize]);
                    start += 1;
                    end += 1;
                }
                (start, end, indel)
            }
        };
        Some((start, end, indel))
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
