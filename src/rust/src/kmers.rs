use anyhow::{bail, Result};
use extendr_api::prelude::*;

fn word_base4(seq: &[u8]) -> Result<usize> {
    let mut acc = 0;
    for &base in seq {
        acc = acc * 4
            + match base {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => bail!("Invalid base"),
            }
    }
    Ok(acc)
}

pub(crate) fn detect_kmers_across_sequences<'a, T: ExactSizeIterator<Item = &'a [u8]>>(
    sequences: T,
    kmer_size: u32,
) -> Array2<bool> {
    let nrows = 4_usize.pow(kmer_size);
    let ncols = sequences.len();
    let mut kmer_counts = Array2::<bool>::from_elem((nrows, ncols), false);
    kmer_counts
        .axis_iter_mut(Axis(1))
        .zip(sequences)
        .for_each(|(mut col, seq)| {
            for kmer in kmers(seq, kmer_size as usize) {
                col[kmer] = true;
            }
        });
    kmer_counts
}

pub(crate) fn kmers(sequence: &[u8], k: usize) -> impl Iterator<Item = usize> + '_ {
    sequence.windows(k).flat_map(word_base4)
}

#[cfg(test)]
mod tests {
    use super::*;
    macro_rules! assert_kmer {
        ($seq:expr, $k:expr, $expected:expr) => {
            let result = kmers($seq, $k).collect::<Vec<usize>>();
            assert_eq!(result, $expected);
        };
    }
    #[test]
    fn sequence_to_kmers_in_base4() {
        assert!(word_base4(b"ACCTggC").is_ok());
        assert!(word_base4(b"NNA").is_err());
        assert_kmer!(b"ACGT", 1, [0, 1, 2, 3]);
        assert_kmer!(b"AAAA", 4, [0]);
        assert_kmer!(b"CAAAA", 4, [64, 0]);
        assert_kmer!(b"ACGT", 2, [1, 6, 11]);
        assert_kmer!(b"NNNN", 2, []);
        assert_kmer!(b"", 1, []);
        assert_kmer!(b"ACGT", 100, []);
        assert_kmer!(b"CAAAAN", 4, [64, 0]);
        assert_kmer!(b"CAAAANACGT", 4, [64, 0, 27]);
    }
}
