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

pub(crate) fn detect_kmers_across_sequences<'a, T: ExactSizeIterator<Item = &'a str>>(
    sequences: T,
    kmer_size: u32,
) -> Array2<u8> {
    let nrows = 4_usize.pow(kmer_size);
    let ncols = sequences.len();
    let mut kmer_counts = Array2::<u8>::from_elem((nrows, ncols), 0);
    kmer_counts
        .axis_iter_mut(Axis(1))
        .zip(sequences)
        .for_each(|(mut col, seq)| {
            for kmer in kmers(seq.as_bytes(), kmer_size as usize) {
                col[kmer] = 1;
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
    #[test]
    fn kmer_matrix() {
        let sequences = vec!["ACGT", "ACGT", "ACGT"];
        let input = sequences.iter().map(|&s| s);
        let kmer_counts = detect_kmers_across_sequences(input, 1);
        assert_eq!(kmer_counts.shape(), [4, 3]);
        assert_eq!(kmer_counts.column(0).iter().sum::<u8>(), 4);
        assert_eq!(kmer_counts.column(1).iter().sum::<u8>(), 4);
        assert_eq!(kmer_counts.column(2).iter().sum::<u8>(), 4);

        let sequences = vec!["ACGT", "ACGT", "ACG"];
        let input = sequences.iter().map(|&s| s);
        let kmer_counts = detect_kmers_across_sequences(input, 4);
        assert_eq!(kmer_counts.column(0).iter().sum::<u8>(), 1);
        assert_eq!(kmer_counts.column(1).iter().sum::<u8>(), 1);
        assert_eq!(kmer_counts.column(2).iter().sum::<u8>(), 0);
        let col_one = kmer_counts.column(0);
        assert_eq!(
            *col_one.iter().nth(word_base4(b"ACGT").unwrap()).unwrap(),
            1
        );
        let col_two = kmer_counts.column(1);
        assert_eq!(
            *col_two.iter().nth(word_base4(b"ACGT").unwrap()).unwrap(),
            1
        );
    }
}
