use anyhow::{bail, Result};
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
pub(crate) fn kmers<'a>(sequence: &'a [u8], k: usize) -> impl Iterator<Item = usize> + 'a {
    sequence.windows(k).flat_map(word_base4)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_words_to_base4() {
        assert_eq!(word_base4(b"AAAA").unwrap(), 0);
        assert_eq!(word_base4(b"CAAA").unwrap(), 64);
        assert_eq!(word_base4(b"ACGT").unwrap(), 27);
        assert_eq!(word_base4(b"ACGT").unwrap(), 27);
    }
    #[test]
    fn sequence_to_kmers_in_base4() {
        let inner = |seq: &[u8], k: usize| -> Vec<usize> { kmers(seq, k).collect::<Vec<usize>>() };
        assert_eq!(inner(b"ACGT", 1), vec![0, 1, 2, 3]);
        assert_eq!(inner(b"CAAAA", 4), vec![64, 0]);
        assert_eq!(inner(b"ACGT", 2), vec![1, 6, 11]);
        assert_eq!(inner(b"NNNN", 2), vec![]);
        assert_eq!(inner(b"CAAAAN", 4), vec![64, 0]);
        assert_eq!(inner(b"CAAAANACGT", 4), vec![64, 0, 27]);
    }
}
