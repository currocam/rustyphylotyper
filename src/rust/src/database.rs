use anyhow::{ensure, Result};
use extendr_api::prelude::*;

use crate::kmers;
use std::collections::HashMap;
pub(crate) struct KmerDatabase<T> {
    pub conditional_probs: T,
    pub genera: Vec<String>,
}

impl<'a> KmerDatabase<Array2<f64>> {
    pub fn build<T: ExactSizeIterator<Item = &'a str>, U: ExactSizeIterator<Item = &'a str>>(
        sequences: T,
        genera: U,
        kmer_size: u32,
    ) -> Result<Self> {
        ensure!(sequences.len() == genera.len());
        ensure!(kmer_size > 0);
        // First, deal with the genera and get bijection genus <-> index, and count how many of each
        let genera = genera.collect::<Vec<&str>>();
        // We don't know how many unique, but less than the number of genera
        let mut genera_count: Vec<usize> = Vec::with_capacity(genera.len());
        let mut unique_genera: Vec<String> = Vec::with_capacity(genera.len());
        let mut dictionary: HashMap<&'a str, usize> = HashMap::new();
        let mut index = 0;
        for genus in genera.iter() {
            let x = dictionary.entry(genus).or_insert_with(|| {
                let current_index = index;
                genera_count.push(1); // pseudocount
                unique_genera.push(genus.to_string());
                index += 1;
                current_index
            });
            genera_count[*x] += 1;
        }
        unique_genera.shrink_to_fit();
        let n_unique_genera = index;
        let genera_count: Array1<f64> = genera_count.iter().map(|&count| count as f64).collect();

        // Allocate a matrix to store the counts of each kmer for each genus
        let mut genus_count = Array2::<f64>::zeros((4_usize.pow(kmer_size), n_unique_genera));
        // Allocate prior vector
        let mut kmers_count = Array1::<f64>::zeros(4_usize.pow(kmer_size));
        let n_sequences = sequences.len();
        for (i, sequence) in sequences.enumerate() {
            let genus = dictionary
                .get(genera[i])
                .ok_or_else(|| anyhow::anyhow!("Genus not found"))?;
            for kmer in kmers::kmers(sequence.as_bytes(), kmer_size as usize) {
                // Add one to the count of the kmer for the genus
                genus_count[[kmer, *genus]] += 1.0;
                // Add one to the count of the kmer for the prior
                kmers_count[kmer] += 1.0;
            }
        }
        // Compute prior probabilities
        let priors = (kmers_count + 0.5) / (n_sequences + 1) as f64;
        // Compute Genus-specific conditional probabilities
        let mut conditional_probs = genus_count; // m(w_i)
        for mut col in conditional_probs.axis_iter_mut(Axis(1)) {
            col += &priors; // m(w_i) + P_i
        }
        for mut row in conditional_probs.axis_iter_mut(Axis(0)) {
            row /= &genera_count;
        }
        Ok(KmerDatabase {
            conditional_probs,
            genera: unique_genera,
        })
    }
}

impl KmerDatabase<ArrayView2<'_, f64>> {
    pub fn read_from_list(list: List) -> Result<Self> {
        let names_and_values = list.into_hashmap();
        let conditional_probs = names_and_values
            .get("conditional_prob")
            .ok_or_else(|| anyhow::anyhow!("Missing conditional_prob"))?;
        let conditional_probs: ArrayView2<f64> = conditional_probs
            .try_into()
            .expect("Failed to convert conditional_prob to matrix");
        let genera = names_and_values
            .get("genera")
            .ok_or_else(|| anyhow::anyhow!("Missing genera"))?
            .as_string_vector()
            .ok_or_else(|| anyhow::anyhow!("Failed to convert genera to vector"))?;
        Ok(KmerDatabase {
            conditional_probs,
            genera,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    macro_rules! assert_eq_arrays {
        ($expected:expr, $actual:expr) => {
            for (expected, actual) in $expected.iter().zip($actual.iter()) {
                // Assert that the values are equal
                assert!((*expected - *actual).abs() < 1e-6);
            }
        };
    }
    #[test]
    fn test_genus_specific_conditional_probabilities() {
        // c("03212130", "03212131", "03212131")
        let sequences = vec!["ATGCGCTA", "ATGCGCTC", "ATGCGCTC"];
        let genera = vec!["Genus1", "Genus2", "Genus2"];
        let kmer_size = 3;
        let result = KmerDatabase::build(
            sequences.iter().map(|&s| s),
            genera.iter().map(|&g| g),
            kmer_size,
        )
        .expect("Failed to build kmer database");
        let p = result.conditional_probs;
        // all 3 = (c(1, 2)+0.5) / (c(1, 2) + 1) = 0.74 & 0.8333333
        assert_eq_arrays!(
            vec![0.9375000, 0.9583333],
            p.index_axis(Axis(0), 25).iter().collect::<Vec<&f64>>()
        );
        // only 1 =(c(1, 0)+0.375) / (c(1, 2) + 1)
        assert_eq_arrays!(
            vec![0.6875, 0.1250],
            p.index_axis(Axis(0), 28).iter().collect::<Vec<&f64>>()
        );
        assert_eq_arrays!(
            vec![0.3125, 0.8750],
            p.index_axis(Axis(0), 29).iter().collect::<Vec<&f64>>()
        );
        assert_eq_arrays!(
            vec![0.06250000, 0.04166667],
            p.index_axis(Axis(0), 63).iter().collect::<Vec<&f64>>()
        );
    }
    #[test]
    fn test_unique_genera() {
        // c("03212130", "03212131", "03212131")
        let sequences = vec!["ATGCGCTA", "ATGCGCTC", "ATGCGCTC", ""];
        let genera = vec!["Genus1", "Genus2", "Genus2", "Genus3"];
        let kmer_size = 1;
        let result = KmerDatabase::build(
            sequences.iter().map(|&s| s),
            genera.iter().map(|&g| g),
            kmer_size,
        )
        .expect("Failed to build kmer database");
        let genera = result.genera;
        assert_eq!(genera.len(), 3);
        dbg!(&genera);
        assert_eq!(genera[0], "Genus1");
        assert_eq!(genera[1], "Genus2");
        assert_eq!(genera[2], "Genus3");
    }
}
