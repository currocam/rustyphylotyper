use anyhow::{bail, ensure, Result};
use extendr_api::prelude::*;

use crate::kmers;
use std::collections::HashMap;
pub(crate) struct KmerDatabase {
    pub conditional_probs: Array2<f64>,
}

impl KmerDatabase {
    pub fn build<'a, T: ExactSizeIterator<Item = &'a str>, U: ExactSizeIterator<Item = &'a str>>(
        sequences: T,
        genera: U,
        kmer_size: u32,
    ) -> Result<Self> {
        ensure!(sequences.len() == genera.len());
        // First, deal with the genera and get bijection genus <-> index, and count how many of each
        let genera = genera.collect::<Vec<&str>>();
        // We don't know how many unique, but less than the number of genera
        let mut genera_count: Vec<usize> = Vec::with_capacity(genera.len());
        let mut dictionary: HashMap<&'a str, usize> = HashMap::new();
        let mut index = 0;
        for genus in genera.iter() {
            let x = dictionary.entry(genus).or_insert_with(|| {
                let current_index = index;
                genera_count.push(1); // pseudocount
                index += 1;
                current_index
            });
            genera_count[*x] += 1;
        }
        let genera_count: Array1<f64> = genera_count.iter().map(|&count| count as f64).collect();

        // Allocate a matrix to store the counts of each kmer for each genus
        let mut genus_count = Array2::<f64>::zeros((4_usize.pow(kmer_size), index));
        // Allocate prior vector
        let mut kmers_count = Array1::<f64>::zeros(4_usize.pow(kmer_size));
        let n_sequences = sequences.len();
        for (i, sequence) in sequences.enumerate() {
            let genus = dictionary.get(genera[i]).expect("Genus not found");
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
                                                 // dbg1 row 25
        for mut col in conditional_probs.axis_iter_mut(Axis(1)) {
            col += &priors; // m(w_i) + P_i
        }
        for mut row in conditional_probs.axis_iter_mut(Axis(0)) {
            row /= &genera_count;
        }
        let _genera = genera
            .into_iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>();
        Ok(KmerDatabase { conditional_probs })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
        for (expected, actual) in vec![0.9375000, 0.9583333]
            .iter()
            .zip(p.index_axis(Axis(0), 25).iter())
        {
            // Assert that the values are equal
            assert!((dbg!(expected) - dbg!(actual)).abs() < 1e-6);
        }
        // only 1 =(c(1, 0)+0.375) / (c(1, 2) + 1)
        for (expected, actual) in vec![0.6875, 0.1250]
            .iter()
            .zip(p.index_axis(Axis(0), 28).iter())
        {
            // Assert that the values are equal
            assert!((dbg!(expected) - dbg!(actual)).abs() < 1e-6);
        }
        // only 1 =(c(1, 0)+0.375) / (c(1, 2) + 1)
        for (expected, actual) in vec![0.3125, 0.8750]
            .iter()
            .zip(p.index_axis(Axis(0), 29).iter())
        {
            // Assert that the values are equal
            assert!((dbg!(expected) - dbg!(actual)).abs() < 1e-6);
        }
        for (expected, actual) in vec![0.06250000, 0.04166667]
            .iter()
            .zip(p.index_axis(Axis(0), 63).iter())
        {
            // Assert that the values are equal
            assert!((dbg!(expected) - dbg!(actual)).abs() < 1e-6);
        }
    }
}
