use extendr_api::prelude::*;

pub(crate) mod kmers;
/// Count k-mers across sequences.
/// @export
#[extendr]
fn detect_kmers_across_sequences(sequences: &[Rstr], kmer_size: u32) -> Robj {
    let sequences = sequences.iter().map(|seq| seq.as_bytes());
    let kmer_counts = kmers::detect_kmers_across_sequences(sequences, kmer_size);
    kmer_counts.try_into().expect("Valid matrix")
}
// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod rustyphylotyper;
    fn detect_kmers_across_sequences;
}
