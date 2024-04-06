use extendr_api::prelude::*;

pub(crate) mod database;
pub(crate) mod kmers;
/// Count k-mers across sequences.
/// @export
#[extendr]
fn detect_kmers_across_sequences(sequences: &[Rstr], kmer_size: u32) -> Robj {
    let sequences = sequences.iter().map(|seq| seq.as_str());
    let kmer_counts = kmers::detect_kmers_across_sequences(sequences, kmer_size);
    kmer_counts.try_into().expect("Valid matrix")
}

/// Construct k-mer database
/// @export
#[extendr]
fn kmer_database(sequences: &[Rstr], genera: &[Rstr], kmer_size: u32) -> List {
    let sequences = sequences.iter().map(|seq| seq.as_str());
    let genera = genera.iter().map(|f| f.as_str());
    let db = database::KmerDatabase::build(sequences, genera, kmer_size).expect("Valid database");
    let probs: Robj = db.conditional_probs.try_into().expect("Valid matrix");
    list!(conditional_prob = probs, genera = db.genera)
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod rustyphylotyper;
    fn kmer_database;
}
