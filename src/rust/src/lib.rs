use extendr_api::prelude::*;

pub(crate) mod kmers;
/// Count k-mers across sequences.
/// @export
#[extendr]
fn detect_kmers_across_sequences(sequences: &[Rstr], kmer_size: u32) -> Robj {
    let nrows = 4_usize.pow(kmer_size);
    let ncols = sequences.len();
    let mut kmer_counts = Array2::<bool>::from_elem((nrows, ncols), false);
    kmer_counts
        .axis_iter_mut(Axis(1))
        .zip(sequences.iter())
        .for_each(|(mut col, seq)| {
            for kmer in kmers::kmers(seq.as_bytes(), kmer_size as usize) {
                col[kmer] = true;
            }
        });
    kmer_counts.try_into().expect("Valid matrix")
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod rustyphylotyper;
    fn detect_kmers_across_sequences;
}
