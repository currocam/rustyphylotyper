use extendr_api::prelude::*;
mod kmers;
/// Temporary function that returns the indexes (in base4) of the kemers.
/// @export
#[extendr]
fn kmers(sequence: &str, k: i32) -> Vec<usize> {
    kmers::kmers(sequence.as_bytes(), k as usize).collect()
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod rustyphylotyper;
    fn kmers;
}
