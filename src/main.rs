use clap::Parser;

use rand::Rng;

mod aldknni;
mod filter_missing;
mod geno;
mod helpers;
mod mvi;
mod phen;
mod structs_and_traits;
mod sync;
mod vcf;
use crate::aldknni::*;
use crate::geno::*;
use crate::mvi::*;
use crate::structs_and_traits::*;
use crate::sync::*;
use crate::vcf::*;

#[derive(Parser, Debug)]
#[clap(
    author = "Jeff Paril",
    version = "1.0.0",
    about = "Impute allele frequencies to reduce sparsity of genotype data from polyploids, pooled individuals, and populations.",
    long_about = "Imputation of genotype data from sequencing of more than 2 sets of genomes, i.e. polyploid individuals, population samples, or pools of individuals. This library can also perform simple genotype data filtering prior to imputation. Two imputation methods are available: (1) mean value imputation which uses the arithmentic mean of the locus across non-missing pools; (2) adaptive linkage-informed k-nearest neighbour imputation. This is an attempt to extend the [LD-kNNi method of Money et al, 2015, i.e. LinkImpute](https://doi.org/10.1534/g3.115.021667), which was an extension of the [kNN imputation of Troyanskaya et al, 2001](https://doi.org/10.1093/bioinformatics/17.6.520). Similar to LD-kNNi, LD is estimated using Pearson's product moment correlation across loci per pair of samples. Mean absolute difference in allele frequencies is used to define genetic distance between samples, instead of taxicab or Manhattan distance in LD-kNNi. Four parameters can be set by the user, (1) minimum loci correlation threshold: dictates the minimum LD between the locus requiring imputation and other loci which will be used to estimate genetic distance between samples; (2) maximum genetic distance threshold: sets the maximum genetic distance between the sample requiring imputation and the samples (i.e. nearest neighbours) to be used in weighted mean imputation of missing allele frequencies; (3) minimum number of loci linked to the locus requiring imputation: overrides minimum loci correlation threshold if this minimum is not met; and (4) minimum k-nearest neighbours: overrides maximum genetic distance threshold if this minimum is not met. The first two parameters (minimum loci correlation and maximum genetic distance thresholds) can be optimised per locus requiring imputation using non-missing samples as replicates simulating missing data to minimum the mean absolute error in imputation."
)]
struct Args {
    /// Filename of the genotype file to be imputed. The format can be an uncompressed [vcf](https://github.com/jeffersonfparil/imputef/tree/main?tab=readme-ov-file#variant-call-format-vcf), [sync](https://github.com/jeffersonfparil/imputef/tree/main?tab=readme-ov-file#synchronised-pileup-sync), or [allele frequency table](https://github.com/jeffersonfparil/imputef/tree/main?tab=readme-ov-file#allele-frequency-table-tsv)
    #[clap(short, long)]
    fname: String,
    /// Imputation method. Use "mean" for mean value imputation or "aldknni" for adaptive LD-kNN imputation
    #[clap(short, long, default_value = "aldknni")]
    method: String,
    /// Minimum coverage per locus, i.e. if at a locus, a pool falls below this value (does not skip missing data, i.e. missing locus has a depth of zero), then the whole locus is omitted.
    /// Set this to zero if the vcf has been filtered and contains missing values, i.e. `./.` or `.|.`.
    #[clap(long, default_value_t = 0)]
    min_coverage: u64,
    /// Minimum allele frequency per locus, i.e. if at a locus, a pool has all its alleles below this value and/or above the additive complement of this value (skipping missing data), then the entire locus is omitted.
    #[clap(long, default_value_t = 0.0001)]
    min_allele_frequency: f64,
    /// Maximum fraction of pools missing per locus, i.e. if at a locus, there were more pools missing than the coverage dictated by this threshold, then the locus is omitted.
    #[clap(long, default_value_t = 1.00)]
    max_missingness_rate_per_locus: f64,
    /// Vector of pool sizes, i.e. the number of individuals included in each pool. Enter the pool sizes separated by commas `,`.
    /// This can also be set to a single arbitrarily large value like 100 for individual polyploids or if allele frequency estimates are expected to be accurate.
    #[clap(long, value_parser, value_delimiter = ',', default_value = "100.0")]
    pool_sizes: Vec<f64>,
    /// Minimum depth at which loci with depth below this threshold are set to missing. Set to one if the input vcf has already been filtered and the loci beyond the depth thresholds have been set to missing, otherwise set to an integer above zero.
    #[clap(long, default_value_t = 1.00)]
    min_depth_below_which_are_missing: f64,
    /// Maximum depth at which loci with depth above this threshold are set to missing. Set to some large arbitrarily large value (e.g. 1000000) if the input vcf has already been filtered and the loci beyond the depth thresholds have been set to missing, otherwise set to an integer above zero.
    #[clap(long, default_value_t = 1000000.0)]
    max_depth_above_which_are_missing: f64,
    /// Fraction of pools with the highest number of missing loci to be omitted. Set to zero if the input vcf has already been filtered and the loci beyond the depth thresholds have been set to missing, otherwise set to a decimal number between zero and one.
    #[clap(long, default_value_t = 0.0)]
    frac_top_missing_pools: f64,
    /// Fraction of loci with the highest number of pools with missing data to be omitted. Set to zero if the input vcf has already been filtered and the loci beyond the depth thresholds have been set to missing, otherwise set to an decimal number between zero and one.
    #[clap(long, default_value_t = 0.0)]
    frac_top_missing_loci: f64,
    /// Minimum correlation (Pearson's correlation) between the locus requiring imputation and other loci deemed to be in linkage with it.
    /// Ranges from 0.0 to 1.0, but use -1.0 or any negative number to perform per locus optimisations to find the best value minimising imputation.
    #[clap(long, default_value_t = 0.9)]
    min_loci_corr: f64,
    /// Maximum genetic distance (mean absolute difference in allele frequencies) between the pool or sample requiring imputation and pools or samples deemed to be the closest neighbours.
    /// Ranges from 0.0 to 1.0, but use -1.0 or any negative number to perform per locus optimisations to find the best value minimising imputation.
    #[clap(long, default_value_t = 0.1)]
    max_pool_dist: f64,
    /// Minimum number of linked loci to be used in estimating genetic distances between the pool or sample requiring imputation and other pools or samples (minimum value of 1).
    /// This argument overrides `min_loci_corr`, i.e. the minimum number of loci will be met regardless of the minimum loci correlation threshold.
    #[clap(long, default_value_t = 20)]
    min_l_loci: usize,
    /// Minimum number of k-nearest neighbours of the pool or sample requiring imputation (minimum value of 1).
    /// This argument overrides `max_pool_dist`, i.e. the minimum number of k-nearest neighbours will be met regardless of the maximum genetic distance threshold.
    #[clap(long, default_value_t = 5)]
    min_k_neighbours: usize,
    /// Restrict the choice of linked loci to within the chromosome the locus requiring imputation belongs to? [default: false]
    #[clap(long, action)]
    restrict_linked_loci_per_chromosome: bool,
    /// Number of replications for the estimation of imputation accuracy in terms of mean absolute error (MAE). It is used to define the number of random non-missing samples to use as replicates for the estimation of MAE and optimisation (minimum value of 1).
    #[clap(long, default_value_t = 10)]
    n_reps: usize,
    /// Number of computing threads or processor cores to use in the computations.
    #[clap(long, default_value_t = 2)]
    n_threads: usize,
    /// Prefix of the output files including the [imputed allele frequency table](#allele-frequency-table-tsv) (`<fname_out_prefix>-<time>-<random_id>-IMPUTED.tsv`).
    #[clap(long, default_value = "")]
    fname_out_prefix: String,
}

/// # imputef: Impute allele frequencies to reduce sparsity of genotype data from polyploids, pooled individuals, and populations.
/// Imputation of genotype data from sequencing of more than 2 sets of genomes, i.e. polyploid individuals, population samples, or pools of individuals.
/// Two imputation methods are available:
/// 1. mean value imputation which uses the arithmetic mean of the locus across non-missing pools (`?imputef::mvi`),
/// 2. adaptive linkage-informed k-nearest neighbour imputation (`?imputef::aldknni`).
/// This is an attempt to extend the [LD-kNNi method of Money et al, 2015, i.e. LinkImpute](https://doi.org/10.1534/g3.115.021667), which was an extension of the [kNN imputation of Troyanskaya et al, 2001](https://doi.org/10.1093/bioinformatics/17.6.520).
/// Similar to LD-kNNi, LD is estimated using Pearson's product moment correlation across loci per pair of samples.
/// Mean absolute difference in allele frequencies is used to define genetic distance between samples, instead of taxicab or Manhattan distance in LD-kNNi.
/// Four parameters can be set by the user:
/// 1. minimum loci correlation threshold: dictates the minimum LD between the locus requiring imputation and other loci which will be used to estimate genetic distance between samples,
/// 2.  maximum genetic distance threshold: sets the maximum genetic distance between the sample requiring imputation and the samples (i.e. nearest neighbours) to be used in weighted mean imputation of missing allele frequencies,
/// 3. minimum number of loci linked to the locus requiring imputation: overrides minimum loci correlation threshold if this minimum is not met, and
/// 4. minimum k-nearest neighbours: overrides maximum genetic distance threshold if this minimum is not met.
/// The first two parameters (minimum loci correlation and maximum genetic distance thresholds) can be optimised per locus requiring imputation using non-missing samples as replicates simulating missing data to minimum the mean absolute error in imputation.
/// This library can also perform simple genotype data filtering prior to imputation.
fn main() {
    let args = Args::parse();
    // Identify the format of the input file
    // 1) vcf - richest (*.vcf)
    // 2) sync - intermediate richness and most preferred (*.sync)
    // 3) geno - least detailed - tab-delimited: chr,pos,allele,sample-1,sample-2,some-name-@#@#$%^&*(+)}:<'?"-with-a-bunch-of-asci-characters,... (*.txt)
    let extension_name: &str = args
        .fname
        .split('.')
        .collect::<Vec<&str>>()
        .last()
        .expect("Error extracting the last character of the input filename in impute().");
    println!("##################################################");
    if extension_name == "vcf" {
        println!("Input uncompressed vcf file:\n{:?}", args.fname);
    } else if extension_name == "sync" {
        println!("Input sync file:\n{:?}", args.fname);
    } else {
        println!("Input tsv file:\n{:?}", args.fname);
    }
    println!(
        "Minimum coverage per locus across pools:\n{:?}",
        args.min_coverage
    );
    println!("Minimum allele frequency:\n{:?}", args.min_allele_frequency);
    println!(
        "Maximum fraction of pools missing at each locus:\n{:?}",
        args.max_missingness_rate_per_locus
    );
    println!(
        "Restrict linked loci per chromosome:\n{:?}",
        args.restrict_linked_loci_per_chromosome
    );
    println!("Pool sizes:\n{:?}", args.pool_sizes);
    println!(
        "Number of threads to use for parallel processing:\n{:?}",
        args.n_threads
    );
    // Create a random ID by concatenating the time in seconds with some trailing random numbers for extra safe write-outs (extra safe because Rust will complain if the file which we want to write on already exists)
    let time = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .expect("Error extracting time in UNIX_EPOCH within impute().")
        .as_secs_f64();
    let mut rng = rand::thread_rng();
    let num = rng.gen::<u32>();
    let rand_id = time.to_string() + "." + &num.to_string();
    println!("rand_id={:?}", rand_id);
    // Prepare filtering variables
    let mut filter_stats = FilterStats {
        remove_ns: true,
        min_quality: 0.01, // for pileups only (i.e. minimum base error rate)
        min_coverage: args.min_coverage,
        min_allele_frequency: args.min_allele_frequency,
        max_missingness_rate: args.max_missingness_rate_per_locus,
        pool_sizes: args.pool_sizes.clone(),
    };
    // Load genotype data
    let (mut genotypes_and_phenotypes, _) = if extension_name == "vcf" {
        load_vcf(
            &args.fname,
            &mut filter_stats,
            &args.fname_out_prefix,
            &rand_id,
            &args.n_threads,
        )
        .expect("Error loading vcf.")
    } else if extension_name == "sync" {
        load_sync(
            &args.fname,
            &mut filter_stats,
            &args.fname_out_prefix,
            &rand_id,
            &args.n_threads,
        )
        .expect("Error loading sync.")
    } else {
        load_geno(
            &args.fname,
            &mut filter_stats,
            &args.fname_out_prefix,
            &rand_id,
            &args.n_threads,
        )
        .expect("Error loading sync.")
    };
    // Define missing data
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .set_missing_by_depth(
            &args.min_depth_below_which_are_missing,
            &args.max_depth_above_which_are_missing,
        )
        .expect("Error calling set_missing_by_depth() method within impute().");
    let end = std::time::SystemTime::now();
    let duration = end
        .duration_since(start)
        .expect("Error measuring the duration of setting missing data within impute().");
    println!(
        "Set loci beyond the minimum ({}X) and maximum ({}X) depth thresholds to missing: {} pools x {} loci | Missingness: {}% | Duration: {} seconds",
        &args.min_depth_below_which_are_missing,
        &args.max_depth_above_which_are_missing,
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        genotypes_and_phenotypes.missing_rate().expect("Error measuring sparsity via missing_rate() method after setting missing by depth within impute()."),
        duration.as_secs_f64()
    );
    // Filter pools
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .filter_out_top_missing_pools(&args.frac_top_missing_pools)
        .expect("Error filtering out top-most missing pools within impute().");
    let end = std::time::SystemTime::now();
    let duration = end
        .duration_since(start)
        .expect("Error measuring the duration of filtering pools within impute().");
    println!(
        "Filtered out sparsest pools: {} pools x {} loci | Missingness: {}% | Duration: {} seconds",
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        genotypes_and_phenotypes.missing_rate().expect("Error measuring sparsity via missing_rate() method after filtering pools within impute()."),
        duration.as_secs_f64()
    );
    // Filter loci
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .filter_out_top_missing_loci(&args.frac_top_missing_loci)
        .expect("Error filtering out top-most missing loci within impute().");
    let end = std::time::SystemTime::now();
    let duration = end
        .duration_since(start)
        .expect("Error measuring the duration of filtering loci within impute().");
    println!(
        "Filtered out sparsest loci: {} pools x {} loci | Missingness: {}% | Duration: {} seconds",
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        genotypes_and_phenotypes.missing_rate().expect("Error measuring sparsity via missing_rate() method after filtering loci within impute()."),
        duration.as_secs_f64()
    );
    // Prepare output file name
    let fname_out = if args.fname_out_prefix == *"" {
        args.fname.to_owned() + "-" + &rand_id + "-IMPUTED.tsv"
    } else {
        args.fname_out_prefix.to_owned() + "-" + &rand_id + "-IMPUTED.tsv"
    };
    let _ = if args.method == "mean" {
        println!("###################################################################################################");
        println!("mvi: mean value imputation");
        println!("###################################################################################################");
        impute_mean(
            genotypes_and_phenotypes,
            &filter_stats,
            &args.n_threads, // not used in mvi
            &fname_out,
        )
        .expect("Error performing mean value imputation via impute_mean() within impute().")
    } else {
        println!("###################################################################################################");
        println!("aldknni: adaptive linkage disequilibrium (LD)-based k-nearest neighbour imputation of genotype data");
        println!("###################################################################################################");
        // Handling NA conversion into Rust's f64:NAN
        let min_loci_corr = if args.min_loci_corr < 0.0 {
            f64::NAN
        } else {
            args.min_loci_corr
        };
        let max_pool_dist = if args.max_pool_dist < 0.0 {
            f64::NAN
        } else {
            args.max_pool_dist
        };
        // Use a single rep for estimating imputation accuracy if not optimising for min_loci_corr and/or max_pool_dist thresholds
        let n_reps = if !min_loci_corr.is_nan() && !max_pool_dist.is_nan() {
            1
        } else {
            args.n_reps
        };
        impute_aldknni(
            genotypes_and_phenotypes,
            &filter_stats,
            (
                &min_loci_corr,
                &max_pool_dist,
                &args.min_l_loci,
                &args.min_k_neighbours,
                &n_reps,
            ),
            args.restrict_linked_loci_per_chromosome,
            &args.n_threads,
            &fname_out,
        )
        .expect("Error performing adaptive LD-kNN imputation via impute_aldknni() within impute().")
    };
    println!(
        "Imputation output in allele frequency table format: {}",
        fname_out
    );
}

// # Running tests
// cargo fmt
// cargo fix --allow-dirty
// cargo clippy
// time cargo test
// rm intermediate_output-* test-* tests/test-* tests/test.*-*.* tests/test_2.*-*.* tests/test_2.*-*.*
// # Compilation and building docs (for x86_64-unknown-linux-gnu)
// cargo build --release
// RUSTDOCFLAGS="--html-in-header res/doc_header_for_maths.html" cargo doc --no-deps --document-private-items ### open: $(pwd)/target/doc/imputef/index.html
// # Compilation on Windows (for x86_64-pc-windows-gnu)
// Install git and rustup on Windows
// git clone https://github.com/jeffersonfparil/imputef.git
// cd imputef
// cargo build --release
// Compile on macOS
// Install Docker Desktop on Windows and activate WSL2 backend
// Open WSL2 and install: sudo apt -y install bridge-utils cpu-checker libvirt-clients libvirt-daemon qemu qemu-kvm
// docker pull sickcodes/docker-osx:latest
// Open Docker Desktop, find the downloaded docker image and run
// git clone https://github.com/jeffersonfparil/imputef.git
// cd imputef
// cargo build --release
