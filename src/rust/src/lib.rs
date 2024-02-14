use extendr_api::prelude::*;
use ndarray::prelude::*;
use rand::Rng;
use std::fs::File;
use std::io::{self, prelude::*, BufReader};

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
use crate::mvi::*;
use crate::structs_and_traits::*;
use crate::vcf::*;

#[extendr]
#[allow(clippy::too_many_arguments)]
fn impute(
    fname: String,
    imputation_method: String,
    min_coverage: u64,
    min_allele_frequency: f64,
    max_missingness_rate_per_locus: f64,
    pool_sizes: Vec<f64>,
    min_depth_below_which_are_missing: f64,
    max_depth_above_which_are_missing: f64,
    frac_top_missing_pools: f64,
    frac_top_missing_loci: f64,
    n_reps: u64,
    min_loci_corr: f64,
    max_pool_dist: f64,
    min_l_loci: u64,
    min_k_neighbours: u64,
    restrict_linked_loci_per_chromosome: bool,
    n_threads: u64,
    fname_out_prefix: String,
    // ) -> String {
) -> Robj {
    // Identify the format of the input file
    // 1) vcf - richest (*.vcf)
    // 2) sync - intermediate richness and most preferred (*.sync)
    // 3) geno - least detailed - tab-delimited: chr,pos,allele,sample-1,sample-2,some-name-@#@#$%^&*(+)}:<'?"-with-a-bunch-of-asci-characters,... (*.txt)
    let extension_name: &str = fname
        .split('.')
        .collect::<Vec<&str>>()
        .last()
        .expect("Error extracting the last character of the input filename in impute().");
    println!("##################################################");
    if extension_name == "vcf" {
        println!("Input uncompressed vcf file:\n{:?}", fname);
    } else if extension_name == "sync" {
        println!("Input sync file:\n{:?}", fname);
    } else {
        println!("Input csv file:\n{:?}", fname);
    }
    println!(
        "Minimum coverage per locus across pools:\n{:?}",
        min_coverage
    );
    println!("Minimum allele frequency:\n{:?}", min_allele_frequency);
    println!(
        "Maximum fraction of pools missing at each locus:\n{:?}",
        max_missingness_rate_per_locus
    );
    println!("Pool sizes:\n{:?}", pool_sizes);
    println!(
        "Number of threads to use for parallel processing:\n{:?}",
        n_threads
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
        min_coverage,
        min_allele_frequency,
        max_missingness_rate: max_missingness_rate_per_locus,
        pool_sizes: pool_sizes.clone(),
    };
    let keep_p_minus_1 = false;
    // Load genotype data
    let mut genotypes_and_phenotypes: GenotypesAndPhenotypes = if extension_name == "vcf" {
        // Extract pool names from the vcf file
        let mut pool_names: Vec<String> = vec![];
        let file = File::open(fname.clone()).expect("Error opening the input vcf file.");
        let reader = BufReader::new(file);
        for l in reader.lines() {
            let mut line = l.expect("Error reading the input vcf file.");
            // Remove trailing newline character in Unix-like (\n) and Windows (\r)
            if line.ends_with('\n') {
                line.pop();
                if line.ends_with('\r') {
                    line.pop();
                }
            }
            let vec_line: Vec<&str> = line.split('\t').collect();
            if vec_line[0] == "#CHROM" {
                pool_names = vec_line[9..vec_line.len()]
                    .iter()
                    .map(|&x| x.to_owned())
                    .collect();
                break;
            }
        }
        let n = pool_names.len();
        // If a single pool size was supplied then we are assuming the same sizes across all pools
        let pool_sizes = if pool_sizes.len() == 1 {
            vec![pool_sizes[0]; n]
        } else {
            pool_sizes
        };
        filter_stats.pool_sizes = pool_sizes;
        // Convert vcf into sync
        let file_vcf = FileVcf {
            filename: fname.to_owned(),
        };
        let fname_sync_out = if fname_out_prefix == *"" {
            fname.to_owned() + "-" + &rand_id + ".sync"
        } else {
            fname_out_prefix.to_owned() + "-" + &rand_id + ".sync"
        };
        let _ = file_vcf
            .read_analyse_write(
                &filter_stats,
                &fname_sync_out,
                &(n_threads as usize),
                vcf_to_sync,
            )
            .expect("Error converting the vcf into sync via read_analyse_write() method within impute().");
        // Initialise dummy phen struct (for other poolgen analyses)
        let file_sync_phen = FileSyncPhen {
            filename_sync: fname_sync_out.to_owned(),
            pool_names,
            pool_sizes: filter_stats.pool_sizes.clone(),
            phen_matrix: Array2::from_elem((n, 1), f64::NAN),
            test: "".to_owned(),
        };
        file_sync_phen
            .convert_into_genotypes_and_phenotypes(&filter_stats, keep_p_minus_1, &(n_threads as usize))
            .expect("Error parsing the input genotype (converted from vcf into sync) and dummy phenotype data via convert_into_genotypes_and_phenotypes() method within impute().")
    } else if extension_name == "sync" {
        // Extract pool names from the sync file
        let mut pool_names: Vec<String> = vec![];
        let file = File::open(fname.clone()).expect("Error reading the input vcf file.");
        let reader = BufReader::new(file);
        for l in reader.lines() {
            let mut line = l.expect("Error reading the input sync file.");
            // Remove trailing newline character in Unix-like (\n) and Windows (\r)
            if line.ends_with('\n') {
                line.pop();
                if line.ends_with('\r') {
                    line.pop();
                }
            }
            let vec_line: Vec<&str> = line.split('\t').collect();
            if vec_line[0] == "#chr" {
                pool_names = vec_line[3..vec_line.len()]
                    .iter()
                    .map(|&x| x.to_owned())
                    .collect();
                break;
            }
        }
        let n = pool_names.len();
        // If a single pool size was supplied then we are assuming the same sizes across all pools
        let pool_sizes = if pool_sizes.len() == 1 {
            vec![pool_sizes[0]; n]
        } else {
            pool_sizes
        };
        filter_stats.pool_sizes = pool_sizes;
        // Initialise dummy phen struct (for other poolgen analyses)
        let file_sync_phen = FileSyncPhen {
            filename_sync: fname.to_owned(),
            pool_names,
            pool_sizes: filter_stats.pool_sizes.clone(),
            phen_matrix: Array2::from_elem((n, 1), f64::NAN),
            test: "".to_owned(),
        };
        file_sync_phen
            .convert_into_genotypes_and_phenotypes(&filter_stats, keep_p_minus_1, &(n_threads as usize))
            .expect("Error parsing the genotype (sync format) and dummy phenotype data via convert_into_genotypes_and_phenotypes() method within impute().")
    } else {
        // Extract pool names from the txt file
        let file: File =
            File::open(fname.clone()).expect("Error reading the allele frequency table file.");
        let reader = io::BufReader::new(file);
        let mut header: String = reader
            .lines()
            .next()
            .expect("Error reading the allele frequency table file.")
            .expect("Please check the format of the allele frequency table text file.");
        if header.ends_with('\n') {
            header.pop();
            if header.ends_with('\r') {
                header.pop();
            }
        }
        let vec_header: Vec<&str> = header.split('\t').collect();
        let vec_header: Vec<&str> = if vec_header.len() == 1 {
            header.split(',').collect()
        } else {
            vec_header
        };
        let vec_header: Vec<&str> = if vec_header.len() == 1 {
            header.split(';').collect()
        } else {
            vec_header
        };
        let pool_names: Vec<String> = vec_header[3..vec_header.len()]
            .iter()
            .map(|&x| x.to_owned())
            .collect();
        let n = pool_names.len();
        // If a single pool size was supplied then we are assuming the same sizes across all pools
        let pool_sizes = if pool_sizes.len() == 1 {
            vec![pool_sizes[0]; n]
        } else {
            pool_sizes
        };
        let file_geno = FileGeno {
            filename: fname.clone(),
        };
        filter_stats.pool_sizes = pool_sizes;
        file_geno
            .convert_into_genotypes_and_phenotypes(&filter_stats, keep_p_minus_1, &(n_threads as usize))
            .expect("Error parsing the genotype data (extracted from allele frequency table text file) via convert_into_genotypes_and_phenotypes() method within impute().")
    };
    // println!("genotypes_and_phenotypes={:?}", genotypes_and_phenotypes);
    // Define missing data
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .set_missing_by_depth(
            &min_depth_below_which_are_missing,
            &max_depth_above_which_are_missing,
        )
        .expect("Error calling set_missing_by_depth() method within impute().");
    let end = std::time::SystemTime::now();
    let duration = end
        .duration_since(start)
        .expect("Error measuring the duration of setting missing data within impute().");
    println!(
        "Set loci beyond the minimum and maximum depth thresholds to missing: {} pools x {} loci | Missingness: {}% | Duration: {} seconds",
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        genotypes_and_phenotypes.missing_rate().expect("Error measuring sparsity via missing_rate() method after setting missing by depth within impute()."),
        duration.as_secs()
    );
    // Filter pools
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .filter_out_top_missing_pools(&frac_top_missing_pools)
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
        duration.as_secs()
    );
    // Filter loci
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .filter_out_top_missing_loci(&frac_top_missing_loci)
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
        duration.as_secs()
    );
    // Prepare output file name
    let fname_out = if fname_out_prefix == *"" {
        fname.to_owned() + "-" + &rand_id + "-IMPUTED.csv"
    } else {
        fname_out_prefix.to_owned() + "-" + &rand_id + "-IMPUTED.csv"
    };
    let _ = if imputation_method == *"mean" {
        println!("###################################################################################################");
        println!("mvi: mean value imputation");
        println!("###################################################################################################");
        impute_mean(
            genotypes_and_phenotypes,
            &filter_stats,
            &(n_threads as usize),
            &fname_out,
        )
        .expect("Error performing mean value imputation via impute_mean() within impute().")
    } else {
        println!("###################################################################################################");
        println!("aldknni: adaptive linkage disequilibrium (LD)-based k-nearest neighbour imputation of genotype data");
        println!("###################################################################################################");
        // Handling NA conversion into Rust's f64:NAN
        let min_loci_corr = if min_loci_corr < 0.0 {
            f64::NAN
        } else {
            min_loci_corr
        };
        let max_pool_dist = if max_pool_dist < 0.0 {
            f64::NAN
        } else {
            max_pool_dist
        };
        impute_aldknni(
            genotypes_and_phenotypes,
            &filter_stats,
            (
                &min_loci_corr,
                &max_pool_dist,
                &(min_l_loci as usize),
                &(min_k_neighbours as usize),
                &(n_reps as usize),
            ),
            restrict_linked_loci_per_chromosome,
            &(n_threads as usize),
            &fname_out,
        )
        .expect("Error performing adaptive LD-kNN imputation via impute_aldknni() within impute().")
    };
    r!(fname_out.to_owned())
    // fname_out.to_owned()
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod imputef;
    fn impute;
}

// Format, fix and ask clippy:
// ```shell
// cargo fmt
// cargo fix --allow-dirty
// cargo clippy
// ```
// Test in R: `time Rscript imputef/tests/tests.R`
