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
mod optim;
mod phen;
mod structs_and_traits;
mod sync;
mod vcf;
use crate::aldknni::*;

use crate::mvi::*;

use crate::structs_and_traits::*;

use crate::vcf::*;

#[extendr]
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
    window_size_bp: u64,
    min_loci_per_window: u64,
    min_loci_corr: f64,
    max_pool_dist: f64,
    optimise_for_thresholds: bool,
    optimise_n_steps_corr: u64,
    optimise_n_steps_dist: u64,
    optimise_n_reps: u64,
    n_threads: u64,
    fname_out_prefix: String,
    // ) -> String {
) -> Robj {
    // Identify the format of the input file
    // 1) vcf - richest (*.vcf)
    // 2) sync - intermediate richness and most preferred (*.sync)
    // 3) geno - least detailed - tab-delimited: chr,pos,allele,sample-1,sample-2,some-name-@#@#$%^&*(+)}:<'?"-with-a-bunch-of-asci-characters,... (*.txt)
    let extension_name: &str = fname.split(".").collect::<Vec<&str>>().last().unwrap();
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
        .unwrap()
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
        let file = File::open(fname.clone()).unwrap();
        let reader = BufReader::new(file);
        for l in reader.lines() {
            let mut line = l.unwrap();
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
        let _fname_sync_out: String = fname.to_owned() + "-" + &rand_id + ".sync";
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
            .unwrap();
        // Initialise dummy phen struct (for other poolgen analyses)
        let file_sync_phen = FileSyncPhen {
            filename_sync: fname_sync_out.to_owned(),
            pool_names,
            pool_sizes: filter_stats.pool_sizes.clone(),
            phen_matrix: Array2::from_elem((n, 1), f64::NAN),
            test: "".to_owned(),
        };
        file_sync_phen
            .into_genotypes_and_phenotypes(&filter_stats, keep_p_minus_1, &(n_threads as usize))
            .unwrap()
    } else if extension_name == "sync" {
        // Extract pool names from the sync file
        let mut pool_names: Vec<String> = vec![];
        let file = File::open(fname.clone()).unwrap();
        let reader = BufReader::new(file);
        for l in reader.lines() {
            let mut line = l.unwrap();
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
            .into_genotypes_and_phenotypes(&filter_stats, keep_p_minus_1, &(n_threads as usize))
            .unwrap()
    } else {
        // Extract pool names from the txt file
        let file: File = File::open(fname.clone()).unwrap();
        let reader = io::BufReader::new(file);
        let mut header: String = reader.lines().next().unwrap().unwrap();
        if header.ends_with('\n') {
            header.pop();
            if header.ends_with('\r') {
                header.pop();
            }
        }
        let vec_header: Vec<&str> = header.split("\t").collect();
        let vec_header: Vec<&str> = if vec_header.len() == 1 {
            header.split(",").collect()
        } else {
            vec_header
        };
        let vec_header: Vec<&str> = if vec_header.len() == 1 {
            header.split(";").collect()
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
            .into_genotypes_and_phenotypes(&filter_stats, keep_p_minus_1, &(n_threads as usize))
            .unwrap()
    };
    // println!("genotypes_and_phenotypes={:?}", genotypes_and_phenotypes);
    // Define missing data
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .set_missing_by_depth(
            &min_depth_below_which_are_missing,
            &max_depth_above_which_are_missing,
        )
        .unwrap();
    let end = std::time::SystemTime::now();
    let duration = end.duration_since(start).unwrap();
    println!(
        "Set loci beyond the minimum and maximum depth thresholds to missing: {} pools x {} loci | Missingness: {}% | Duration: {} seconds",
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        genotypes_and_phenotypes.missing_rate().unwrap(),
        duration.as_secs()
    );
    // Filter pools
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .filter_out_top_missing_pools(&frac_top_missing_pools)
        .unwrap();
    let end = std::time::SystemTime::now();
    let duration = end.duration_since(start).unwrap();
    println!(
        "Filtered out sparsest pools: {} pools x {} loci | Missingness: {}% | Duration: {} seconds",
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        genotypes_and_phenotypes.missing_rate().unwrap(),
        duration.as_secs()
    );
    // Filter loci
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .filter_out_top_missing_loci(&frac_top_missing_loci)
        .unwrap();
    let end = std::time::SystemTime::now();
    let duration = end.duration_since(start).unwrap();
    println!(
        "Filtered out sparsest loci: {} pools x {} loci | Missingness: {}% | Duration: {} seconds",
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        genotypes_and_phenotypes.missing_rate().unwrap(),
        duration.as_secs()
    );
    // Determine if the input data is diploid biallelic then we use LinkImpute's weighted modal imputation
    let mut do_linkimpute_weighted_mode = true;
    for x in genotypes_and_phenotypes.intercept_and_allele_frequencies.iter() {
        if !(*x).is_nan() {
            if (*x!=0.0) & (*x!=0.5) & (*x!=1.0) {
                do_linkimpute_weighted_mode = false;
                break;
            } else {
                continue
            }
        } else {
            continue;
        }
    }
    if do_linkimpute_weighted_mode {
        println!("The input genotype data is biallelic diploid and will be using weighted modal imputation.");
    } else {
        println!("The input genotype data is not biallelic diploid and will be using weighted mean imputation.");
    }
    // println!("genotypes_and_phenotypes.intercept_and_allele_frequencies={:?}", genotypes_and_phenotypes.intercept_and_allele_frequencies);
    // println!("do_linkimpute_weighted_mode={:?}", do_linkimpute_weighted_mode);
    // Prepare output file name
    let fname_out = if fname_out_prefix == *"" {
        fname.to_owned() + "-" + &rand_id + "-IMPUTED.csv"
    } else {
        fname_out_prefix.to_owned() + "-" + &rand_id + "-IMPUTED.csv"
    };
    let _ = if &imputation_method == &"mean".to_owned() {
        println!("###################################################################################################");
        println!("mvi: mean value imputation");
        println!("###################################################################################################");
        impute_mean(
            genotypes_and_phenotypes,
            &filter_stats,
            &min_depth_below_which_are_missing,
            &max_depth_above_which_are_missing,
            &frac_top_missing_pools,
            &frac_top_missing_loci,
            &(n_threads as usize),
            &fname_out,
        )
        .unwrap()
    } else {
        println!("###################################################################################################");
        println!("aldknni: adaptive linkage disequilibrium (LD)-based k-nearest neighbour imputation of genotype data");
        println!("###################################################################################################");
        impute_aldknni(
            genotypes_and_phenotypes,
            &filter_stats,
            &window_size_bp,
            &window_size_bp, // set the slide size as the window size hence non-overlapping windows
            &min_loci_per_window,
            &min_loci_corr,
            &max_pool_dist,
            &optimise_for_thresholds,
            &(optimise_n_steps_corr as usize),
            &(optimise_n_steps_dist as usize),
            &(optimise_n_reps as usize),
            do_linkimpute_weighted_mode,
            &(n_threads as usize),
            &fname_out,
        )
        .unwrap()
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

// #[cfg(test)]
// mod tests {
//     use super::*;
//     #[test]
//     fn test_lib() {
//         let out = impute(
//             "tests/test.vcf".to_owned(),
//             "mvi".to_owned(),
//             0,
//             0.00001,
//             1.0,
//             vec![100.0],
//             1.0,
//             1_000_000.0,
//             0.0,
//             0.0,
//             1_000_000,
//             1,
//             0.9,
//             0.1,
//             true,
//             2,
//             "".to_owned());
//         assert_eq!(0, 1);
//     }
// }
