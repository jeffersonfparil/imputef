use ndarray::{prelude::*, Zip};
use std::cmp::Ordering;
use std::io;

use crate::helpers::*;
use crate::optim::*;
use crate::structs_and_traits::*;

fn calculate_genomewide_ld(
    intercept_and_allele_frequencies: &Array2<f64>,
) -> io::Result<Vec<Vec<f64>>> {
    // Errors and f64::NAN are all converted into 0.0 for simplicity
    println!("Calculating genomewide correlations between pairs of loci.");
    let (_n, p) = intercept_and_allele_frequencies.dim();
    let mut corr: Vec<Vec<f64>> = vec![vec![]; p - 1];
    let vec_idx: Vec<usize> = (0..(p - 1)).collect();
    Zip::from(&mut corr).and(&vec_idx).par_for_each(|c, &idx| {
        for j in idx..(p - 1) {
            let corr = match pearsons_correlation_pairwise_complete(
                &intercept_and_allele_frequencies.column(j),
                &intercept_and_allele_frequencies.column(j + 1),
            ) {
                Ok(x) => {
                    if x.0.is_nan() {
                        0.0
                    } else {
                        x.0
                    }
                }
                Err(_) => 0.0,
            };
            c.push(corr);
        }
        assert_eq!(
            c.len(),
            (p - (1 + idx)),
            "Error calculating genomewide allele frequency correlations at idx={}.",
            idx
        )
    });
    println!("Finished calculating genomewide correlations between pairs of loci.");
    Ok(corr)
}

fn calculate_genetic_distances_between_pools(
    idx_row: usize,
    linked_loci_idx: &Vec<usize>,
    intercept_and_allele_frequencies: &Array2<f64>,
) -> io::Result<Array1<f64>> {
    // Errors and f64::NAN are all converted into the maximum possible distance of 1.00 for simplicity
    let (n, _p) = intercept_and_allele_frequencies.dim();
    let q: Array1<f64> = intercept_and_allele_frequencies
        .row(idx_row)
        .select(Axis(0), linked_loci_idx);
    let idx: Vec<usize> = (0..n).collect();
    let mut distances: Array1<f64> = Array1::from_elem(n, 1.0);
    Zip::from(&mut distances).and(&idx).par_for_each(|d, &i| {
        *d = if i != idx_row {
            let q1: Array1<f64> = intercept_and_allele_frequencies
                .row(i)
                .select(Axis(0), linked_loci_idx);
            match (&q - &q1)
                .into_iter()
                .filter(|&x| !x.is_nan())
                .collect::<Array1<f64>>()
                .map(|x| x.abs())
                .mean()
            {
                Some(x) => x,
                None => 1.00,
            }
        } else {
            1.00
        };
    });
    Ok(distances)
}

fn find_l_linked_loci(
    idx_col: usize,
    corr: &Vec<Vec<f64>>,
    min_loci_corr: &f64,
    misc_min_l: usize,
) -> io::Result<(Vec<usize>, Vec<f64>)> {
    assert!(
        misc_min_l > 0,
        "Error: the minimum number of linked loci need to be greater than 0."
    );
    let p = corr.len();
    let mut vec_corr: Vec<f64> = vec![];
    let mut vec_idx: Vec<usize> = vec![];
    // Across rows and at the idx_col of the triangular matrix of correlations
    for i in 0..idx_col {
        let c = corr[i][idx_col - i]; // Less the current index as the length of the vectors decreases by 1 each time
        if c >= *min_loci_corr {
            vec_idx.push(i);
        }
        if misc_min_l > 0 {
            vec_corr.push(c);
        }
    }
    // Across columns and at the idx_col row of the triangular matrix of correlations
    for j in 0..(p - idx_col) {
        let c = corr[idx_col][j];
        if c >= *min_loci_corr {
            vec_idx.push(idx_col + j); // Add the current index as the length of the vectors decreases by 1 each time
        }
        if misc_min_l > 0 {
            vec_corr.push(c);
        }
    }
    // If less than the minimum number of loci passed the threshold, then we sort the correlations (decreasing) and pick the top-most correlated loci
    let (vec_idx, vec_corr) = if vec_idx.len() < misc_min_l {
        let mut indices: Vec<usize> = (0..p).collect();
        indices.sort_by(
            |&a, &b| match (vec_corr[b].is_nan(), vec_corr[a].is_nan()) {
                (true, true) => Ordering::Equal,
                (true, false) => Ordering::Greater,
                (false, true) => Ordering::Less,
                (false, false) => vec_corr[b].partial_cmp(&vec_corr[a]).unwrap(),
            },
        );
        let mut indices = indices[0..misc_min_l].to_vec();
        indices.sort();
        let mut corr: Vec<f64> = vec![];
        for i in indices.iter() {
            corr.push(vec_corr[*i]);
        }
        (indices, corr)
    } else {
        (vec_idx, vec_corr)
    };
    Ok((vec_idx, vec_corr))
}

fn find_k_nearest_neighbours(
    distances: &Array1<f64>,
    max_pool_dist: &f64,
    misc_min_k: usize,
    idx_allele_ini: usize,
    idx_allele_fin: usize,
    intercept_and_allele_frequencies: &Array2<f64>,
) -> io::Result<(Vec<usize>, Vec<f64>, Array2<f64>)> {
    assert!(
        misc_min_k > 0,
        "Error: the minimum number of k-nearest neighbours need to be greater than 0."
    );
    let (n, _p) = intercept_and_allele_frequencies.dim();
    assert_eq!(n, distances.len(), "Error - incompatible number of pools in the genotype matrix and computed distances in find_k_nearest_neighbours().");
    let mut idx: Vec<usize> = vec![];
    let mut dist: Vec<f64> = vec![];
    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(
        |&a, &b| match (distances[a].is_nan(), distances[b].is_nan()) {
            (true, true) => Ordering::Equal,
            (true, false) => Ordering::Greater,
            (false, true) => Ordering::Less,
            (false, false) => distances[a].partial_cmp(&distances[b]).unwrap(),
        },
    );
    for i in indices.clone().into_iter() {
        let bool_dist_nan = distances[i].is_nan() == false;
        let bool_dist_thresh = distances[i] <= *max_pool_dist;
        let bool_freq_nan = intercept_and_allele_frequencies[(i, idx_allele_ini)].is_nan() == false;
        if bool_dist_nan {
            if bool_dist_thresh & bool_freq_nan {
                idx.push(i);
                dist.push(distances[i]);
            }
        }
    }
    // If less than the minimum number of neighbours passed the threshold, then we sort the distances (increasing) and pick the nearest neighbours
    if idx.len() < misc_min_k {
        idx = vec![];
        dist = vec![];
        for i in indices.clone().into_iter() {
            let bool_dist_nan = distances[i].is_nan() == false;
            let bool_freq_nan =
                intercept_and_allele_frequencies[(i, idx_allele_ini)].is_nan() == false;
            let bool_dist_cnt = idx.len() < misc_min_k;
            if bool_dist_cnt == false {
                break;
            }
            if bool_dist_nan {
                if bool_freq_nan {
                    idx.push(i);
                    dist.push(distances[i]);
                }
            }
        }
    }
    // Extract non-missing allele frequencies at the locus requiring imputation of the k-nearest neighbours
    let mut freq: Array2<f64> =
        Array2::from_elem((idx.len(), idx_allele_fin - idx_allele_ini), f64::NAN);
    for i in 0..idx.len() {
        for j in 0..freq.ncols() {
            let idx_row = idx[i];
            let idx_col = idx_allele_ini + j;
            freq[(i, j)] = intercept_and_allele_frequencies[(idx_row, idx_col)];
        }
    }
    Ok((idx, dist, freq))
}

fn impute_allele_frequencies(
    frequencies: &Array2<f64>,
    distances: &Vec<f64>,
    do_linkimpute_weighted_mode: bool,
) -> io::Result<Vec<f64>> {
    let (n, p) = frequencies.dim();
    for i in 0..n {
        for j in 0..p {
            assert!(
                frequencies[(i, j)].is_nan() == false,
                "Error: There should be no missing allele frequencies here."
            );
        }
    }
    let mut imputed_freqs = vec![0.0; p];
    // LD-kNN imputations (weighted mode and mean)
    if do_linkimpute_weighted_mode {
        // Perform weighted modal imputation as in LinkImpute for biallelic diploids - the only 2 differences are that we are performing this per chromosome and the distance metric is MAE rather than Manhattan distance
        assert_eq!(frequencies.ncols(), 1, "Error in the number of alleles per locus. We expect a biallelic locus, please remove the alternative allele or set do_linkimpute_weighted_mode to false.");
        let vec_geno = vec![0.0, 0.5, 1.0];
        let mut max_score = 0.0;
        let mut weighted_mode = 0.0;
        for j in 0..vec_geno.len() {
            let a = vec_geno[j];
            let mut score = 0.0;
            for i in 0..frequencies.column(0).len() {
                let f = 1.00 / (&distances[i] + f64::EPSILON);
                let g = if frequencies.column(0)[i] == a {
                    1.0
                } else {
                    0.0
                };
                score += f * g;
            }
            if score > max_score {
                max_score = score;
                weighted_mode = vec_geno[j];
            }
        }
        // println!("weighted_mode={:?}", weighted_mode);
        imputed_freqs.push(weighted_mode);
    } else {
        // Perform weighted mean allele frequencies
        let additive_inverse_plus_epsilon: Array1<f64> =
            Array1::from_iter(distances.iter().map(|&x| (1.0 - x) + f64::EPSILON)); // Add a small value if all distances are 1.0 and allow equal contributions in such cases
        let weights_sum: f64 = additive_inverse_plus_epsilon
            .iter()
            .fold(0.0, |sum, x| sum + x);
        let weights: Array1<f64> = additive_inverse_plus_epsilon / weights_sum;
        for i in 0..n {
            for j in 0..p {
                imputed_freqs[j] += weights[i] * frequencies[(i, j)];
            }
        }
        if imputed_freqs[0].is_nan() {
            println!("frequencies={:?}", frequencies);
            println!("distances={:?}", distances);
        }
        // Correct allele frequencies so that they sum up to 1, if we have more than 1 allele present
        if p > 1 {
            let sum = imputed_freqs.iter().fold(0.0, |sum, &x| sum + x);
            if sum != 1.0 {
                for j in 0..p {
                    imputed_freqs[j] = imputed_freqs[j] / sum;
                }
            }
        }
    }
    Ok(imputed_freqs)
}

impl GenotypesAndPhenotypes {
    pub fn adaptive_ld_knn_imputation(
        &mut self,
        window_size_bp: &u64,
        _window_slide_size_bp: &u64,
        _min_loci_per_window: &u64,
        min_loci_corr: &f64,
        max_pool_dist: &f64,
        _show_mvi_revert_stats: bool,
        do_linkimpute_weighted_mode: bool,
        misc_min_l: &u64,
        misc_min_k: &u64,
    ) -> io::Result<&mut Self> {
        self.check().expect("Error self.check() within adaptive_ld_knn_imputation() method of GenotypesAndPhenotypes struct.");
        // We are assuming that all non-zero alleles across pools are kept, i.e. biallelic loci have 2 columns, triallelic have 3, and so on.
        let (n, p) = self.intercept_and_allele_frequencies.dim();
        let misc_min_l = if *misc_min_l >= (p as u64) {
            p - 1
        } else {
            *misc_min_l as usize
        };
        let misc_min_k = if *misc_min_k >= (n as u64) {
            n - 1
        } else {
            *misc_min_k as usize
        };
        // Define sliding windows
        let _window_size_bp = if *window_size_bp == 0 {
            // If window size is not set, then use the length of the largest chromosome or scaffold
            self.position
                .iter()
                .fold(0, |max_len, &x| if x > max_len { x } else { max_len })
        } else {
            *window_size_bp
        };
        // Extract loci indices
        let (loci_idx, _loci_chr, _loci_pos) = self.count_loci().expect("Error calling count_loci() method within adaptive_ld_knn_imputation() method for GenotypesAndPhenotypes struct.");
        // Calculate LD across the entire genome
        let corr = calculate_genomewide_ld(&self.intercept_and_allele_frequencies)
            .expect("Error estimating pairwise linkage between loci across the entire genome.");
        // Noting actual correlation and distance parameters used
        let mut actual_corr: Vec<f64> = vec![];
        let mut actual_dist: Vec<f64> = vec![];
        let mut actual_l: Vec<usize> = vec![];
        let mut actual_k: Vec<usize> = vec![];
        // Iterative imputation per locus
        for idx_locus_major_allele in 0..(loci_idx.len() - 1) {
            // Index of the major allele of the current locus
            let j = loci_idx[idx_locus_major_allele];
            // Index of the last allele of the current locus
            let j1 = loci_idx[idx_locus_major_allele + 1];
            // Identify pools requiring imputation at the current locus
            let mut pools_idx: Vec<usize> = vec![];
            for i in 0..n {
                if self.intercept_and_allele_frequencies[(i, j)].is_nan() {
                    pools_idx.push(i);
                }
            }
            if pools_idx.len() == 0 {
                continue;
            }
            // Find loci most correlated to the major allele of the current locus, i.e. the first allele of the locus as they were sorted by decreasing allele frequency (see Sort trait)
            let (linked_loci_idx, correlations) =
                find_l_linked_loci(j, &corr, min_loci_corr, misc_min_l).expect("Error calling find_l_linked_loci() within adaptive_ld_knn_imputation() method for GenotypesAndPhenotypes trait.");
            // Iterate across pools requiring imputation
            for i in pools_idx.iter() {
                // Using the linked loci, estimate the pairwise genetic distance between the current pool and the other pools
                let distances_all_loci = calculate_genetic_distances_between_pools(
                    *i,
                    &linked_loci_idx,
                    &self.intercept_and_allele_frequencies)
                    .expect("Error calling calculate_genetic_distances_between_pools() within adaptive_ld_knn_imputation() method for GenotypesAndPhenotypes trait.");
                // Find the k-nearest neighbours given the maximum distance and/or minimum k-neighbours (shadowing the distances across all pools with distances across k-nearest neighbours)
                let (_idx_neighbours, distances, frequencies) =
                    find_k_nearest_neighbours(
                        &distances_all_loci,
                        max_pool_dist,
                        misc_min_k,
                        j,
                        j1,
                        &self.intercept_and_allele_frequencies,
                    )
                    .expect("Error calling find_k_nearest_neighbours() within adaptive_ld_knn_imputation() method for GenotypesAndPhenotypes trait.");
                // if self.position[j] == 2309794 {
                //     println!("pools_idx={:?}", pools_idx); 
                //     println!("distances_all_loci={:?}", distances_all_loci); 
                //     println!("distances={:?}", distances); 
                //     println!("self.intercept_and_allele_frequencies[(*i, j)]={:?}", self.intercept_and_allele_frequencies[(*i, j)]); 
                //     println!("self.chromosome[*i]={:?}", self.chromosome[*i]); 
                //     println!("self.position[*i]={:?}", self.position[*i]); 
                // }
                // Impute missing allele frequencies at the current locus
                let imputed_freq = impute_allele_frequencies(&frequencies, &distances, do_linkimpute_weighted_mode).expect("Error calling impute_allele_frequencies() within adaptive_ld_knn_imputation() method for GenotypesAndPhenotypes trait.");
                let idx_alleles: Vec<usize> = (j..j1).collect();
                assert_eq!(idx_alleles.len(), imputed_freq.len(), "Error: the number of allele expected to be imputed and number of imputed allele frequencies do not match.");
                for imputed_idx in 0..imputed_freq.len() {
                    let j = idx_alleles[imputed_idx];
                    // if self.position[j] == 2309794 {
                    //     println!("idx_alleles={:?}", idx_alleles); 
                    //     println!("pools_idx={:?}", pools_idx); 
                    //     println!("distances_all_loci={:?}", distances_all_loci); 
                    //     println!("distances={:?}", distances); 
                    //     println!("imputed_freq={:?}", imputed_freq); 
                    //     println!("self.intercept_and_allele_frequencies[(*i, j)]={:?}", self.intercept_and_allele_frequencies[(*i, j)]); 
                    //     println!("self.chromosome[j]={:?}", self.chromosome[j]); 
                    //     println!("self.position[j]={:?}", self.position[j]); 
                    // }
                    self.intercept_and_allele_frequencies[(*i, j)] = imputed_freq[imputed_idx];
                }
                // Convert imputed coverages to f64::INFINITY
                self.coverages[(*i, idx_locus_major_allele)] = f64::INFINITY;
                // Noting actual correlation and distance parameters used
                actual_corr.push(
                    correlations
                        .iter()
                        .fold(0.0, |max, &x| if x > max { x } else { max }),
                );
                actual_dist.push(
                    distances
                        .iter()
                        .fold(0.0, |min, &x| if x > min { x } else { min }),
                );
                actual_l.push(linked_loci_idx.len());
                actual_k.push(distances.len());
            }
        }
        // Report actual correlation and distance parameters used
        let (actual_corr_min, actual_corr_mean, actual_corr_max) = match summary_stats(&actual_corr)
        {
            Ok(x) => x,
            Err(_) => (f64::NAN, f64::NAN, f64::NAN),
        };
        let (actual_dist_min, actual_dist_mean, actual_dist_max) = match summary_stats(&actual_dist)
        {
            Ok(x) => x,
            Err(_) => (f64::NAN, f64::NAN, f64::NAN),
        };
        let (actual_l_min, actual_l_mean, actual_l_max) =
            match summary_stats(&(actual_l.into_iter().map(|x| x as f64).collect::<Vec<f64>>())) {
                Ok(x) => x,
                Err(_) => (f64::NAN, f64::NAN, f64::NAN),
            };
        let (actual_k_min, actual_k_mean, actual_k_max) =
            match summary_stats(&(actual_k.into_iter().map(|x| x as f64).collect::<Vec<f64>>())) {
                Ok(x) => x,
                Err(_) => (f64::NAN, f64::NAN, f64::NAN),
            };
        println!(
            "Actual minimum correlation threshold statistics | min={}; mean={}; max={}",
            actual_corr_min, actual_corr_mean, actual_corr_max
        );
        println!(
            "Actual maximum distance threshold statistics | min={}; mean={}; max={}",
            actual_dist_min, actual_dist_mean, actual_dist_max
        );
        println!(
            "Actual minimum l-linked loci statistics | min={}; mean={}; max={}",
            actual_l_min, actual_l_mean, actual_l_max
        );
        println!(
            "Actual minimum k-nearest neighbours statistics | min={}; mean={}; max={}",
            actual_k_min, actual_k_mean, actual_k_max
        );
        Ok(self)
    }
}

// Impute using adaptive linkage disequilibrium (estimated using correlations within a window) k-nearest neighbour weighted allele frequencies imputation
pub fn impute_aldknni(
    mut genotypes_and_phenotypes: GenotypesAndPhenotypes,
    filter_stats: &FilterStats,
    window_size_bp: &u64,
    window_slide_size_bp: &u64,
    min_loci_per_window: &u64,
    min_loci_corr: &f64,
    max_pool_dist: &f64,
    optimise_for_thresholds: &bool,
    optimise_n_steps_corr: &usize,
    optimise_n_steps_dist: &usize,
    optimise_n_reps: &usize,
    do_linkimpute_weighted_mode: bool,
    misc_min_l: &u64,
    misc_min_k: &u64,
    n_threads: &usize,
    out: &String,
) -> io::Result<String> {

    // Will need to remove window-related parameters and rename misc_min_l and misc_min_l to be main parameters
    // Also note that the we are no longer reverting to MVI as setting and optimising for min_l and min_k clashes with that idea

    let (min_loci_corr, max_pool_dist, mae) = optimise_params_and_estimate_accuracy(
        &genotypes_and_phenotypes,
        min_loci_corr,
        max_pool_dist,
        optimise_for_thresholds,
        optimise_n_steps_corr,
        optimise_n_steps_dist,
        optimise_n_reps,
        window_size_bp,
        window_slide_size_bp,
        min_loci_per_window,
        do_linkimpute_weighted_mode,
        misc_min_l,
        misc_min_k,
    )
    .expect("Error calling optimise_params_and_estimate_accuracy() in impute_aldknni().");
    println!("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    if min_loci_corr > 1.0 {
        println!("Number of linked loci included = {}", min_loci_corr);
    } else {
        println!(
            "Minimum loci correlation considered in LD = {}",
            min_loci_corr
        );
    }
    if max_pool_dist > 1.0 {
        println!("Number of nearest neighbours = {}", max_pool_dist);
    } else {
        println!(
            "Maximum genetic distance between nearest neighbours = {}",
            max_pool_dist
        );
    }
    println!(
        "Expected imputation accuracy in terms of mean absolute error: {}",
        mae
    );
    println!("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .adaptive_ld_knn_imputation(
            window_size_bp,
            window_slide_size_bp,
            min_loci_per_window,
            &min_loci_corr,
            &max_pool_dist,
            true,
            do_linkimpute_weighted_mode,
            misc_min_l,
            misc_min_k,
        )
        .expect("Error calling adaptive_ld_knn_imputation() within impute_aldknni().");
    let end = std::time::SystemTime::now();
    let duration = end.duration_since(start).expect("Error measuring the duration of running adaptive_ld_knn_imputation() within impute_aldknni().");
    println!(
        "Adaptive LD-kNN imputation: {} pools x {} loci | Missingness: {}% | Duration: {} seconds",
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        genotypes_and_phenotypes.missing_rate().expect("Error measuring sparsity of the data using missing_rate() method within impute_aldknni()."),
        duration.as_secs()
    );
    // Remove 100% of the loci with missing data
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .filter_out_top_missing_loci(&1.00)
        .expect("Error calling filter_out_top_missing_loci() method within impute_aldknni().");
    let end = std::time::SystemTime::now();
    let duration = end.duration_since(start).expect("Error measuring the duration of running filter_out_top_missing_loci() within impute_aldknni().");
    println!(
        "Missing data removed, i.e. loci which cannot be imputed because of extreme sparsity: {} pools x {} loci | Missingness: {}% | Duration: {} seconds",
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        genotypes_and_phenotypes.missing_rate().expect("Error measuring sparsity of the data using missing_rate() method after filtering for missing top loci within impute_aldknni()."),
        duration.as_secs()
    );
    // Output
    let out = genotypes_and_phenotypes
        .write_csv(filter_stats, false, out, n_threads)
        .expect(
            "Error writing the output file using the write_csv() method within impute_aldknni().",
        );
    Ok(out)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_aldknni() {
        let file_sync = FileSync {
            filename: "./tests/test.sync".to_owned(),
            test: "load".to_owned(),
        };
        let file_phen = FilePhen {
            filename: "./tests/test_pheno.csv".to_owned(),
            delim: ",".to_owned(),
            names_column_id: 0,
            sizes_column_id: 1,
            trait_values_column_ids: vec![2, 3],
            format: "default".to_owned(),
        };
        let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
        let filter_stats = FilterStats {
            remove_ns: true,
            min_quality: 0.005,
            min_coverage: 1,
            min_allele_frequency: 0.005,
            max_missingness_rate: 0.0,
            pool_sizes: vec![20., 20., 20., 20., 20.],
        };
        let n_threads = 2;
        let mut frequencies_and_phenotypes = file_sync_phen
            .into_genotypes_and_phenotypes(&filter_stats, false, &n_threads)
            .unwrap();
        println!(
            "Before simulating missing data:\n{:?}",
            frequencies_and_phenotypes.intercept_and_allele_frequencies
        );
        let min_depth_below_which_are_missing = 5.0;
        let max_depth_above_which_are_missing = 100.0;
        let _ = frequencies_and_phenotypes
            .set_missing_by_depth(
                &min_depth_below_which_are_missing,
                &max_depth_above_which_are_missing,
            )
            .unwrap();
        // frequencies_and_phenotypes.intercept_and_allele_frequencies[(0, 1)] = f64::NAN;
        // frequencies_and_phenotypes.intercept_and_allele_frequencies[(0, 2)] = f64::NAN;
        // frequencies_and_phenotypes.intercept_and_allele_frequencies[(3, 5)] = f64::NAN;
        // frequencies_and_phenotypes.intercept_and_allele_frequencies[(3, 6)] = f64::NAN;
        println!(
            "Before imputation:\n{:?}",
            frequencies_and_phenotypes.intercept_and_allele_frequencies
        );
        let n_nan = frequencies_and_phenotypes
            .intercept_and_allele_frequencies
            .iter()
            .fold(0, |n_nan, &x| if x.is_nan() { n_nan + 1 } else { n_nan });
        println!("n_nan={}", n_nan);
        assert_eq!(n_nan, 22_997);

        // CORRELATION MATRIX CALCULATION
        let corr =
            calculate_genomewide_ld(&frequencies_and_phenotypes.intercept_and_allele_frequencies)
                .unwrap();
        assert_eq!(
            corr.len(),
            frequencies_and_phenotypes
                .intercept_and_allele_frequencies
                .ncols()
                - 1
        );
        assert_eq!(
            corr[0].len(),
            frequencies_and_phenotypes
                .intercept_and_allele_frequencies
                .ncols()
                - 1
        );
        assert_eq!(
            corr[1].len(),
            frequencies_and_phenotypes
                .intercept_and_allele_frequencies
                .ncols()
                - 2
        );
        assert_eq!(
            corr[2].len(),
            frequencies_and_phenotypes
                .intercept_and_allele_frequencies
                .ncols()
                - 3
        );
        // LINKED LOCI IDENTIFICATION
        let (linked_loci_idx, correlations) = find_l_linked_loci(7020, &corr, &0.99, 500).unwrap();
        // println!("linked_loci_idx={:?}", linked_loci_idx);
        // println!("linked_loci_idx[0]={:?}", linked_loci_idx[0]);
        // println!("frequencies_and_phenotypes.intercept_and_allele_frequencies.dim()={:?}", frequencies_and_phenotypes.intercept_and_allele_frequencies.dim());
        // println!("linked_loci_idx.len()={:?}", linked_loci_idx.len());
        assert_eq!(linked_loci_idx.len(), correlations.len());
        assert_eq!(linked_loci_idx.len(), 500);
        // GENETIC DISTANCES BETWEEN POOLS USING LINKED LOCI
        let idx_pool = 0;
        let distances = calculate_genetic_distances_between_pools(
            idx_pool,
            &linked_loci_idx,
            &frequencies_and_phenotypes.intercept_and_allele_frequencies,
        )
        .unwrap();
        // println!("distances={:?}", distances);
        assert_eq!(distances.len(), 5);
        // K-NEAREST NEIGHBOURS
        let idx_locus_ini = 1;
        let idx_locus_fin = 3;
        let (_idx_neighbours, distances, frequencies) = find_k_nearest_neighbours(
            &distances,
            &0.5,
            3,
            idx_locus_ini,
            idx_locus_fin,
            &frequencies_and_phenotypes.intercept_and_allele_frequencies,
        )
        .unwrap();
        // println!("idx_neighbours={:?}", idx_neighbours);
        // println!("distances={:?}", distances);
        // println!("frequencies={:?}", frequencies);
        assert_eq!(frequencies.dim(), (3, 2));
        // IMPUTATION
        let imputed_freq = impute_allele_frequencies(&frequencies, &distances, false).unwrap();
        println!("imputed_freq={:?}", imputed_freq);
        assert_eq!(imputed_freq.len(), 2);
        assert_eq!(imputed_freq[0] + imputed_freq[1], 1.0);

        let _frac_top_missing_pools = 0.0;
        let _frac_top_missing_loci = 0.2;
        let window_size_bp = 1e6 as u64;
        let window_slide_size_bp = window_size_bp;
        let min_loci_per_window = 1;
        let min_loci_corr = 0.9;
        let max_pool_dist = 0.1;
        let _optimise_for_thresholds = true;
        let _optimise_n_steps_corr = 10;
        let _optimise_n_steps_dist = 10;
        let _optimise_n_reps = 3;
        let _do_linkimpute_weighted_mode = false;
        let misc_min_l = 1;
        let misc_min_k = 1;
        let _ = frequencies_and_phenotypes
            .adaptive_ld_knn_imputation(
                &window_size_bp,
                &window_slide_size_bp,
                &min_loci_per_window,
                &min_loci_corr,
                &max_pool_dist,
                true,
                false,
                &misc_min_l,
                &misc_min_k,
            )
            .unwrap();
        println!(
            "After imputation:\n{:?}",
            frequencies_and_phenotypes.intercept_and_allele_frequencies
        );
        let n_nan = frequencies_and_phenotypes
            .intercept_and_allele_frequencies
            .iter()
            .fold(0, |n_nan, &x| if x.is_nan() { n_nan + 1 } else { n_nan });
        println!("n_nan={}", n_nan);
        assert_eq!(n_nan, 7_945); // corresponds to the 1_589 alleles completely missing across all pools

        let keep_p_minus_1 = false;
        let _start = std::time::SystemTime::now();
        let _genotypes_and_phenotypes = file_sync_phen
            .into_genotypes_and_phenotypes(&filter_stats, keep_p_minus_1, &n_threads)
            .unwrap();

        // let outname = impute_aldknni(
        //     genotypes_and_phenotypes,
        //     &filter_stats,
        //     &window_size_bp,
        //     &window_slide_size_bp,
        //     &min_loci_per_window,
        //     &min_loci_corr,
        //     &max_pool_dist,
        //     &optimise_for_thresholds,
        //     &optimise_n_steps_corr,
        //     &optimise_n_steps_dist,
        //     &optimise_n_reps,
        //     do_linkimpute_weighted_mode,
        //     &misc_min_l,
        //     &misc_min_k,
        //     &n_threads,
        //     &"test-impute_aldknni.csv".to_owned(),
        // )
        // .unwrap();
        // assert_eq!(outname, "test-impute_aldknni.csv".to_owned()); // Do better!!! Load data - thus working on improving load_table()

        // println!("frequencies_and_phenotypes.intercept_and_allele_frequencies.slice(s![0..5, 39..42])={:?}", frequencies_and_phenotypes.intercept_and_allele_frequencies.slice(s![0..5, 39..42]));
        // assert_eq!(
        //     frequencies_and_phenotypes
        //         .intercept_and_allele_frequencies
        //         .slice(s![0..5, 1..3])
        //         .sum_axis(Axis(1))
        //         .map(|x| sensible_round(*x, 2)),
        //     Array1::from_elem(5, 1.0)
        // );
        // assert_eq!(
        //     frequencies_and_phenotypes
        //         .intercept_and_allele_frequencies
        //         .slice(s![0..5, 39..42])
        //         .sum_axis(Axis(1))
        //         .map(|x| sensible_round(*x, 2)),
        //     Array1::from_elem(5, 1.0)
        // );
        // assert_eq!(
        //     frequencies_and_phenotypes
        //         .intercept_and_allele_frequencies
        //         .slice(s![0..5, 119..121])
        //         .sum_axis(Axis(1))
        //         .map(|x| sensible_round(*x, 2)),
        //     Array1::from_elem(5, 1.0)
        // );
        // assert_eq!(
        //     frequencies_and_phenotypes
        //         .intercept_and_allele_frequencies
        //         .slice(s![0..5, 400..402])
        //         .sum_axis(Axis(1))
        //         .map(|x| sensible_round(*x, 2)),
        //     Array1::from_elem(5, 1.0)
        // );
        // assert_eq!(0, 1);
    }
}
