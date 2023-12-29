use ndarray::{prelude::*, Zip};
use std::io;

use crate::helpers::*;
use crate::optim::*;
use crate::structs_and_traits::*;

// // For benchmarking
// use std::time::Instant;

fn calculate_mean_absolute_distances(
    window_freqs: ArrayView2<f64>,
    corr: ArrayView1<f64>,
    min_loci_corr: &f64,
) -> io::Result<(Array2<f64>, Vec<usize>)> {
    let (n, p) = window_freqs.dim();
    assert!(p == corr.len());
    assert!(*min_loci_corr >= 0.0);
    // If the correlation between the locus requiring imputation and other loci in the window is less than the minimum loci correlation or number of requested loci,
    // then we include all the loci to estimate genetic distance.
    let idx_linked_alleles = if corr
        .iter()
        .filter(|&x| !x.is_nan())
        .collect::<Vec<&f64>>()
        .len()
        > *min_loci_corr as usize
    {
        // If min_loci_corr > 1.0, then we assume this specifies the number of loci to include and not the minimum correlation
        let min_loci_corr = if *min_loci_corr > 1.0 {
            let mut vec_corr = corr
                .to_owned()
                .into_iter()
                .filter(|&x| !x.is_nan())
                .collect::<Vec<f64>>();
            // println!("vec_corr={:?}", vec_corr);
            // println!("min_loci_corr={:?}", min_loci_corr);
            vec_corr.sort_by(|a, b| b.partial_cmp(a).unwrap());
            let m = if *min_loci_corr as usize > vec_corr.len() {
                vec_corr.len()
            } else {
                *min_loci_corr as usize
            };
            vec_corr[m - 1]
        } else {
            *min_loci_corr
        };
        // Identify loci to be used for distance estimation, i.e. with correlations at least min_loci_corr
        let mut idx_linked_alleles: Vec<usize> = vec![];
        for i in 0..p {
            if corr[i] >= min_loci_corr {
                idx_linked_alleles.push(i);
            }
        }
        idx_linked_alleles
    } else {
        (0..p).collect::<Vec<usize>>()
    };
    // Estimate pairwise distances between pools using mean absolute differences
    let mut dist: Array2<f64> = Array2::from_elem((n, n), 1.0); // Start with equally farthest distances across all pool pairs
    let window_freqs_linked_alleles = window_freqs.select(Axis(1), &idx_linked_alleles);
    for i in 0..n {
        let pool_i = window_freqs_linked_alleles.row(i);
        for j in i..n {
            let pool_j = window_freqs_linked_alleles.row(j);
            let mut d = 0.0;
            let mut m = 0.0;
            for k in 0..pool_i.len() {
                if !pool_i[k].is_nan() & !pool_j[k].is_nan() {
                    d += (pool_i[k] - pool_j[k]).abs();
                    m += 1.00;
                }
            }
            dist[(i, j)] = if (d == 0.0) & (m == 0.0) { 1.0 } else { d / m };
            dist[(j, i)] = dist[(i, j)];
        }
    }
    Ok((dist, idx_linked_alleles))
}

fn find_k_nearest_neighbours(
    window_freqs_col: ArrayView1<f64>,
    dist: ArrayView1<f64>,
    max_pool_dist: &f64,
) -> io::Result<(Array1<f64>, Array1<f64>)> {
    let n = window_freqs_col.len();
    assert!(n == dist.len());
    // if max_pool_dist > 1.0, then we assume this specifies the number of pools to include and not the minimum distance
    let max_pool_dist = if *max_pool_dist > 1.0 {
        // The distance matrix do not have missing values, if the distance cannot be calculated then the default distance is set to the maximum which is 1.00.
        // let mut vec_dist = dist.to_owned().into_iter().filter(|&x| !x.is_nan()).collect::<Vec<f64>>();
        let mut vec_dist = dist.to_vec();
        // println!("dist={:?}", dist);
        // println!("0.0/0.0={:?}", 0.0/0.0);
        vec_dist.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let m = if *max_pool_dist as usize > vec_dist.len() {
            vec_dist.len()
        } else {
            *max_pool_dist as usize
        };
        // println!("dist={:?}", dist);
        // println!("vec_dist={:?}", vec_dist);
        // println!("m={:?}", m);
        // println!("vec_dist[m - 1]={:?}", vec_dist[m - 1]);
        vec_dist[m - 1]
    } else {
        *max_pool_dist
    };
    // println!("max_pool_dist={:?}", max_pool_dist);
    // Identiy the k-nearest neighbours, i.e. with distance at most max_pool_dist
    let mut freqs_k_neighbours: Vec<f64> = vec![];
    let mut dist_k_neighbours: Vec<f64> = vec![];
    for i in 0..n {
        if dist[i] <= max_pool_dist {
            if !window_freqs_col[i].is_nan() {
                freqs_k_neighbours.push(window_freqs_col[i]);
                dist_k_neighbours.push(dist[i]);
            }
        }
    }
    Ok((
        Array1::from_vec(freqs_k_neighbours),
        Array1::from_vec(dist_k_neighbours),
    ))
}

fn correct_allele_frequencies_per_locus<'w>(
    a: usize,
    i: usize,
    j: usize,
    idx_ini: usize,
    window_freqs: &'w mut Array2<f64>,
    idx_window_head: &Vec<usize>,
    idx_window_tail: &Vec<usize>,
    loci_idx: &Vec<usize>,
) -> io::Result<&'w mut Array2<f64>> {
    // println!("@@@@@@@@@@@@@@@@@@@correct_allele_frequencies_per_locus@@@@@@@@@@@@@@@@@@@");
    // let start = Instant::now();
    // Include the start of the next window, i.e. the marker for the end of the last locus in the current window
    let loci_start_indexes_within_the_current_window =
        loci_idx[idx_window_head[a]..(idx_window_tail[a] + 2)].to_vec();
    for j_ in 1..loci_start_indexes_within_the_current_window.len() {
        // Are we at the last allele of the locus?
        if (loci_start_indexes_within_the_current_window[j_] - 1) == (idx_ini + j) {
            // If we are then we find the start of this locus, i.e. its local index
            let j_ini = loci_start_indexes_within_the_current_window[j_ - 1] - idx_ini;
            let freqs_sum = window_freqs
                .slice(s![i, j_ini..(j + 1)])
                .fold(0.0, |sum, &x| if !x.is_nan() { sum + x } else { sum });
            if freqs_sum != 1.0 {
                for j_ in j_ini..(j + 1) {
                    window_freqs[(i, j_)] = if freqs_sum == 0.0 {
                        window_freqs[(i, j_)]
                    } else {
                        window_freqs[(i, j_)] / freqs_sum
                    };
                }
            }
            break;
        }
    }
    // let duration = start.elapsed();
    // println!("duration = {:?}", duration);
    // println!("window_freqs = {:?}", window_freqs);
    Ok(window_freqs)
}

impl GenotypesAndPhenotypes {
    pub fn adaptive_ld_knn_imputation(
        &mut self,
        window_size_bp: &u64,
        window_slide_size_bp: &u64,
        min_loci_per_window: &u64,
        min_loci_corr: &f64,
        max_pool_dist: &f64,
        show_mvi_revert_stats: bool
    ) -> io::Result<&mut Self> {
        self.check().unwrap();
        // We are assuming that all non-zero alleles across pools are kept, i.e. biallelic loci have 2 columns, triallelic have 3, and so on.
        let (n, _p) = self.intercept_and_allele_frequencies.dim();
        // Define sliding windows
        let window_size_bp = if *window_size_bp == 0 {
            // If window size is not set, then use the length of the largest chromosome or scaffold
            self.position
                .iter()
                .fold(0, |max_len, &x| if x > max_len { x } else { max_len })
        } else {
            *window_size_bp
        };
        let (loci_idx, loci_chr, loci_pos) = self.count_loci().unwrap();
        let mut loci_chr_no_redundant_tail = loci_chr.to_owned();
        loci_chr_no_redundant_tail.pop();
        let mut loci_pos_no_redundant_tail = loci_pos.to_owned();
        loci_pos_no_redundant_tail.pop();
        let (idx_window_head, idx_window_tail) = define_sliding_windows(
            &loci_chr_no_redundant_tail,
            &loci_pos_no_redundant_tail,
            &window_size_bp,
            window_slide_size_bp,
            min_loci_per_window,
        )
        .unwrap();
        let w = idx_window_head.len();
        // Parallel processing per window
        let mut vec_windows_freqs: Vec<Array2<f64>> = vec![Array2::from_elem((1, 1), f64::NAN); w];
        let idx_windows: Vec<usize> = (0..w).collect();
        let mut params_ld_k: Vec<Vec<ParametersOfLinkedLociAndNearestNeighbours>> = vec![vec![]; w];
        let mut aldknni_instances: Vec<f64> = vec![0.0; w];
        let mut mvi_instances: Vec<f64> = vec![0.0; w];
        Zip::from(&mut vec_windows_freqs)
            .and(&idx_windows)
            .and(&mut params_ld_k)
            .and(&mut aldknni_instances)
            .and(&mut mvi_instances)
            .par_for_each(
                |window_freqs, &a, params, n_aldknni, n_mvi| {
                    // for a in 0..idx_windows.len() {
                    let idx_ini = loci_idx[idx_window_head[a]];
                    let idx_fin = loci_idx[idx_window_tail[a] + 1]; // add one so that we include the last part of the window!
                    let p = idx_fin - idx_ini;
                    *window_freqs = self
                        // let mut window_freqs = self
                        .intercept_and_allele_frequencies
                        .slice(s![.., idx_ini..idx_fin])
                        .to_owned();
                    // Vector of loci and neighbour parameters
                    *params = vec![];
                    // Calculate correlations between alleles across all loci within the window (calculate only for the upper tringular and just mirror the results into the lower tringular for efficiency)
                    let mut corr: Array2<f64> = Array2::from_elem((p, p), f64::NAN);
                    for j0 in 0..p {
                        for j1 in j0..p {
                            // Note that we are perfectly fine using pairwise complete observations to compute correlations between loci as we are interested in the most correlated loci per focul locus, and not in comparing correlations between pairs of loci
                            match pearsons_correlation_pairwise_complete(
                                &window_freqs.column(j0),
                                &window_freqs.column(j1),
                            ) {
                                Ok(x) => {
                                    corr[(j0, j1)] = x.0;
                                    corr[(j1, j0)] = x.0;
                                }
                                Err(_) => continue,
                            };
                        }
                    }
                    for j in 0..p {
                        let locus_sparsity =
                            window_freqs
                                .column(j)
                                .fold(0, |sum, &x| if x.is_nan() { sum + 1 } else { sum });
                        if locus_sparsity == 0 {
                            // No missing data on the locus
                            continue;
                        } else {
                            let (dist, idx_linked_alleles) = calculate_mean_absolute_distances(
                                window_freqs.view(),
                                corr.column(j),
                                min_loci_corr,
                            )
                            .unwrap();
                            // Now, let's find the pools needing imputation and impute them using the k-nearest neighbours
                            for i in 0..n {
                                if !window_freqs[(i, j)].is_nan() {
                                    continue;
                                } else {
                                    let (freqs_k_neighbours, dist_k_neighbours) =
                                        find_k_nearest_neighbours(
                                            window_freqs.column(j),
                                            dist.column(i),
                                            max_pool_dist,
                                        )
                                        .unwrap();
                                    // Impute
                                    // println!("dist={:?}", dist);
                                    // println!("window_freqs.column(j)={:?}", window_freqs.column(j));
                                    // println!("freqs_k_neighbours={:?}", freqs_k_neighbours);
                                    // println!("dist_k_neighbours={:?}", dist_k_neighbours);
                                    params.push(ParametersOfLinkedLociAndNearestNeighbours {
                                        min_loci_corr: corr
                                            .column(j)
                                            .select(Axis(0), &idx_linked_alleles)
                                            .iter()
                                            .fold(1.0, |min, &x| if x < min { x } else { min }),
                                        max_pool_dist: dist_k_neighbours
                                            .iter()
                                            .fold(0.0, |max, &x| if x > max { x } else { max }),
                                        l_linked_loci: idx_linked_alleles.len() as f64,
                                        k_nearest_neighbours: freqs_k_neighbours.len() as f64,
                                    });
                                    if freqs_k_neighbours.is_empty() {
                                        window_freqs[(i, j)] =
                                            mean_array1_ignore_nan(&window_freqs.column(j))
                                                .unwrap();
                                        *n_mvi = *n_mvi + 1.0;
                                    } else {
                                        let weights = 1.00 - &dist_k_neighbours;
                                        let weights_sum = weights.iter().fold(0.0, |sum, x| sum + x);
                                        let weights = if weights_sum==0.0 {
                                            weights
                                        } else {
                                            weights / weights_sum
                                        };
                                        window_freqs[(i, j)] = (&freqs_k_neighbours * &weights).sum();
                                        *n_aldknni = *n_aldknni + 1.0;
                                        // if window_freqs[(i, j)].is_nan() {
                                        //     println!("corr.column(j).select(Axis(0), &idx_linked_alleles)={:?}", corr.column(j).select(Axis(0), &idx_linked_alleles));
                                        //     println!("dist_k_neighbours={:?}", dist_k_neighbours);
                                        //     println!("weights={:?}", weights);
                                        //     println!("freqs_k_neighbours={:?}", freqs_k_neighbours);
                                        // }
                                    };
                                }
                                // Need to correct for when the imputed allele frequencies do not add up to one!
                                if j > 0 {
                                    correct_allele_frequencies_per_locus(
                                        a,
                                        i,
                                        j,
                                        idx_ini,
                                        window_freqs,
                                        // &mut window_freqs,
                                        &idx_window_head,
                                        &idx_window_tail,
                                        &loci_idx,
                                    )
                                    .unwrap();
                                }
                            } // Impute across pools with missing data
                        } // Impute if we have missing data
                    } // Iterate across alleles across loci within the window
                }, // Parallel processing across windows
            );
        // }
        // println!("params_ld_k:\n{:?}", params_ld_k);
        // println!("vec_windows_freqs[0]:\n{:?}", vec_windows_freqs[0]);
        // Write-out the imputed data
        // println!("@@@@@@@@@@@@@@@@@@@Writing imputed data@@@@@@@@@@@@@@@@@@@");
        // let start = Instant::now();
        let mut vec_corr: Vec<f64> = vec![];
        let mut vec_dist: Vec<f64> = vec![];
        let mut vec_l: Vec<f64> = vec![];
        let mut vec_k: Vec<f64> = vec![];
        let mut n_aldknni_instances = 0.0;
        let mut n_mvi_instances = 0.0;
        for a in 0..w {
            // Use the indexes of each locus
            let idx_ini = loci_idx[idx_window_head[a]];
            let idx_fin = loci_idx[idx_window_tail[a] + 1]; // add one so that we include the last part of the window!
            let p = idx_fin - idx_ini;
            for i in 0..n {
                for j in 0..p {
                    self.intercept_and_allele_frequencies[(i, idx_ini + j)] =
                        vec_windows_freqs[a][(i, j)];
                }
            }
            let par = params_ld_k[a].clone();
            for i in 0..par.len() {
                vec_corr.push(par[i].min_loci_corr);
                vec_dist.push(par[i].max_pool_dist);
                vec_l.push(par[i].l_linked_loci);
                vec_k.push(par[i].k_nearest_neighbours);
            }
            n_aldknni_instances = n_aldknni_instances + aldknni_instances[a];
            n_mvi_instances = n_mvi_instances + mvi_instances[a];
        }
        let corr_mean = vec_corr.iter().fold(0.0, |sum, &x| sum + x) / (vec_corr.len() as f64);
        let dist_mean = vec_dist.iter().fold(0.0, |sum, &x| sum + x) / (vec_dist.len() as f64);
        let l_mean = vec_l.iter().fold(0.0, |sum, &x| sum + x) / (vec_l.len() as f64);
        let l_ge5 = vec_l.iter().fold(0, |n_greater_than_or_equal_5, &x| {
            if x >= 5.0 {
                n_greater_than_or_equal_5 + 1
            } else {
                n_greater_than_or_equal_5
            }
        });
        let k_mean = vec_k.iter().fold(0.0, |sum, &x| sum + x) / (vec_k.len() as f64);
        let k_ge5 = vec_k.iter().fold(0, |n_greater_than_or_equal_5, &x| {
            if x >= 5.0 {
                n_greater_than_or_equal_5 + 1
            } else {
                n_greater_than_or_equal_5
            }
        });
        if show_mvi_revert_stats {
            println!("Missing sites reverted to mean value imputation = {}% ( {} aldknni instances; {} mvi instances)", n_mvi_instances * 100.0 / (n_aldknni_instances + n_mvi_instances), n_aldknni_instances, n_mvi_instances);
            println!("Mean minimum correlation for loci to be considered in LD = {}; Mean maximum genetic distance for pools to be considered nearest neighbours = {}", corr_mean, dist_mean);
            println!("Mean number of loci included in the genetic distance estimation = {}; Number of loci greater than 5 included in the genetic distance estimation = {}", l_mean, l_ge5);
            println!("Mean number of neighbours included in the imputation = {} (if zero, then MVI is generally used); Number of neighbours greater than 5 included in the imputation = {}", k_mean, k_ge5);
        }
        // let duration = start.elapsed();
        // println!("duration = {:?}", duration);
        // Set missing coverages to infinity to mark imputed data
        // println!("@@@@@@@@@@@@@@@@@@@Set missing coverages to infinity to mark imputed data@@@@@@@@@@@@@@@@@@@");
        // let start = Instant::now();
        for j in 0..self.coverages.ncols() {
            // Mark only the imputed loci, i.e. loci which were not completely missing across all pools
            let n_non_missing = self.coverages.select(Axis(1), &[j]).fold(0, |sum, &x| {
                if !x.is_nan() {
                    sum + 1
                } else {
                    sum
                }
            });
            if n_non_missing > 0 {
                for i in 0..self.coverages.nrows() {
                    if self.coverages[(i, j)].is_nan() {
                        self.coverages[(i, j)] = f64::INFINITY
                    };
                }
            }
        }
        // let duration = start.elapsed();
        // println!("duration = {:?}", duration);
        // println!("self.coverages = {:?}", self.coverages);
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
    n_threads: &usize,
    out: &String,
) -> io::Result<String> {
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
    )
    .unwrap();
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
            true
        )
        .unwrap();
    let end = std::time::SystemTime::now();
    let duration = end.duration_since(start).unwrap();
    println!(
        "Adaptive LD-kNN imputation: {} pools x {} loci | Missingness: {}% | Duration: {} seconds",
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        genotypes_and_phenotypes.missing_rate().unwrap(),
        duration.as_secs()
    );
    // Remove 100% of the loci with missing data
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .filter_out_top_missing_loci(&1.00)
        .unwrap();
    let end = std::time::SystemTime::now();
    let duration = end.duration_since(start).unwrap();
    println!(
        "Missing data removed, i.e. loci which cannot be imputed because of extreme sparsity: {} pools x {} loci | Missingness: {}% | Duration: {} seconds",
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        genotypes_and_phenotypes.missing_rate().unwrap(),
        duration.as_secs()
    );
    // Output
    let out = genotypes_and_phenotypes
        .write_csv(filter_stats, false, out, n_threads)
        .unwrap();
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
            filename: "./tests/test.csv".to_owned(),
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
        let _frac_top_missing_pools = 0.0;
        let _frac_top_missing_loci = 0.2;
        let window_size_bp = 1e6 as u64;
        let window_slide_size_bp = window_size_bp;
        let min_loci_per_window = 1;
        let min_loci_corr = 0.9;
        let max_pool_dist = 0.1;
        let optimise_for_thresholds = true;
        let optimise_n_steps_corr = 10;
        let optimise_n_steps_dist = 10;
        let optimise_n_reps = 3;
        let _ = frequencies_and_phenotypes
            .adaptive_ld_knn_imputation(
                &window_size_bp,
                &window_slide_size_bp,
                &min_loci_per_window,
                &min_loci_corr,
                &max_pool_dist,
                true
            )
            .unwrap();
        println!(
            "After imputation:\n{:?}",
            frequencies_and_phenotypes.intercept_and_allele_frequencies
        );
        let keep_p_minus_1 = false;
        let _start = std::time::SystemTime::now();
        let genotypes_and_phenotypes = file_sync_phen
            .into_genotypes_and_phenotypes(&filter_stats, keep_p_minus_1, &n_threads)
            .unwrap();

        let outname = impute_aldknni(
            genotypes_and_phenotypes,
            &filter_stats,
            &window_size_bp,
            &window_slide_size_bp,
            &min_loci_per_window,
            &min_loci_corr,
            &max_pool_dist,
            &optimise_for_thresholds,
            &optimise_n_steps_corr,
            &optimise_n_steps_dist,
            &optimise_n_reps,
            &n_threads,
            &"test-impute_aldknni.csv".to_owned(),
        )
        .unwrap();
        assert_eq!(outname, "test-impute_aldknni.csv".to_owned()); // Do better!!! Load data - thus working on improving load_table()

        println!("frequencies_and_phenotypes.intercept_and_allele_frequencies.slice(s![0..5, 39..42])={:?}", frequencies_and_phenotypes.intercept_and_allele_frequencies.slice(s![0..5, 39..42]));
        assert_eq!(
            frequencies_and_phenotypes
                .intercept_and_allele_frequencies
                .slice(s![0..5, 1..3])
                .sum_axis(Axis(1))
                .map(|x| sensible_round(*x, 2)),
            Array1::from_elem(5, 1.0)
        );
        assert_eq!(
            frequencies_and_phenotypes
                .intercept_and_allele_frequencies
                .slice(s![0..5, 39..42])
                .sum_axis(Axis(1))
                .map(|x| sensible_round(*x, 2)),
            Array1::from_elem(5, 1.0)
        );
        assert_eq!(
            frequencies_and_phenotypes
                .intercept_and_allele_frequencies
                .slice(s![0..5, 119..121])
                .sum_axis(Axis(1))
                .map(|x| sensible_round(*x, 2)),
            Array1::from_elem(5, 1.0)
        );
        assert_eq!(
            frequencies_and_phenotypes
                .intercept_and_allele_frequencies
                .slice(s![0..5, 400..402])
                .sum_axis(Axis(1))
                .map(|x| sensible_round(*x, 2)),
            Array1::from_elem(5, 1.0)
        );
        // assert_eq!(0, 1);
    }
}
