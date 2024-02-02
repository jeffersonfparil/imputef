use ndarray::{prelude::*, Zip};
use std::cmp::Ordering;
use std::io;

use crate::helpers::*;
use crate::optim::*;
use crate::structs_and_traits::*;

pub fn calculate_genomewide_ld(
    genotypes_and_phenotypes: &GenotypesAndPhenotypes,
    restrict_linked_loci_per_chromosome: bool,
) -> io::Result<Vec<Vec<u8>>> {
    // Errors and f64::NAN are all converted into 0.0 for simplicity
    let (_n, p) = genotypes_and_phenotypes
        .intercept_and_allele_frequencies
        .dim();
    let mut corr: Vec<Vec<u8>> = vec![vec![]; p - 1];
    Zip::indexed(&mut corr).par_for_each(|idx, c| {
        let current_chromosome = genotypes_and_phenotypes.chromosome[idx].clone();
        for j in (idx + 1)..p {
            if restrict_linked_loci_per_chromosome
                & (genotypes_and_phenotypes.chromosome[j] != current_chromosome)
            {
                // We are breaking because we assume that the loci are sorted by chromosomes.
                break;
            } else {
                let corr = match pearsons_correlation_pairwise_complete(
                    &genotypes_and_phenotypes
                        .intercept_and_allele_frequencies
                        .column(idx),
                    &genotypes_and_phenotypes
                        .intercept_and_allele_frequencies
                        .column(j),
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
                c.push((corr * 255.0).round() as u8);
            }
        }
    });
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
    let mut distances: Array1<f64> = Array1::from_elem(n, 1.0);
    // for i in 0..n {
    //     if i != idx_row {
    //         let q1: Array1<f64> = intercept_and_allele_frequencies
    //             .row(i)
    //             .select(Axis(0), linked_loci_idx);
    //         distances[i] = match (&q - &q1)
    //             .into_iter()
    //             .filter(|&x| !x.is_nan())
    //             .collect::<Array1<f64>>()
    //             .map(|x| x.abs())
    //             .mean()
    //         {
    //             Some(x) => x,
    //             None => 1.00,
    //         };
    //     }
    // }
    Zip::indexed(&mut distances).par_for_each(|i, d| {
        if i != idx_row {
            let q1: Array1<f64> = intercept_and_allele_frequencies
                .row(i)
                .select(Axis(0), linked_loci_idx);
            *d = (&q - &q1)
                .into_iter()
                .filter(|&x| !x.is_nan())
                .collect::<Array1<f64>>()
                .map(|x| x.abs())
                .mean()
                .unwrap_or(1.00);
        }
    });
    Ok(distances)
}

fn find_l_linked_loci(
    idx_col: usize,
    corr: &Vec<Vec<u8>>,
    min_loci_corr: &f64,
    min_l_loci: usize,
    restrict_linked_loci_per_chromosome: bool,
    current_chromosome: &String,
    chromosomes: &Vec<String>,
) -> io::Result<(Vec<usize>, Vec<f64>)> {
    // Find loci passing the min_loci_corr threshold, then sort and add more loci if it did not pass the min_l_loci threshold
    //      - to optimise for the thresholds, set min_l_loci to 0
    //      - to optimise for counts, set min_loci_corr to 0.0
    assert!(
        min_l_loci > 0,
        "Error: the minimum number of linked loci need to be greater than 0."
    );
    let p = corr.len();
    let mut vec_corr: Vec<f64> = vec![];
    let mut vec_idx: Vec<usize> = vec![];
    // Across rows and at the idx_col of the triangular matrix of correlations
    for i in 0..idx_col {
        if restrict_linked_loci_per_chromosome & (chromosomes[i] != *current_chromosome) {
            // Break if we are not on the same chromosome and we are restricting linked loci to be only within chromosomes. Furthermore, we are assuming loci are sorted by chromosome and position.
            break;
        }
        let c = (corr[i][idx_col - (i + 1)] as f64) / 255.0; // Less the current index as the length of the vectors decreases by 1 each time in addition to 1 less element representing the diagonal elements, i.e. itself
        if c >= *min_loci_corr {
            vec_idx.push(i);
        }
        vec_corr.push(c);
    }
    // Across columns and at the idx_col row of the triangular matrix of correlations
    for j in 0..(p - idx_col) {
        if restrict_linked_loci_per_chromosome
            & (chromosomes[idx_col + j + 1] != *current_chromosome)
        {
            // Break if we are not on the same chromosome and we are restricting linked loci to be only within chromosomes. Furthermore, we are assuming loci are sorted by chromosome and position.
            // See comment below for the explanation for chromosomes[idx_col + j + 1]
            break;
        }
        let c = (corr[idx_col][j] as f64) / 255.0;
        if c >= *min_loci_corr {
            vec_idx.push(idx_col + j + 1); // Add the current index as the length of the vectors decreases by 1 each time plus skipping the locus requiring imputation as we did not calculate the correlation of each locus with itself
        }
        vec_corr.push(c);
    }
    // If less than the minimum number of loci passed the threshold,
    // or if we are not filtering by minimum loci correlation (i.e. min_loci_corr = 0.0),
    // then we sort the correlations (decreasing) and pick the top-most correlated loci
    let (vec_idx, vec_corr) = if (vec_idx.len() < min_l_loci) | (*min_loci_corr == 0.0) {
        // Extract indices to reach min_l_loci or if we do not have enough loci, then just vec_corr.len()
        let mut indices: Vec<usize> = (0..vec_corr.len()).collect();
        indices.sort_by(
            |&a, &b| match (vec_corr[b].is_nan(), vec_corr[a].is_nan()) {
                (true, true) => Ordering::Equal,
                (true, false) => Ordering::Greater,
                (false, true) => Ordering::Less,
                (false, false) => vec_corr[b].partial_cmp(&vec_corr[a]).unwrap(),
            },
        );
        let l_linked_loci = if vec_corr.len() < min_l_loci {
            vec_corr.len()
        } else {
            min_l_loci
        };
        let indices = indices[0..l_linked_loci].to_vec();
        // Extract the correlations corresponding to the extracted indices
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
    min_k_neighbours: usize,
    idx_col: usize,
    intercept_and_allele_frequencies: &Array2<f64>,
) -> io::Result<(Vec<f64>, Array1<f64>)> {
    // Find neighbours passing the max_pool_dist threshold, then sort and add more neighbours if it did not pass the min_k_neighbours threshold
    //      - to optimise for the thresholds, set min_k_neighbours to 0
    //      - to optimise for counts, set max_pool_dist to 1.0 (maximum possible distance, i.e. mean absolute difference of 1.0 in allele frequency)
    assert!(
        min_k_neighbours > 0,
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
    // Are we filtering by maximum distance?
    if *max_pool_dist < 1.00 {
        for i in indices.clone().into_iter() {
            let bool_dist_nan = !distances[i].is_nan();
            let bool_dist_thresh = distances[i] <= *max_pool_dist;
            let bool_freq_nan = !intercept_and_allele_frequencies[(i, idx_col)].is_nan();
            if bool_dist_nan & bool_dist_thresh & bool_freq_nan {
                idx.push(i);
                dist.push(distances[i]);
            }
        }
    }
    // If less than the minimum number of neighbours passed the threshold,
    // or if we are not filtering by maximum pool distance (i.e. max_pool_dist = 1.00),
    // then we sort the distances (increasing) and pick the nearest neighbours
    if (idx.len() < min_k_neighbours) | (*max_pool_dist == 1.00) {
        idx = vec![];
        dist = vec![];
        for i in indices.clone().into_iter() {
            let bool_dist_nan = !distances[i].is_nan();
            let bool_freq_nan = !intercept_and_allele_frequencies[(i, idx_col)].is_nan();
            if idx.len() == min_k_neighbours {
                break;
            }
            if bool_dist_nan & bool_freq_nan {
                idx.push(i);
                dist.push(distances[i]);
            }
        }
    }
    // Extract non-missing allele frequencies at the locus requiring imputation of the k-nearest neighbours
    let mut freq: Array1<f64> = Array1::from_elem(idx.len(), f64::NAN);
    for i in 0..idx.len() {
        freq[i] = intercept_and_allele_frequencies[(idx[i], idx_col)];
    }
    Ok((dist, freq))
}

fn impute_allele_frequencies(
    frequencies: &Array1<f64>,
    // frequencies: &Array2<f64>,
    distances: &Vec<f64>,
) -> io::Result<f64> {
    // ) -> io::Result<Vec<f64>> {
    let n = frequencies.len();
    for i in 0..n {
        assert!(
            !frequencies[i].is_nan(),
            "Error: There should be no missing allele frequencies here."
        );
    }
    let mut imputed_freqs = 0.0;
    // Perform weighted mean allele frequencies
    let additive_inverse_plus_epsilon: Array1<f64> =
        Array1::from_iter(distances.iter().map(|&x| (1.0 - x) + f64::EPSILON)); // Add a small value if all distances are 1.0 and allow equal contributions in such cases
    let weights_sum: f64 = additive_inverse_plus_epsilon
        .iter()
        .fold(0.0, |sum, x| sum + x);
    let weights: Array1<f64> = additive_inverse_plus_epsilon / weights_sum;
    for i in 0..n {
        imputed_freqs += weights[i] * frequencies[i];
    }
    if imputed_freqs.is_nan() {
        println!("frequencies={:?}", frequencies);
        println!("distances={:?}", distances);
    }
    Ok(imputed_freqs)
}

impl GenotypesAndPhenotypes {
    pub fn adaptive_ld_knn_imputation(
        &mut self,
        self_clone: &GenotypesAndPhenotypes,
        min_loci_corr: &f64,
        max_pool_dist: &f64,
        min_l_loci: &u64,
        min_k_neighbours: &u64,
        loci_idx: &Vec<usize>,
        corr: &Vec<Vec<u8>>,
        restrict_linked_loci_per_chromosome: bool,
    ) -> io::Result<&mut Self> {
        self.check().expect("Error self.check() within adaptive_ld_knn_imputation() method of GenotypesAndPhenotypes struct.");
        // We are assuming that all non-zero alleles across pools are kept, i.e. biallelic loci have 2 columns, triallelic have 3, and so on.
        let (n, p) = self.intercept_and_allele_frequencies.dim();
        let min_l_loci = if *min_l_loci >= (p as u64) {
            p - 1
        } else {
            *min_l_loci as usize
        };
        let min_k_neighbours = if *min_k_neighbours >= (n as u64) {
            n - 1
        } else {
            *min_k_neighbours as usize
        };
        Zip::indexed(&mut self.intercept_and_allele_frequencies)
        .par_for_each(|(i, j), q| {
            if q.is_nan() {
                let current_chromosome = self_clone.chromosome[j].to_owned();
                // Find loci most correlated to the major allele of the current locus, i.e. the first allele of the locus as they were sorted by decreasing allele frequency (see Sort trait)
                let (linked_loci_idx, _correlations) =
                    find_l_linked_loci(j, corr, min_loci_corr, min_l_loci,
                        restrict_linked_loci_per_chromosome,
                        &current_chromosome,
                        &self_clone.chromosome).expect("Error calling find_l_linked_loci() within adaptive_ld_knn_imputation() method for GenotypesAndPhenotypes trait.");
                // Using the linked loci, estimate the pairwise genetic distance between the current pool and the other pools
                let distances_all_loci = calculate_genetic_distances_between_pools(
                    i,
                    &linked_loci_idx,
                    &self_clone.intercept_and_allele_frequencies)
                    .expect("Error calling calculate_genetic_distances_between_pools() within adaptive_ld_knn_imputation() method for GenotypesAndPhenotypes trait.");
                // Find the k-nearest neighbours given the maximum distance and/or minimum k-neighbours (shadowing the distances across all pools with distances across k-nearest neighbours)
                let (distances, frequencies) =
                    find_k_nearest_neighbours(
                        &distances_all_loci,
                        max_pool_dist,
                        min_k_neighbours,
                        j,
                        &self_clone.intercept_and_allele_frequencies,
                    )
                    .expect("Error calling find_k_nearest_neighbours() within adaptive_ld_knn_imputation() method for GenotypesAndPhenotypes trait.");
                // Impute missing allele frequencies at the current locus
                *q = impute_allele_frequencies(&frequencies, &distances).expect("Error calling impute_allele_frequencies() within adaptive_ld_knn_imputation() method for GenotypesAndPhenotypes trait.");
            }
        });
        // Correct for allele frequency over- and under-flows if we have more than 1 allele representing each locus
        for j in 0..(loci_idx.len() - 1) {
            let idx_locus_ini = loci_idx[j];
            let idx_locus_fin = loci_idx[j + 1];
            for i in 0..n {
                if self_clone.intercept_and_allele_frequencies[(i, idx_locus_ini)].is_nan() {
                    let sum = self
                        .intercept_and_allele_frequencies
                        .slice(s![i, idx_locus_ini..idx_locus_fin])
                        .sum();
                    if sum != 1.0 {
                        if (idx_locus_fin - idx_locus_ini) == 2 {
                            self.intercept_and_allele_frequencies[(i, idx_locus_ini + 1)] =
                                1.00 - self.intercept_and_allele_frequencies[(i, idx_locus_ini)]
                        } else {
                            for k in idx_locus_ini..idx_locus_fin {
                                self.intercept_and_allele_frequencies[(i, k)] /= sum;
                            }
                        }
                    }
                    self.coverages[(i, j)] = f64::INFINITY;
                }
            }
        }
        Ok(self)
    }
}

pub fn impute_aldknni(
    mut genotypes_and_phenotypes: GenotypesAndPhenotypes,
    filter_stats: &FilterStats,
    min_loci_corr: &f64,
    max_pool_dist: &f64,
    min_l_loci: &u64,
    min_k_neighbours: &u64,
    restrict_linked_loci_per_chromosome: bool,
    optimise_n_steps_min_loci_corr: &usize,
    optimise_n_steps_max_pool_dist: &usize,
    optimise_max_l_loci: &u64,
    optimise_max_k_neighbours: &u64,
    optimise_n_reps: &usize,
    n_threads: &usize,
    out: &String,
) -> io::Result<String> {
    // Will need to remove window-related parameters and rename min_l_loci and min_l_loci to be main parameters
    // Also note that the we are no longer reverting to MVI as setting and optimising for min_l and min_k clashes with that idea

    // Extract loci indices
    let (loci_idx, _loci_chr, _loci_pos) = genotypes_and_phenotypes.count_loci().expect("Error calling count_loci() method within adaptive_ld_knn_imputation() method for GenotypesAndPhenotypes struct.");

    // Clone
    let self_clone = genotypes_and_phenotypes.clone();

    // Calculate LD across the entire genome
    println!("Estimating linkage between loci across the entire genome.");
    let corr = calculate_genomewide_ld(
        &genotypes_and_phenotypes,
        restrict_linked_loci_per_chromosome,
    )
    .expect("Error estimating pairwise linkage between loci across the entire genome.");
    println!("Optimising and/or estimating imputation accuracy.");
    let (
        optimum_min_loci_corr,
        optimum_max_pool_dist,
        optimum_min_l_loci,
        optimum_min_k_neighbours,
        optimum_mae,
    ) = optimise_params_and_estimate_accuracy(
        &genotypes_and_phenotypes,
        &self_clone,
        min_loci_corr,
        max_pool_dist,
        min_l_loci,
        min_k_neighbours,
        &loci_idx,
        &corr,
        restrict_linked_loci_per_chromosome,
        optimise_n_steps_min_loci_corr,
        optimise_n_steps_max_pool_dist,
        optimise_max_l_loci,
        optimise_max_k_neighbours,
        optimise_n_reps,
    )
    .expect("Error calling optimise_params_and_estimate_accuracy() in impute_aldknni().");
    println!("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    println!(
        "Minimum loci correlation threshold: {}",
        optimum_min_loci_corr
    );
    println!(
        "Maximum neighbour distance threshold: {}",
        optimum_max_pool_dist
    );
    println!("Minimum number of linked loci: {}", optimum_min_l_loci);
    println!(
        "Minimum number of k-nearest neighbours: {}",
        optimum_min_k_neighbours
    );
    println!(
        "Estimated imputation accuracy in terms of mean absolute error: {}",
        optimum_mae
    );
    println!("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .adaptive_ld_knn_imputation(
            &self_clone,
            &optimum_min_loci_corr,
            &optimum_max_pool_dist,
            &optimum_min_l_loci,
            &optimum_min_k_neighbours,
            &loci_idx,
            &corr,
            restrict_linked_loci_per_chromosome,
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
        let corr = calculate_genomewide_ld(&frequencies_and_phenotypes, false).unwrap();
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
        let (linked_loci_idx, correlations) = find_l_linked_loci(
            7020,
            &corr,
            &0.99,
            500,
            false,
            &"chr1".to_owned(),
            &frequencies_and_phenotypes.chromosome,
        )
        .unwrap();
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
        let _idx_locus_fin = 3;
        let (distances, frequencies) = find_k_nearest_neighbours(
            &distances,
            &0.5,
            3,
            idx_locus_ini,
            &frequencies_and_phenotypes.intercept_and_allele_frequencies,
        )
        .unwrap();
        // println!("idx_neighbours={:?}", idx_neighbours);
        // println!("distances={:?}", distances);
        // println!("frequencies={:?}", frequencies);
        assert_eq!(frequencies.len(), 3);
        // IMPUTATION
        let imputed_freq = impute_allele_frequencies(&frequencies, &distances).unwrap();
        println!("imputed_freq={:?}", imputed_freq);
        assert_eq!(imputed_freq, 0.2249532554929301);

        let _frac_top_missing_pools = 0.0;
        let _frac_top_missing_loci = 0.2;
        let min_loci_corr = 0.9;
        let max_pool_dist = 0.1;
        let min_l_loci = 10;
        let min_k_neighbours = 10;
        let restrict_linked_loci_per_chromosome = true;

        let optimise_n_steps_min_loci_corr = 1;
        let optimise_n_steps_max_pool_dist = 1;
        let optimise_max_l_loci = 100;
        let optimise_max_k_neighbours = 50;

        let optimise_n_reps = 1;

        let (loci_idx, _loci_chr, _loci_pos) = frequencies_and_phenotypes.count_loci().unwrap();
        let frequencies_and_phenotypes_clone = frequencies_and_phenotypes.clone();

        let _ = frequencies_and_phenotypes
            .adaptive_ld_knn_imputation(
                &frequencies_and_phenotypes_clone,
                &min_loci_corr,
                &max_pool_dist,
                &min_l_loci,
                &min_k_neighbours,
                &loci_idx,
                &corr,
                false,
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
        assert_eq!(n_nan, 1_915); // corresponds to the 1_915 alleles completely missing across all pools

        let keep_p_minus_1 = false;
        let genotypes_and_phenotypes = file_sync_phen
            .into_genotypes_and_phenotypes(&filter_stats, keep_p_minus_1, &n_threads)
            .unwrap();

        let outname = impute_aldknni(
            genotypes_and_phenotypes,
            &filter_stats,
            &min_loci_corr,
            &max_pool_dist,
            &min_l_loci,
            &min_k_neighbours,
            restrict_linked_loci_per_chromosome,
            &optimise_n_steps_min_loci_corr,
            &optimise_n_steps_max_pool_dist,
            &optimise_max_l_loci,
            &optimise_max_k_neighbours,
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
