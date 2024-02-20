use ndarray::{prelude::*, Zip};
use rand::Rng;
use std::cmp::Ordering;
use std::fs::{remove_file, OpenOptions};
use std::io;
use std::io::Write;
use std::time::{SystemTime, UNIX_EPOCH};

use crate::helpers::*;

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
                && (genotypes_and_phenotypes.chromosome[j] != current_chromosome)
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
    idx_row: &usize,
    linked_loci_idx: &[usize],
    intercept_and_allele_frequencies: &Array2<f64>,
) -> io::Result<Array1<f64>> {
    // Errors and f64::NAN are all converted into the maximum possible distance of 1.00 for simplicity
    let (n, _p) = intercept_and_allele_frequencies.dim();
    let q: Array1<f64> = intercept_and_allele_frequencies
        .row(*idx_row)
        .select(Axis(0), linked_loci_idx);
    let mut distances: Array1<f64> = Array1::from_elem(n, 1.0);
    Zip::indexed(&mut distances).par_for_each(|i, d| {
        if i != *idx_row {
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
    corr: &[Vec<u8>],
    min_loci_corr: &f64,
    min_l_loci: &usize,
    restrict_linked_loci_per_chromosome: bool,
    current_chromosome: &str,
    chromosomes: &[String],
) -> io::Result<(Vec<usize>, Vec<f64>)> {
    // Find loci passing the min_loci_corr threshold, then sort and add more loci if it did not pass the min_l_loci threshold
    //      - to optimise for the thresholds, set min_l_loci to 0
    //      - to optimise for counts, set min_loci_corr to 0.0
    assert!(
        *min_l_loci > 0,
        "Error: the minimum number of linked loci need to be greater than 0."
    );
    let p = corr.len();
    let mut vec_corr: Vec<f64> = vec![];
    let mut vec_idx: Vec<usize> = vec![];
    // Across rows and at the idx_col of the triangular matrix of correlations
    for i in 0..idx_col {
        if restrict_linked_loci_per_chromosome && (chromosomes[i] != *current_chromosome) {
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
            && (chromosomes[idx_col + j + 1] != *current_chromosome)
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
    let (vec_idx, vec_corr) = if (vec_idx.len() < *min_l_loci) || (*min_loci_corr == 0.0) {
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
        let l_linked_loci = if vec_corr.len() < *min_l_loci {
            vec_corr.len()
        } else {
            *min_l_loci
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
    min_k_neighbours: &usize,
    idx_col: &usize,
    intercept_and_allele_frequencies: &Array2<f64>,
) -> io::Result<(Vec<f64>, Array1<f64>)> {
    // Find neighbours passing the max_pool_dist threshold, then sort and add more neighbours if it did not pass the min_k_neighbours threshold
    //      - to optimise for the thresholds, set min_k_neighbours to 0
    //      - to optimise for counts, set max_pool_dist to 1.0 (maximum possible distance, i.e. mean absolute difference of 1.0 in allele frequency)
    assert!(
        *min_k_neighbours > 0,
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
            let bool_freq_nan = !intercept_and_allele_frequencies[(i, *idx_col)].is_nan();
            if bool_dist_nan && bool_dist_thresh && bool_freq_nan {
                idx.push(i);
                dist.push(distances[i]);
            }
        }
    }
    // If less than the minimum number of neighbours passed the threshold,
    // or if we are not filtering by maximum pool distance (i.e. max_pool_dist = 1.00),
    // then we sort the distances (increasing) and pick the nearest neighbours
    if (idx.len() < *min_k_neighbours) || (*max_pool_dist == 1.00) {
        idx = vec![];
        dist = vec![];
        for i in indices.clone().into_iter() {
            let bool_dist_nan = !distances[i].is_nan();
            let bool_freq_nan = !intercept_and_allele_frequencies[(i, *idx_col)].is_nan();
            if idx.len() == *min_k_neighbours {
                break;
            }
            if bool_dist_nan && bool_freq_nan {
                idx.push(i);
                dist.push(distances[i]);
            }
        }
    }
    // Extract non-missing allele frequencies at the locus requiring imputation of the k-nearest neighbours
    let mut freq: Array1<f64> = Array1::from_elem(idx.len(), f64::NAN);
    for i in 0..idx.len() {
        freq[i] = intercept_and_allele_frequencies[(idx[i], *idx_col)];
    }
    Ok((dist, freq))
}

fn impute_allele_frequencies(frequencies: &Array1<f64>, distances: &[f64]) -> io::Result<f64> {
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
    fn per_chunk_aldknni(
        &self,
        loci_arguments: (&usize, &usize, &[usize]),
        optimisation_arguments: (&[f64], &[f64], &usize, &usize, &usize),
        corr_matrix_arguments: (&[Vec<u8>], bool),
    ) -> io::Result<(Array2<f64>, f64, f64)> {
        // Parse arguments
        let idx_loci_idx_ini: &usize = loci_arguments.0;
        let idx_loci_idx_fin: &usize = loci_arguments.1;
        let loci_idx: &[usize] = loci_arguments.2;
        let vec_min_loci_corr: &[f64] = optimisation_arguments.0;
        let vec_max_pool_dist: &[f64] = optimisation_arguments.1;
        let min_l_loci: &usize = optimisation_arguments.2;
        let min_k_neighbours: &usize = optimisation_arguments.3;
        let n_reps: &usize = optimisation_arguments.4;
        let corr: &[Vec<u8>] = corr_matrix_arguments.0;
        let restrict_linked_loci_per_chromosome: bool = corr_matrix_arguments.1;
        // Define loci indices
        let idx_ini = loci_idx[*idx_loci_idx_ini];
        let idx_fin = loci_idx[*idx_loci_idx_fin];
        let mut allele_frequencies: Array2<f64> = self
            .intercept_and_allele_frequencies
            .slice(s![.., idx_ini..idx_fin])
            .to_owned();
        // Define the 2D array used for storing the minimum mean absolute error (MAE) estimates across all missing data across pools and loci.
        // We encode the MAE as u8 for memory efficiency, where we set 0: u8 as missing.
        let mut mae_u8: Array2<u8> = Array::from_elem(allele_frequencies.dim(), 0);
        Zip::indexed(&mut allele_frequencies)
        .and(&mut mae_u8)
        .par_for_each(|(i, j_local), q, mu8| {
            // Define global locus index
            let j = j_local + idx_ini;
            // Define the allele frequencies at the current allele/locus
            let vec_q: ArrayView1<f64> = self.intercept_and_allele_frequencies.column(j);
            let n_non_missing = vec_q.fold(0, |t, &x| if !x.is_nan() {t+1}else{t});
            let n_reps = if *n_reps <= n_non_missing {
                *n_reps
            } else {
                n_non_missing
            };
            if q.is_nan() && (n_reps > 0) {
                // Define the current chromosome
                let current_chromosome = self.chromosome[j].to_owned();
                // Select the n_reps most correlated non-missing pools which will be used for imputation/optimisation/estimation of imputation error using all the loci
                let distances_from_all_other_pools = calculate_genetic_distances_between_pools(
                    &i,
                    &(0..self.intercept_and_allele_frequencies.ncols()).collect::<Vec<usize>>(),
                    &self.intercept_and_allele_frequencies)
                    .expect("Error getting distances of all the pools using the 10 most linked loci above.");
                let mut idx_pools: Vec<usize> = (0..distances_from_all_other_pools.len()).collect();
                idx_pools
                    .sort_by(|&a, &b| distances_from_all_other_pools[a]
                        .partial_cmp(&distances_from_all_other_pools[b])
                        .expect("Error sorting indexes"));
                let idx_n_reps_nearest_pools: Vec<usize> = idx_pools[0..n_reps].to_vec();
                // Optimum mae, and parameters
                let mut optimum_mae: f64 = 1.0;
                let mut optimum_min_loci_corr: f64 = vec_min_loci_corr[0];
                let mut optimum_max_pool_dist: f64 = vec_max_pool_dist[0];
                // Grid search optimisation to find the optimal min_loci_corr and max_pool_dist which minimise imputation error (MAE: mean absolute error)
                // Across minimum loci correlation thresholds
                for min_loci_corr in vec_min_loci_corr.iter() {
                    // Find loci most correlated to the major allele of the current locus, i.e. the first allele of the locus as they were sorted by decreasing allele frequency (see Sort trait)
                    let (linked_loci_idx, _correlations) =
                        find_l_linked_loci(j, corr, min_loci_corr, min_l_loci,
                            restrict_linked_loci_per_chromosome,
                            &current_chromosome,
                            &self.chromosome).expect("Error calling find_l_linked_loci() within per_chunk_aldknni() method for GenotypesAndPhenotypes trait.");
                    // Across maximum pool distance thresholds
                    for max_pool_dist in vec_max_pool_dist.iter() {
                        // Across reps
                        let mut mae = 0.0;
                        for idx_i in idx_n_reps_nearest_pools.iter() {
                            // Using the linked loci, estimate the pairwise genetic distance between the current pool and the other pools
                            let distances_from_all_other_pools = calculate_genetic_distances_between_pools(
                                idx_i,
                                &linked_loci_idx,
                                &self.intercept_and_allele_frequencies)
                                .expect("Error calling calculate_genetic_distances_between_pools() within per_chunk_aldknni() method for GenotypesAndPhenotypes trait.");
                            // Find the k-nearest neighbours given the maximum distance and/or minimum k-neighbours (shadowing the distances across all pools with distances across k-nearest neighbours)
                            let (distances, frequencies) =
                                find_k_nearest_neighbours(
                                    &distances_from_all_other_pools,
                                    max_pool_dist,
                                    min_k_neighbours,
                                    &j,
                                    &self.intercept_and_allele_frequencies,
                                )
                                .expect("Error calling find_k_nearest_neighbours() within per_chunk_aldknni() method for GenotypesAndPhenotypes trait.");
                            // Impute and find the error
                            mae += (vec_q[*idx_i] - impute_allele_frequencies(&frequencies, &distances).expect("Error calling impute_allele_frequencies() within per_chunk_aldknni() method for GenotypesAndPhenotypes trait.")
                            ).abs();
                        }
                        mae /= n_reps as f64;
                        if mae < optimum_mae {
                            optimum_mae = mae;
                            optimum_min_loci_corr = *min_loci_corr;
                            optimum_max_pool_dist = *max_pool_dist;
                        }
                    }
                }
                // Take note of the minimum MAE which we code are u8 for memory efficiency
                // We use 254 instead of 255 and add 1 because we set 0 as missing above
                *mu8 = (optimum_mae * 254.0).round() as u8 + 1;
                // Impute actual missing data point (ith pool and jth locus)
                // Find loci most correlated to the major allele of the current locus, i.e. the first allele of the locus as they were sorted by decreasing allele frequency (see Sort trait)
                let (linked_loci_idx, _correlations) =
                    find_l_linked_loci(j, corr, &optimum_min_loci_corr, min_l_loci,
                        restrict_linked_loci_per_chromosome,
                        &current_chromosome,
                        &self.chromosome).expect("Error calling find_l_linked_loci() within per_chunk_aldknni() method for GenotypesAndPhenotypes trait.");
                // Using the linked loci, estimate the pairwise genetic distance between the current pool and the other pools
                let distances_from_all_other_pools = calculate_genetic_distances_between_pools(
                    &i,
                    &linked_loci_idx,
                    &self.intercept_and_allele_frequencies)
                    .expect("Error calling calculate_genetic_distances_between_pools() within per_chunk_aldknni() method for GenotypesAndPhenotypes trait.");
                // Find the k-nearest neighbours given the maximum distance and/or minimum k-neighbours (shadowing the distances across all pools with distances across k-nearest neighbours)
                let (distances, frequencies) =
                    find_k_nearest_neighbours(
                        &distances_from_all_other_pools,
                        &optimum_max_pool_dist,
                        min_k_neighbours,
                        &j,
                        &self.intercept_and_allele_frequencies,
                    )
                    .expect("Error calling find_k_nearest_neighbours() within per_chunk_aldknni() method for GenotypesAndPhenotypes trait.");
                // Impute missing allele frequencies at the current locus
                *q = impute_allele_frequencies(&frequencies, &distances).expect("Error calling impute_allele_frequencies() within per_chunk_aldknni() method for GenotypesAndPhenotypes trait.");
            }
        });
        // Extract average minimum MAE across loci and pools (Note that we used 0: u8 as missing)
        let mut sum_mae = 0.0;
        let mut n_missing = 0.0;
        for mu8 in mae_u8.into_iter() {
            if mu8 > 0 {
                sum_mae += (mu8 as f64 - 1.0) / 254.0;
                n_missing += 1.0;
            }
        }
        // Correct for allele frequency over- and under-flows, as we are assuming all loci are represented by all of its alleles (see assumptions above)
        // Local positions in the chunk (add loci_idx[*idx_loci_idx_ini] to get the global index)
        let mut vec_chunk_loci_idx: Vec<usize> = vec![];
        for idx in loci_idx
            .iter()
            .take(*idx_loci_idx_fin + 1)
            .skip(*idx_loci_idx_ini)
        {
            vec_chunk_loci_idx.push(idx - idx_ini);
        }
        let n = allele_frequencies.nrows();
        for j in 0..(vec_chunk_loci_idx.len() - 1) {
            let idx_locus_ini = vec_chunk_loci_idx[j];
            let idx_locus_fin = vec_chunk_loci_idx[j + 1];
            for i in 0..n {
                if self.intercept_and_allele_frequencies
                    [(i, (idx_locus_ini + loci_idx[*idx_loci_idx_ini]))]
                    .is_nan()
                {
                    let sum = allele_frequencies
                        .slice(s![i, idx_locus_ini..idx_locus_fin])
                        .sum();
                    if sum != 1.0 {
                        for k in idx_locus_ini..idx_locus_fin {
                            allele_frequencies[(i, k)] /= sum;
                        }
                    }
                }
            }
        }
        Ok((allele_frequencies, sum_mae, n_missing))
    }

    pub fn adaptive_ld_knn_imputation(
        &self,
        loci_idx: &[usize],
        optimisation_arguments: (&f64, &f64, &usize, &usize, &usize),
        corr_matrix_arguments: (&[Vec<u8>], bool),
        n_threads: &usize,
        prefix: &str,
    ) -> io::Result<(String, f64)> {
        // We are assuming that all non-zero alleles across pools are kept, i.e. biallelic loci have 2 columns, triallelic have 3, and so on.
        //      - Input vcf file will have all alleles per locus extracted.
        //      - Similarly, input sync file will have all alleles per locus extracted.
        //      - Finally, input allele frequency table file which can be represented by all alleles or one less allele per locus will be appended with the alternative allele if the sum of allele frequencies per locus do not add up to one.
        self.check().expect("Error self.check() within adaptive_ld_knn_imputation() method of GenotypesAndPhenotypes struct.");
        // Parse arguments
        let min_loci_corr: &f64 = optimisation_arguments.0;
        let max_pool_dist: &f64 = optimisation_arguments.1;
        let min_l_loci: &usize = optimisation_arguments.2;
        let min_k_neighbours: &usize = optimisation_arguments.3;
        let n_reps: &usize = optimisation_arguments.4;
        let corr: &[Vec<u8>] = corr_matrix_arguments.0;
        let restrict_linked_loci_per_chromosome: bool = corr_matrix_arguments.1;
        // Define fixed linkage and distance thresholds, i.e. the minimum number of loci to use in estimating genetic distances, and the minimum number of nearest neighbours to include in imputation
        let (n, p) = self.intercept_and_allele_frequencies.dim();
        let min_l_loci = if *min_l_loci >= p { p - 1 } else { *min_l_loci };
        let min_k_neighbours = if *min_k_neighbours >= n {
            n - 1
        } else {
            *min_k_neighbours
        };
        // Define the range of minimum loci correlation and maximum pool distance thresholds to be used for optimisation.
        // However, if the user-supplied values are non-missing then that user-supplied fixed threshold will be used without optimisation
        // We have set a hard-coded number of steps which can be explored across the 2 parameter spaces, i.e. 10 steps provides a reasonable balance between accuracy and computational efficiency
        let vec_min_loci_corr: Vec<f64> = if min_loci_corr.is_nan() {
            vec![0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0, 1.0]
        } else {
            vec![*min_loci_corr]
        };
        let vec_max_pool_dist: Vec<f64> = if max_pool_dist.is_nan() {
            vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.0]
        } else {
            vec![*max_pool_dist]
        };
        // Define chunks which respect loci groupings
        let mut n_chunks = *n_threads; // Should be equal to the number of threads because std::thread::scope will wait for other threads to finish before starting with another thread once it finishes
        let l = loci_idx.len();
        let chunk_size = l / n_chunks;
        // Define the indices of the indices of loci
        let vec_idx_loci_idx_ini: Vec<usize> = (0..(l - chunk_size)).step_by(chunk_size).collect();
        let mut vec_idx_loci_idx_fin: Vec<usize> = (chunk_size..l).step_by(chunk_size).collect();
        let idx_fin = vec_idx_loci_idx_fin[vec_idx_loci_idx_fin.len() - 1];
        if idx_fin != (l - 1) {
            vec_idx_loci_idx_fin.pop();
            vec_idx_loci_idx_fin.push(l - 1);
        }
        assert_eq!(vec_idx_loci_idx_ini.len(), vec_idx_loci_idx_fin.len());
        n_chunks = vec_idx_loci_idx_ini.len();
        // Instantiate thread object for parallel execution
        let mut thread_objects: Vec<(String, f64, f64)> = vec![];
        for i in 0..n_chunks {
            let idx_loci_idx_ini = vec_idx_loci_idx_ini[i];
            let idx_loci_idx_fin = vec_idx_loci_idx_fin[i];
            let thread = std::thread::scope(|_scope| {
                let (allele_frequencies, sum_mae, n_missing) = self
                    .per_chunk_aldknni(
                        (&idx_loci_idx_ini, &idx_loci_idx_fin, loci_idx),
                        (
                            &vec_min_loci_corr,
                            &vec_max_pool_dist,
                            &min_l_loci,
                            &min_k_neighbours,
                            n_reps,
                        ),
                        (corr, restrict_linked_loci_per_chromosome),
                    )
                    .expect("Error executing per_chunk_aldknni() method on self_clone.");
                // Write-out intermediate file
                let idx_ini = loci_idx[idx_loci_idx_ini];
                let idx_fin = loci_idx[idx_loci_idx_fin];
                let time = SystemTime::now()
                        .duration_since(UNIX_EPOCH)
                        .expect("Error extracting time in UNIX_EPOCH within write_csv() method for GenotypesAndPhenotypes struct.")
                        .as_secs_f64();
                let mut rng = rand::thread_rng();
                let random_number = rng.gen_range(1_000_000..10_000_000);
                let n_digits: usize = loci_idx
                    .last()
                    .expect("Error extracting the last element of loci_idx. Probably empty.")
                    .to_owned()
                    .to_string()
                    .len();
                let mut start_index = idx_ini.to_string();
                let mut end_index = idx_fin.to_string();
                for _ in 0..(n_digits - start_index.len()) {
                    start_index = "0".to_owned() + &start_index;
                }
                for _ in 0..(n_digits - end_index.len()) {
                    end_index = "0".to_owned() + &end_index;
                }
                let fname_intermediate_file: String = prefix.to_owned()
                    + "-"
                    + &start_index
                    + "-"
                    + &end_index
                    + "-"
                    + &time.to_string()
                    + &random_number.to_string()
                    + ".tmp";
                let chromosome = self.chromosome[idx_ini..idx_fin].to_owned();
                let position = self.position[idx_ini..idx_fin].to_owned();
                let allele = self.allele[idx_ini..idx_fin].to_owned();
                let p = allele_frequencies.ncols();
                assert_eq!(chromosome.len(), p, "Error, the number of chromosome names and the total number of loci are not equal.");
                // Instantiate output file
                println!(
                    "--> {}: Writing out intermediate file with expected MAE of {}: {}",
                    i,
                    sensible_round(sum_mae / n_missing, 4),
                    &fname_intermediate_file
                );
                let error_writing_file =
                    "Unable to create file: ".to_owned() + &fname_intermediate_file;
                let mut file_out = OpenOptions::new()
                    .create_new(true)
                    .write(true)
                    .append(false)
                    .open(&fname_intermediate_file)
                    .expect(&error_writing_file);
                // Write the header only for the first chunk
                if i == 0 {
                    file_out
                        .write_all(
                            ("#chr,pos,allele,".to_owned() + &self.pool_names.join(",") + "\n").as_bytes(),
                        )
                        .expect("Error calling write_all() within the write_csv() method for GenotypesAndPhenotypes struct.");
                }
                // Write allele frequencies line by line
                for j in 0..p {
                    let line = [
                        chromosome[j].to_owned(),
                        position[j].to_string(),
                        allele[j].to_owned(),
                        allele_frequencies
                            .column(j)
                            .iter()
                            .map(|&x| parse_f64_roundup_and_own(x, 6))
                            .collect::<Vec<String>>()
                            .join(","),
                    ]
                    .join(",")
                        + "\n";
                    file_out.write_all(line.as_bytes()).expect("Error calling write_all() per line of the output file within the write_csv() method for GenotypesAndPhenotypes struct.");
                }
                (fname_intermediate_file, sum_mae, n_missing)
            });
            thread_objects.push(thread);
        }
        // Concatenate output per chunk
        let mut vec_fname_intermediate_files_and_mae: Vec<(String, f64, f64)> = vec![];
        for thread in thread_objects {
            vec_fname_intermediate_files_and_mae.push(thread);
        }
        // Sort by intermediate output filenames (named according to indices)
        vec_fname_intermediate_files_and_mae.sort_by(|a, b| {
            a.0.partial_cmp(&b.0)
                .expect("Error sorting intermediate output filenames")
        });
        // Concatenate intermediate output filenames together and calculate mae
        let mut file_0 = OpenOptions::new()
            .append(true)
            .open(&vec_fname_intermediate_files_and_mae[0].0)
            .expect("Error opening the intermediate file of the first chunk.");
        let mut sum_mae = vec_fname_intermediate_files_and_mae[0].1;
        let mut n_missing = vec_fname_intermediate_files_and_mae[0].2;
        for name_and_mae in vec_fname_intermediate_files_and_mae.iter().skip(1) {
            let mut file_1 = OpenOptions::new()
                .read(true)
                .open(&name_and_mae.0)
                .expect("Error opening the intermediate file of a chunk.");
            io::copy(&mut file_1, &mut file_0)
                .expect("Error concatenating intermediate output files.");
            remove_file(&name_and_mae.0).expect("Error removing the intermediate file of a chunk.");
            sum_mae += name_and_mae.1;
            n_missing += name_and_mae.2;
        }
        let mae = sum_mae / n_missing;
        Ok((vec_fname_intermediate_files_and_mae[0].0.to_owned(), mae))
    }
}

/// # `aldknni`: adaptive linkage disequilibrium (LD)-based k-nearest neighbour imputation of genotype data
///
/// This is an attempt to extend the [LD-kNNi method of Money et al, 2015, i.e. LinkImpute](https://doi.org/10.1534/g3.115.021667), which was an extension of the [kNN imputation of Troyanskaya et al, 2001](https://doi.org/10.1093/bioinformatics/17.6.520). Similar to LD-kNNi, linkage disequilibrium (LD) is estimated using Pearson's product moment correlation per pair of loci, which is computed per chromosome by default, but can be computed across the entire genome. We use the mean absolute difference/error (MAE) between allele frequencies among linked loci as an estimate of genetic distance between samples. Fixed values for the minimum correlation to identify loci used in distance estimation, and maximum genetic distance to select the k-nearest neighbours can be defined. Additionally, minimum number of loci to include in distance estimation, and minimum number of nearest neighbours can be set. Moreover, all four parameters can be optimised, i.e. the minimum correlation and/or maximum distance and/or minimum number of loci and/or minimum number of nearest neighbours which minimises the MAE between predicted and expected allele frequencies after simulating 10% missing data are identified.
///
/// The allele depth information (`AD`), i.e. the unfiltered allele depth which includes the reads which did not pass the variant caller filters are used to calculate allele frequencies. If the `GT` field is present but the `AD` field is absent, then each sample is assumed to be an individual diploid, i.e., neither a polyploid nor a pool. Optional filtering steps based on minimum depth, minimum allele frequency, and maximum sparsity are available. Genotype data are not imported into R, LD estimation and imputation per se are multi-threaded, and imputation output is written into disk as an [allele frequency table](#allele-frequency-table-csv). The structs, traits, methods, and functions defined in this library are subsets of [poolgen](https://github.com/jeffersonfparil/poolgen), and will eventually be merged.
///
/// The imputed allele frequency is computed as:
///
/// $$
/// \hat q_{r,j} = { \sum_{i \ne r}^{k} q_{i,j} (1 - \delta_{i,r}) }
/// $$
///
/// with:
///
/// $$
/// \delta_{i,r} = { {1 \over \sum d_{i,r}} d_{i,r} }
/// $$
///
/// and
///
/// $$
/// d_{i,r} = { {1 \over c} { \sum_{j=1}^{c} |q_{i,j} - q_{r,j}| } }
/// $$
///
/// where:
///
/// - $\hat q_{r,j}$ is the imputed allele frequency of sample $r$ at the $j^{\text {th}}$ locus,
/// - $n$ is the total number of samples,
/// - $m$ is the number of samples which are missing data at the $j^{\text {th}}$ locus,
/// - $q_{i,j}$ is the known allele frequency of the $i^{\text {th}}$ sample at the $j^{\text {th}}$ locus,
/// - $k$ is the number of nearest neighbours or the samples most closely related to the sample requiring imputation, i.e. sample $r$ at locus $j$, and
/// - $\delta_{i,r}$ is scaled $d_{i,r}$ which is the genetic distance between the $i^{\text {th}}$ sample and sample $r$. This distance is the mean absolute difference in allele frequencies between the two samples across $c$ linked loci.
///
/// The variables $k$ and $c$ are proportional to the user inputs `max_pool_dist` (default=0.1) and `min_loci_corr` (default=0.9), respectively. The former defines the maximum distance of samples to be considered as one of the k-nearest neighbours, while the latter refers to the minimum correlation with the locus requiring imputation to be included in the estimation of the genetic distance.
pub fn impute_aldknni(
    genotypes_and_phenotypes: GenotypesAndPhenotypes,
    filter_stats: &FilterStats,
    optimisation_arguments: (&f64, &f64, &usize, &usize, &usize),
    restrict_linked_loci_per_chromosome: bool,
    n_threads: &usize,
    out: &str,
) -> io::Result<String> {
    // Parse arguments
    let min_loci_corr: &f64 = optimisation_arguments.0;
    let max_pool_dist: &f64 = optimisation_arguments.1;
    // Extract loci indices
    let (loci_idx, _loci_chr, _loci_pos) = genotypes_and_phenotypes.count_loci().expect("Error calling count_loci() method within adaptive_ld_knn_imputation() method for GenotypesAndPhenotypes struct.");
    // Calculate LD across the entire genome
    println!("Estimating linkage between loci across the entire genome.");
    let corr = calculate_genomewide_ld(
        &genotypes_and_phenotypes,
        restrict_linked_loci_per_chromosome,
    )
    .expect("Error estimating pairwise linkage between loci across the entire genome.");
    if min_loci_corr.is_nan() || max_pool_dist.is_nan() {
        println!("Optimising and estimating imputation accuracy.");
    } else {
        println!("Estimating imputation accuracy.");
    }
    let start = std::time::SystemTime::now();
    let (fname_imputed, mae) = genotypes_and_phenotypes
        .adaptive_ld_knn_imputation(
            &loci_idx,
            optimisation_arguments,
            (&corr, restrict_linked_loci_per_chromosome),
            n_threads,
            &(out.replace(".csv", "")),
        )
        .expect("Error calling adaptive_ld_knn_imputation() within impute_aldknni().");
    let end = std::time::SystemTime::now();
    let duration = end.duration_since(start).expect("Error measuring the duration of running adaptive_ld_knn_imputation() within impute_aldknni().");
    println!("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    println!(
        "Expected imputation accuracy in terms of mean absolute error: {}",
        mae
    );
    println!("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    let mut genotypes_and_phenotypes = FileGeno {
        filename: fname_imputed.to_owned(),
    }
    .convert_into_genotypes_and_phenotypes(filter_stats, false, n_threads)
    .expect("Error loading the imputed genotype file.");
    remove_file(fname_imputed).expect("Error removing concatenated intermediate file.");
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
            .convert_into_genotypes_and_phenotypes(&filter_stats, false, &n_threads)
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
            &500,
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
            &idx_pool,
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
            &3,
            &idx_locus_ini,
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

        let (loci_idx, _loci_chr, _loci_pos) = frequencies_and_phenotypes.count_loci().unwrap();
        // let frequencies_and_phenotypes_clone = frequencies_and_phenotypes.clone();
        let n_reps = 5;
        let min_loci_corr = f64::NAN;
        let max_pool_dist = f64::NAN;
        let min_l_loci = 10;
        let min_k_neighbours = 10;
        let restrict_linked_loci_per_chromosome = false;
        let n_threads = 8;

        let (_fname_imputed, mae) = frequencies_and_phenotypes
            .adaptive_ld_knn_imputation(
                &loci_idx,
                (
                    &min_loci_corr,
                    &max_pool_dist,
                    &min_l_loci,
                    &min_k_neighbours,
                    &n_reps,
                ),
                (&corr, restrict_linked_loci_per_chromosome),
                &n_threads,
                &"intermediate_output".to_owned(),
            )
            .unwrap();
        println!(
            "After imputation:\n{:?}",
            frequencies_and_phenotypes.intercept_and_allele_frequencies
        );
        println!("Estimated MAE={}", mae);
        // let n_nan = frequencies_and_phenotypes
        //     .intercept_and_allele_frequencies
        //     .iter()
        //     .fold(0, |n_nan, &x| if x.is_nan() { n_nan + 1 } else { n_nan });
        // println!("n_nan={}", n_nan);
        // assert_eq!(n_nan, 1_915); // corresponds to the 1_915 alleles completely missing across all pools

        let keep_p_minus_1 = false;
        let genotypes_and_phenotypes = file_sync_phen
            .convert_into_genotypes_and_phenotypes(&filter_stats, keep_p_minus_1, &n_threads)
            .unwrap();

        let outname = impute_aldknni(
            genotypes_and_phenotypes,
            &filter_stats,
            (
                &min_loci_corr,
                &max_pool_dist,
                &min_l_loci,
                &min_k_neighbours,
                &n_reps,
            ),
            restrict_linked_loci_per_chromosome,
            &n_threads,
            &"test-impute_aldknni.csv".to_owned(),
        )
        .unwrap();
        assert_eq!(outname, "test-impute_aldknni.csv".to_owned()); // Do better!!! Load data - thus working on improving load_table()

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
