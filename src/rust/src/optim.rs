use rand::prelude::IteratorRandom;
use std::io;
use ndarray::prelude::*;

use crate::structs_and_traits::*;

impl GenotypesAndPhenotypes {
    pub fn simulate_missing(
        &mut self,
        max_sparsity: &f64,
    ) -> io::Result<(
        &mut Self,
        f64,
        Vec<usize>,
        Vec<usize>,
        Vec<f64>,
        Vec<Vec<f64>>,
    )> {
        self.check().expect("Error calling check() method within simulate_missing() method for GenotypesAndPhenotypes struct.");
        let (n, l) = self.coverages.dim();
        let (loci_idx, _loci_chr, _loci_pos) = self.count_loci().expect("Error defining loci indexes and identities via count_loci() method within simulate_missing() method for GenotypesAndPhenotypes struct.");
        // Limit max_sparsity to 90%
        let max_sparsity = if *max_sparsity > 0.9 {
            0.9
        } else {
            *max_sparsity
        };
        // Define the number of loci to mask
        let n_masked = (((n * l) as f64) * max_sparsity).floor() as usize;
        // Sample random loci
        let mut rng = rand::thread_rng();
        let mut vec_masked_loci_idx_tmp = (0..(n * l)).choose_multiple(&mut rng, n_masked);
        vec_masked_loci_idx_tmp.sort();
        // Mask and extract the masked information
        let mut vec_masked_loci_idx: Vec<usize> = vec![];
        let mut vec_masked_coverages: Vec<f64> = vec![];
        let mut vec_vec_masked_alleles_freqs: Vec<Vec<f64>> = vec![];
        let mut idx_counter = 0;
        let mut idx = 0;
        'outer: for i in 0..n {
            for j in 0..l {
                if vec_masked_loci_idx_tmp[idx_counter] == idx {
                    if !self.coverages[(i, j)].is_nan() {
                        // Extract index of the locus
                        vec_masked_loci_idx.push(idx);
                        // Extract mask the coverages
                        vec_masked_coverages.push(self.coverages[(i, j)]);
                        self.coverages[(i, j)] = f64::NAN;
                        // Use the indexes of the locus extract masked frequencies and to set missing values to all alleles in the locus
                        let idx_ini = loci_idx[j];
                        let idx_fin = loci_idx[j + 1];
                        let mut vec_freqs: Vec<f64> = vec![];
                        for k in idx_ini..idx_fin {
                            vec_freqs.push(self.intercept_and_allele_frequencies[(i, k)]);
                            self.intercept_and_allele_frequencies[(i, k)] = f64::NAN;
                        }
                        vec_vec_masked_alleles_freqs.push(vec_freqs);
                    }
                    idx_counter += 1;
                    if idx_counter == vec_masked_loci_idx_tmp.len() {
                        break 'outer;
                    }
                }
                idx += 1;
            }
        }
        self.check().expect("Error calling check() method within simulate_missing() method for GenotypesAndPhenotypes struct.");
        Ok((
            self,
            max_sparsity,
            loci_idx,
            vec_masked_loci_idx,
            vec_masked_coverages,
            vec_vec_masked_alleles_freqs,
        ))
    }

    pub fn extract_imputed_mask(
        &self,
        loci_idx: &Vec<usize>,
        vec_masked_loci_idx: &Vec<usize>,
    ) -> io::Result<Vec<Vec<f64>>> {
        let (n, _p) = self.intercept_and_allele_frequencies.dim();
        let (_n, l) = self.coverages.dim();
        let mut vec_vec_imputed_mask: Vec<Vec<f64>> = vec![];
        let mut idx_counter = 0;
        let mut idx = 0;
        'outer: for i in 0..n {
            for j in 0..l {
                if vec_masked_loci_idx[idx_counter] == idx {
                    let idx_ini = loci_idx[j];
                    let idx_fin = loci_idx[j + 1];
                    let mut vec_freqs: Vec<f64> = vec![];
                    for k in idx_ini..idx_fin {
                        vec_freqs.push(self.intercept_and_allele_frequencies[(i, k)]);
                    }
                    vec_vec_imputed_mask.push(vec_freqs);
                    idx_counter += 1;
                    if idx_counter == vec_masked_loci_idx.len() {
                        break 'outer;
                    }
                }
                idx += 1;
            }
        }
        Ok(vec_vec_imputed_mask)
    }

    pub fn estimate_expected_mae_in_aldknni(
        &mut self,
        min_loci_corr: &f64,
        max_pool_dist: &f64,
        min_l_loci: &u64,
        min_k_neighbours: &u64,
        restrict_linked_loci_per_chromosome: bool,
        loci_idx: &Vec<usize>,
        vec_masked_loci_idx: &Vec<usize>,
        vec_vec_masked_alleles_freqs: &Vec<Vec<f64>>,
    ) -> io::Result<f64> {
        self.adaptive_ld_knn_imputation(
            min_loci_corr,
            max_pool_dist,
            min_l_loci,
            min_k_neighbours,
            restrict_linked_loci_per_chromosome,
        )
        .expect("Error calling adaptive_ld_knn_imputation() within estimate_expected_mae_in_aldknni() method for GenotypesAndPhenotypes struct.");
        // Extract imputed frequencies
        let vec_vec_imputed_mask = self
            .extract_imputed_mask(&loci_idx, &vec_masked_loci_idx)
            .expect("Error calling extract_imputed_mask() within estimate_expected_mae_in_aldknni() method for GenotypesAndPhenotypes struct.");
        let m = vec_vec_masked_alleles_freqs.len(); // total number of loci
        let mean_absolute_error = if m == 0 {
            // If no loci were simulated to be missing
            1.00
        } else {
            let mut sum_abs = 0.0;
            let mut a = 0.0; // total number of alleles across all loci
            for i in 0..m {
                let w = vec_vec_masked_alleles_freqs[i].len();
                for j in 0..w {
                    if (vec_vec_masked_alleles_freqs[i][j].is_nan())
                        | (vec_vec_imputed_mask[i][j].is_nan())
                    {
                        // Skip if we if masked already missing loci or if we cannot impute because the whole locus is missing
                        continue;
                    } else {
                        sum_abs +=
                            (vec_vec_masked_alleles_freqs[i][j] - vec_vec_imputed_mask[i][j]).abs();
                        // (vec_vec_masked_alleles_freqs[i][j] - vec_vec_imputed_mask[i][j]).powf(2.0); // MSE seem to have the same behaviour as MAE
                        a += 1.00;
                    }
                }
            }
            if a == 0.0 {
                1.00
            } else {
                sum_abs / a
            }
        };
        Ok(mean_absolute_error)
    }

    pub fn estimate_expected_mae_in_mvi(&self) -> io::Result<f64> {
        let mut genotype_data_for_optimisation = self.clone();
        let (n, p) = genotype_data_for_optimisation.intercept_and_allele_frequencies.dim();
        let missing_rate_sim = if (n * p) < 20_000 {
            0.01
        } else {
            10_000.0 / ((n * p) as f64)   
        };
        let (
            _,
            max_sparsity,
            loci_idx,
            vec_masked_loci_idx,
            _vec_masked_coverages,
            vec_vec_masked_alleles_freqs,
        ) = genotype_data_for_optimisation
            .simulate_missing(&missing_rate_sim)
            .expect("Error calling simulate_missing() within estimate_expected_mae_in_mvi() method for GenotypesAndPhenotypes struct.");
        assert_eq!(
            max_sparsity, missing_rate_sim,
            "Ooops! Missing unexpected simulated sparsity!"
        );
        let _ = genotype_data_for_optimisation.mean_imputation().expect("Error calling mean_imputation() within estimate_expected_mae_in_mvi() method for GenotypesAndPhenotypes struct.");
        let vec_vec_imputed_mask = genotype_data_for_optimisation
            .extract_imputed_mask(&loci_idx, &vec_masked_loci_idx)
            .expect("Error calling extract_imputed_mask() within estimate_expected_mae_in_mvi() method for GenotypesAndPhenotypes struct.");
        // Compare imputed and expected frequencies
        let m = vec_vec_masked_alleles_freqs.len(); // total number of loci
        let mean_absolute_error = if m == 0 {
            // If no loci were simulated to be missing
            1.00
        } else {
            let mut sum_abs = 0.0;
            let mut a = 0.0; // total number of alleles across all loci
            for i in 0..m {
                let w = vec_vec_masked_alleles_freqs[i].len();
                for j in 0..w {
                    if (vec_vec_masked_alleles_freqs[i][j].is_nan())
                        | (vec_vec_imputed_mask[i][j].is_nan())
                    {
                        // Skip if we if masked already missing loci or if we cannot impute because the whole locus is missing
                        continue;
                    } else {
                        sum_abs +=
                            (vec_vec_masked_alleles_freqs[i][j] - vec_vec_imputed_mask[i][j]).abs();
                        // (vec_vec_masked_alleles_freqs[i][j] - vec_vec_imputed_mask[i][j]).powf(2.0); // MSE seem to have the same behaviour as MAE
                        a += 1.00;
                    }
                }
            }
            if a == 0.0 {
                1.00
            } else {
                sum_abs / a
            }
        };
        Ok(mean_absolute_error)
    }
}

pub fn optimise_params_and_estimate_accuracy(
    genotypes_and_phenotypes: &GenotypesAndPhenotypes,
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
) -> io::Result<(f64, f64, u64, u64, f64)> {
    // Defining the ranges of min_loci_corr, max_pool_dist , min_l_loci and min_k_neighbours to test
    // Note: The number of steps may vary from the input depending on how we can evenly divide the range
    let vec_min_loci_corr: Vec<f64> = if *optimise_n_steps_min_loci_corr > 1 {
        let step_size = (101.0 / *optimise_n_steps_min_loci_corr as f64).round() as usize;
        (0..=100)
            .step_by(step_size)
            .map(|x| x as f64 / 100.0)
            .collect::<Vec<f64>>()
            .into_iter()
            .rev()
            .collect()
    } else {
        vec![*min_loci_corr]
    };
    let vec_max_pool_dist: Vec<f64> = if *optimise_n_steps_max_pool_dist > 1 {
        let step_size = (101.0 / *optimise_n_steps_max_pool_dist as f64).round() as usize;
        (0..=100)
            .step_by(step_size)
            .map(|x| x as f64 / 100.0)
            .collect()
    } else {
        vec![*max_pool_dist]
    };
    let vec_min_l_loci: Vec<u64> = if *optimise_max_l_loci > *min_l_loci {
        let step_size = 1;
        (1..=*optimise_max_l_loci)
            .step_by(step_size)
            .map(|x| x as u64)
            .collect()
    } else {
        vec![*min_l_loci]
    };
    let vec_min_k_neighbours: Vec<u64> = if *optimise_max_k_neighbours > *min_k_neighbours {
        let k_max = if *optimise_max_k_neighbours as usize
            > genotypes_and_phenotypes
                .intercept_and_allele_frequencies
                .nrows()
        {
            genotypes_and_phenotypes
                .intercept_and_allele_frequencies
                .nrows()
        } else {
            *optimise_max_k_neighbours as usize
        };
        let step_size = 1;
        (1..=k_max).step_by(step_size).map(|x| x as u64)
            .collect()
    } else {
        vec![*min_k_neighbours]
    };
    // Optimisation via coordinate-descent-like algorithm starting with the strictest parameters
    // Parameter optimisation order:
    //       (1) k-neighbours (from low to high)
    //       (2) maximum pool distances (from low to high)
    //       (3) l loci (from low to high)
    //       (4) minimum loci correlation (from high to low)
    // As we explore the parameter space one at a time at decreasing stringency,
    // we replace the optimum parameter if mae decreases, and
    // break if mae increases so that we are more computationally efficient.
    let all_parameters_are_fixed = (vec_min_loci_corr.len() == 1)
        & (vec_max_pool_dist.len() == 1)
        & (vec_min_l_loci.len() == 1)
        & (vec_min_k_neighbours.len() == 1);
    let optimise_n_reps = if all_parameters_are_fixed & (*optimise_n_reps > 1) {
        1
    } else {
        *optimise_n_reps
    };
    let mut min_loci_corr = vec![];
    let mut max_pool_dist = vec![];
    let mut min_l_loci = vec![];
    let mut min_k_neighbours = vec![];
    let mut mae = vec![];

    // let optimise_n_reps = 3;

    if all_parameters_are_fixed == false {
        println!("-----------------------------------------------");
        if optimise_n_reps > 1 {
            println!("rep\tmin_k_neighbours\tmax_pool_dist\tmin_l_loci\tmin_loci_corr\tmae");
        } else {
            println!("min_k_neighbours\tmax_pool_dist\tmin_l_loci\tmin_loci_corr\tmae");
        }
    }
    // Assumes smooth unimodal parameter spaces
    for r in 0..optimise_n_reps {
        // Simulate 10% sparsity to determine accuracy (in terms of MAE) and to optimise for the best k and l, if applicable, i.e. n_loci_to_estimate_distance==0 or k_neighbours==0
        let mut genotype_data_for_optimisation = genotypes_and_phenotypes.clone();
        let (n, p) = genotype_data_for_optimisation.intercept_and_allele_frequencies.dim();
        let missing_rate_sim = if (n * p) < 20_000 {
            0.01
        } else {
            10_000.0 / ((n * p) as f64)   
        };
        let (
            _,
            max_sparsity,
            loci_idx,
            vec_masked_loci_idx,
            _vec_masked_coverages,
            vec_vec_masked_alleles_freqs,
        ) = genotype_data_for_optimisation
            .simulate_missing(&missing_rate_sim)
            .expect("Error calling simulate_missing() method within optimise_params_and_estimate_accuracy().");
        assert_eq!(
            max_sparsity, missing_rate_sim,
            "Ooops! Missing unexpected simulated sparsity!"
        );
        // Coordinate-descent-like optimisation starting with the strictest parameters
        let mut optimum_min_loci_corr = vec_min_loci_corr
            .last()
            .expect("Error: vec_min_loci_corr is empty.")
            .to_owned();
        let mut optimum_max_pool_dist = vec_max_pool_dist[0];
        let mut optimum_min_l_loci = vec_min_l_loci[0];
        let mut optimum_min_k_neighbours = vec_min_k_neighbours[0];
        let mut optimum_mae = genotype_data_for_optimisation
            .clone()
            .estimate_expected_mae_in_aldknni(
                &optimum_min_loci_corr,
                &optimum_max_pool_dist,
                &optimum_min_l_loci,
                &optimum_min_k_neighbours,
                restrict_linked_loci_per_chromosome,
                &loci_idx,
                &vec_masked_loci_idx,
                &vec_vec_masked_alleles_freqs,
            )
            .expect("Error calling estimate_expected_mae_in_aldknni() method within optimise_params_and_estimate_accuracy().");
        if all_parameters_are_fixed == false {
            for l in 0..vec_min_k_neighbours.len() {
                let mae = genotype_data_for_optimisation
                                .clone()
                                .estimate_expected_mae_in_aldknni(
                                    &optimum_min_loci_corr,
                                    &optimum_max_pool_dist,
                                    &optimum_min_l_loci,
                                    &vec_min_k_neighbours[l],
                                    restrict_linked_loci_per_chromosome,
                                    &loci_idx,
                                    &vec_masked_loci_idx,
                                    &vec_vec_masked_alleles_freqs,
                                )
                                .expect("Error calling estimate_expected_mae_in_aldknni() method within optimise_params_and_estimate_accuracy().");
                if optimise_n_reps > 1 {
                    println!(
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        r,
                        vec_min_k_neighbours[l],
                        optimum_max_pool_dist,
                        optimum_min_l_loci,
                        optimum_min_loci_corr,
                        mae
                    );
                } else {
                    println!(
                        "{}\t{}\t{}\t{}\t{}",
                        vec_min_k_neighbours[l],
                        optimum_max_pool_dist,
                        optimum_min_l_loci,
                        optimum_min_loci_corr,
                        mae
                    );
                }
                if mae < optimum_mae {
                    optimum_min_k_neighbours = vec_min_k_neighbours[l];
                    optimum_mae = mae;
                } else if (mae - optimum_mae) > f64::EPSILON {
                    break;
                }
            }
            for j in 0..vec_max_pool_dist.len() {
                let mae = genotype_data_for_optimisation
                                .clone()
                                .estimate_expected_mae_in_aldknni(
                                    &optimum_min_loci_corr,
                                    &vec_max_pool_dist[j],
                                    &optimum_min_l_loci,
                                    &optimum_min_k_neighbours,
                                    restrict_linked_loci_per_chromosome,
                                    &loci_idx,
                                    &vec_masked_loci_idx,
                                    &vec_vec_masked_alleles_freqs,
                                )
                                .expect("Error calling estimate_expected_mae_in_aldknni() method within optimise_params_and_estimate_accuracy().");
                if optimise_n_reps > 1 {
                    println!(
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        r,
                        optimum_min_k_neighbours,
                        vec_max_pool_dist[j],
                        optimum_min_l_loci,
                        optimum_min_loci_corr,
                        mae
                    );
                } else {
                    println!(
                        "{}\t{}\t{}\t{}\t{}",
                        optimum_min_k_neighbours,
                        vec_max_pool_dist[j],
                        optimum_min_l_loci,
                        optimum_min_loci_corr,
                        mae
                    );
                }
                if mae < optimum_mae {
                    optimum_max_pool_dist = vec_max_pool_dist[j];
                    optimum_mae = mae;
                } else if (mae - optimum_mae) > f64::EPSILON {
                    break;
                }
            }
            for k in 0..vec_min_l_loci.len() {
                let mae = genotype_data_for_optimisation
                                .clone()
                                .estimate_expected_mae_in_aldknni(
                                    &optimum_min_loci_corr,
                                    &optimum_max_pool_dist,
                                    &vec_min_l_loci[k],
                                    &optimum_min_k_neighbours,
                                    restrict_linked_loci_per_chromosome,
                                    &loci_idx,
                                    &vec_masked_loci_idx,
                                    &vec_vec_masked_alleles_freqs,
                                )
                                .expect("Error calling estimate_expected_mae_in_aldknni() method within optimise_params_and_estimate_accuracy().");
                if optimise_n_reps > 1 {
                    println!(
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        r,
                        optimum_min_k_neighbours,
                        optimum_max_pool_dist,
                        vec_min_l_loci[k],
                        optimum_min_loci_corr,
                        mae
                    );
                } else {
                    println!(
                        "{}\t{}\t{}\t{}\t{}",
                        optimum_min_k_neighbours,
                        optimum_max_pool_dist,
                        vec_min_l_loci[k],
                        optimum_min_loci_corr,
                        mae
                    );
                }
                if mae < optimum_mae {
                    optimum_min_l_loci = vec_min_l_loci[k];
                    optimum_mae = mae;
                } else if (mae - optimum_mae) > f64::EPSILON {
                    break;
                }
            }
            for i in 0..vec_min_loci_corr.len() {
                let mae = genotype_data_for_optimisation
                                .clone()
                                .estimate_expected_mae_in_aldknni(
                                    &vec_min_loci_corr[i],
                                    &optimum_max_pool_dist,
                                    &optimum_min_l_loci,
                                    &optimum_min_k_neighbours,
                                    restrict_linked_loci_per_chromosome,
                                    &loci_idx,
                                    &vec_masked_loci_idx,
                                    &vec_vec_masked_alleles_freqs,
                                )
                                .expect("Error calling estimate_expected_mae_in_aldknni() method within optimise_params_and_estimate_accuracy().");
                if optimise_n_reps > 1 {
                    println!(
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        r,
                        optimum_min_k_neighbours,
                        optimum_max_pool_dist,
                        optimum_min_l_loci,
                        vec_min_loci_corr[i],
                        mae
                    );
                } else {
                    println!(
                        "{}\t{}\t{}\t{}\t{}",
                        optimum_min_k_neighbours,
                        optimum_max_pool_dist,
                        optimum_min_l_loci,
                        vec_min_loci_corr[i],
                        mae
                    );
                }
                if mae < optimum_mae {
                    optimum_min_loci_corr = vec_min_loci_corr[i];
                    optimum_mae = mae;
                } else if (mae - optimum_mae) > f64::EPSILON {
                    break;
                }
            }
        }
        // Append parameters
        mae.push(optimum_mae);
        min_loci_corr.push(optimum_min_loci_corr);
        max_pool_dist.push(optimum_max_pool_dist);
        min_l_loci.push(optimum_min_l_loci);
        min_k_neighbours.push(optimum_min_k_neighbours);
    }
    if all_parameters_are_fixed == false {
        println!("-----------------------------------------------");
    }
    Ok((
        min_loci_corr.iter().fold(0.0, |sum, &x| sum + x) / (optimise_n_reps as f64),
        max_pool_dist.iter().fold(0.0, |sum, &x| sum + x) / (optimise_n_reps as f64),
        min_l_loci.iter().fold(0, |sum, &x| sum + x) / (optimise_n_reps as u64),
        min_k_neighbours.iter().fold(0, |sum, &x| sum + x) / (optimise_n_reps as u64),
        mae.iter().fold(0.0, |sum, &x| sum + x) / (optimise_n_reps as f64),
    ))
}

// Make tests
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_optim() {
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
        println!("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
        let n_threads = 2;
        let mut frequencies_and_phenotypes = file_sync_phen
            .into_genotypes_and_phenotypes(&filter_stats, false, &n_threads)
            .unwrap();
        println!(
            "Before simulating missing data:\n{:?}",
            frequencies_and_phenotypes.intercept_and_allele_frequencies
        );
        println!("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
        let (
            _,
            _max_sparsity,
            loci_idx,
            vec_masked_loci_idx,
            _vec_masked_coverages,
            vec_vec_masked_alleles_freqs,
        ) = frequencies_and_phenotypes.simulate_missing(&0.25).unwrap();
        println!(
            "After simulating missing data:\n{:?}",
            frequencies_and_phenotypes.intercept_and_allele_frequencies
        );
        println!("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
        let vec_vec_imputed_mask = frequencies_and_phenotypes
            .extract_imputed_mask(&loci_idx, &vec_masked_loci_idx)
            .unwrap();

        println!("loci_idx.len()={}", loci_idx.len());
        println!("vec_masked_loci_idx.len()={}", vec_masked_loci_idx.len());
        println!(
            "vec_vec_masked_alleles_freqs.len()={}",
            vec_vec_masked_alleles_freqs.len()
        );
        println!("vec_vec_imputed_mask.len()={}", vec_vec_imputed_mask.len());

        println!(
            "vec_masked_loci_idx[vec_masked_loci_idx.len()-1]={}",
            vec_masked_loci_idx[vec_masked_loci_idx.len() - 1]
        );
        println!(
            "vec_vec_masked_alleles_freqs[vec_vec_masked_alleles_freqs.len()-1]={:?}",
            vec_vec_masked_alleles_freqs[vec_vec_masked_alleles_freqs.len() - 1]
        );

        println!("After simlating missing data:");
        println!(
            "vec_vec_imputed_mask[0..10]={:?}",
            &vec_vec_imputed_mask[0..10]
        );
        println!(
            "vec_vec_imputed_mask[vec_vec_imputed_mask.len()-1]={:?}",
            vec_vec_imputed_mask[vec_vec_imputed_mask.len() - 1]
        );
        assert!(vec_vec_imputed_mask[0][0].is_nan());
        assert!(vec_vec_imputed_mask[vec_vec_imputed_mask.len() - 1][0].is_nan());
        let _vec_vec_imputed_mask = frequencies_and_phenotypes
            .extract_imputed_mask(&loci_idx, &vec_masked_loci_idx)
            .unwrap();
        // assert_eq!(0, 1);
    }
}
