use ndarray::prelude::*;
use rand::prelude::IteratorRandom;
use std::io;

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
        self.check().unwrap();
        let (n, l) = self.coverages.dim();
        let (loci_idx, _loci_chr, _loci_pos) = self.count_loci().unwrap();
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
        // println!("loci_idx={:?}", loci_idx);
        // println!("n_masked={:?}", n_masked);
        // println!("vec_masked_loci_idx_tmp[0..10]={:?}", &vec_masked_loci_idx_tmp[0..10]);
        // println!("vec_masked_loci_idx_tmp[100..110]={:?}", &vec_masked_loci_idx_tmp[100..110]);
        // println!("vec_masked_loci_idx.len()={:?}", vec_masked_loci_idx.len());
        // println!("vec_masked_coverages.len()={:?}", vec_masked_coverages.len());
        // println!("vec_vec_masked_alleles_freqs.len()={:?}", vec_vec_masked_alleles_freqs.len());
        self.check().unwrap();
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
                        // if self.intercept_and_allele_frequencies[(i, k)].is_nan() {
                        //     println!("self.intercept_and_allele_frequencies={:?}", self.intercept_and_allele_frequencies);
                        // }
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
        window_size_bp: &u64,
        window_slide_size_bp: &u64,
        min_loci_per_window: &u64,
        min_loci_corr: &f64,
        max_pool_dist: &f64,
        loci_idx: &Vec<usize>,
        vec_masked_loci_idx: &Vec<usize>,
        vec_vec_masked_alleles_freqs: &Vec<Vec<f64>>,
    ) -> io::Result<f64> {
        // println!("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
        // let q = genotype_data_for_optimisation.intercept_and_allele_frequencies.fold_axis(Axis(0), 0, |&sum, &x| if x.is_nan() {sum + 1} else {sum})
        //     .iter()
        //     .fold(0, |max, &x| if x > max {x} else {max});
        // let r = genotype_data_for_optimisation.intercept_and_allele_frequencies.fold_axis(Axis(1), 0, |&sum, &x| if x.is_nan() {sum + 1} else {sum})
        //     .iter()
        //     .fold(0, |max, &x| if x > max {x} else {max});
        // println!("q={:?}; r={:?}", q, r);
        // println!("genotype_data_for_optimisation.intercept_and_allele_frequencies={:?}", genotype_data_for_optimisation.intercept_and_allele_frequencies);
        self.adaptive_ld_knn_imputation(
            window_size_bp,
            window_slide_size_bp,
            min_loci_per_window,
            min_loci_corr,
            max_pool_dist,
            false
        )
        .unwrap();
        // Extract imputed frequencies
        let vec_vec_imputed_mask = self
            .extract_imputed_mask(&loci_idx, &vec_masked_loci_idx)
            .unwrap();
        let m = vec_vec_masked_alleles_freqs.len(); // total number of loci
        let mean_absolute_error = if m == 0 {
            // If no loci were simulated to be missing
            1.00
        } else {
            // println!(
            //     "vec_vec_masked_alleles_freqs={:?}",
            //     vec_vec_masked_alleles_freqs
            // );
            // println!("vec_vec_imputed_mask={:?}", vec_vec_imputed_mask);
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
                            // (vec_vec_masked_alleles_freqs[i][j] - vec_vec_imputed_mask[i][j]).abs();
                            (vec_vec_masked_alleles_freqs[i][j] - vec_vec_imputed_mask[i][j]).powf(2.0);
                        a += 1.00;
                    }
                }
            }
            if a == 0.0 {
                1.00
            } else {
                // sum_abs / a
                (sum_abs / a).sqrt()
            }
        };
        // println!("mean_absolute_error={:?}", mean_absolute_error);
        Ok(mean_absolute_error)
    }

    pub fn estimate_expected_mae_in_mvi(&self) -> io::Result<f64> {
        let mut genotype_data_for_optimisation = self.clone();
        let (
            _,
            max_sparsity,
            loci_idx,
            vec_masked_loci_idx,
            _vec_masked_coverages,
            vec_vec_masked_alleles_freqs,
        ) = genotype_data_for_optimisation
            .simulate_missing(&0.1)
            .unwrap();
        assert_eq!(
            max_sparsity, 0.1,
            "Ooops! Missing unexpected simulated sparsity!"
        );
        let _ = genotype_data_for_optimisation.mean_imputation().unwrap();
        let vec_vec_imputed_mask = genotype_data_for_optimisation
            .extract_imputed_mask(&loci_idx, &vec_masked_loci_idx)
            .unwrap();
        // Compare imputed and expected frequencies
        let m = vec_vec_masked_alleles_freqs.len(); // total number of loci
        let mean_absolute_error = if m == 0 {
            // If no loci we're simulated to be missing
            f64::INFINITY
        } else {
            let mut sum_abs = 0.0;
            let mut a = 0.0; // total number of alleles across all loci
            for i in 0..m {
                let w = vec_vec_masked_alleles_freqs[i].len();
                for j in 0..w {
                    // sum_abs += (vec_vec_masked_alleles_freqs[i][j] - vec_vec_imputed_mask[i][j]).abs();
                    sum_abs += (vec_vec_masked_alleles_freqs[i][j] - vec_vec_imputed_mask[i][j]).powf(2.0);
                    a += 1.00;
                }
            }
            // sum_abs / a
            (sum_abs / a).sqrt()
        };
        Ok(mean_absolute_error)
    }
}

pub fn optimise_params_and_estimate_accuracy(
    genotypes_and_phenotypes: &GenotypesAndPhenotypes,
    min_loci_corr: &f64,
    max_pool_dist: &f64,
    optimise_for_thresholds: &bool,
    optimise_n_steps_corr: &usize,
    optimise_n_steps_dist: &usize,
    optimise_n_reps: &usize,
    window_size_bp: &u64,
    window_slide_size_bp: &u64,
    min_loci_per_window: &u64,
) -> io::Result<(f64, f64, f64)> {
    let (min_loci_corr, max_pool_dist, mae) = if (*min_loci_corr == 0.0) | (*max_pool_dist == 0.0) {
        // If we are optimising for l or corr and k or dist
        let (vec_min_corr, vec_max_dist): (Vec<f64>, Vec<f64>) = if *optimise_for_thresholds {
            println!("Optimising for the optimum minimum loci correlation and maximum genetic distance to identify the loci in LD to be used in genetic distance estimation and the nearest neighbours for weighted mean allele frequency calculation.");
            let vec_min_corr: Vec<f64> = if *min_loci_corr == 0.0 {
                (10..=100).step_by(*optimise_n_steps_corr).map(|x| x as f64 / 100.0).collect()
            } else {
                vec![*min_loci_corr]
            };
            let vec_max_dist: Vec<f64> = if *max_pool_dist == 0.0 {
                (10..=100).step_by(*optimise_n_steps_dist).map(|x| x as f64 / 100.0).collect()
            } else {
                vec![*max_pool_dist]
            };
            println!("Number of parameter levels to test | min_cor={} | max_dist={} | reps={}", vec_min_corr.len(), vec_max_dist.len(), *optimise_n_reps);
            (vec_min_corr, vec_max_dist)
        } else {
            println!("Optimising for the optimum number of loci to include in genetic distance estimation and number of nearest neighbours to include in weighted mean allele frequency calculation.");
            let n: usize = genotypes_and_phenotypes.coverages.nrows();
            let step_size_l = if (n - 1) < *optimise_n_steps_corr {
                1
            } else {
                ((n as f64) / *optimise_n_steps_corr as f64).ceil() as usize
            };
            let step_size_k = if (n - 1) < *optimise_n_steps_dist {
                1
            } else {
                ((n as f64) / *optimise_n_steps_dist as f64).ceil() as usize
            };
            let vec_l: Vec<f64> = if *min_loci_corr == 0.0 {
                (1..=n)
                    .step_by(step_size_l)
                    .map(|x| x as f64)
                    .collect()
            } else {
                vec![*min_loci_corr]
            };
            let vec_k: Vec<f64> = if *max_pool_dist == 0.0 {
                (1..=n)
                    .step_by(step_size_k)
                    .map(|x| x as f64)
                    .collect()
            } else {
                vec![*max_pool_dist]
            };
            println!("Number of parameter levels to test | l-loci={} | k-neighbours={} | reps={}", vec_l.len(), vec_k.len(), *optimise_n_reps);
            (vec_l, vec_k)
        };
        // println!("vec_min_corr={:?}", vec_min_corr);
        // println!("vec_max_dist={:?}", vec_max_dist);
        // let optimise_n_reps = 10;
        println!("-----------------------------------------");
        if *optimise_n_reps > 1 {
            if *optimise_for_thresholds {
                println!("rep\tmin_corr\tmax_dist\tmae");
            } else {
                println!("rep\tl-loci\tk-neighbours\tmae");
            }
        } else {
            if *optimise_for_thresholds {
                println!("min_corr\tmax_dist\tmae");
            } else {
                println!("l-loci\tk-neighbours\tmae");
            }
        }
        // Use the same simulated missing data across the range of corr/l and dist/k...
        let mut array3_mae: Array3<f64> = Array3::from_elem(
            (*optimise_n_reps, vec_min_corr.len(), vec_max_dist.len()),
            f64::NAN,
        );
        for r in 0..*optimise_n_reps {
            // Simulate 10% sparsity to determine accuracy (in terms of MAE) and to optimise for the best k and l, if applicable, i.e. n_loci_to_estimate_distance==0 or k_neighbours==0
            let mut genotype_data_for_optimisation = genotypes_and_phenotypes.clone();
            let (
                _,
                max_sparsity,
                loci_idx,
                vec_masked_loci_idx,
                _vec_masked_coverages,
                vec_vec_masked_alleles_freqs,
            ) = genotype_data_for_optimisation
                .simulate_missing(&0.1)
                .unwrap();
            assert_eq!(
                max_sparsity, 0.1,
                "Ooops! Missing unexpected simulated sparsity!"
            );
            for i in 0..vec_min_corr.len() {
                for j in 0..vec_max_dist.len() {
                    let mae = genotype_data_for_optimisation
                        .clone()
                        .estimate_expected_mae_in_aldknni(
                            window_size_bp,
                            window_slide_size_bp,
                            min_loci_per_window,
                            &vec_min_corr[i],
                            &vec_max_dist[j],
                            &loci_idx,
                            &vec_masked_loci_idx,
                            &vec_vec_masked_alleles_freqs,
                        )
                        .unwrap();
                    array3_mae[(r, i, j)] = mae;
                    if *optimise_n_reps > 1 {
                        println!("{}\t{}\t{}\t{}", r, vec_min_corr[i], vec_max_dist[j], mae);
                    } else {
                        println!("{}\t{}\t{}", vec_min_corr[i], vec_max_dist[j], mae);
                    }
                }
            }
        }
        println!("-----------------------------------------");
        // Calculate means per corr or l x dist or k pair
        if *optimise_n_reps > 1 {
            println!("-----------------------------------------");
            if *optimise_for_thresholds {
                println!("min_corr\tmax_dist\tmae");
            } else {
                println!("l-loci\tk-neighbours\tmae");
            }
        }
        let t = vec_min_corr.len() * vec_max_dist.len();
        let mut vec_corr: Vec<f64> = Vec::with_capacity(t);
        let mut vec_dist: Vec<f64> = Vec::with_capacity(t);
        let mut vec_mae: Vec<f64> = Vec::with_capacity(t);
        let (mut optim_corr, mut optim_dist, mut optim_mae) = (0.0, 0.0, 1.0);
        for i in 0..vec_min_corr.len() {
            for j in 0..vec_max_dist.len() {
                let (mut mae, mut n) = (0.0, 0.0);
                for r in 0..*optimise_n_reps {
                    if array3_mae[(r, i, j)].is_nan() == false {
                        mae += array3_mae[(r, i, j)];
                        n += 1.0;
                    }
                }
                mae = mae / n;
                // variance across reps
                let mut var = 0.0;
                if *optimise_n_reps > 1 {
                    for r in 0..*optimise_n_reps {
                        var = var
                            + if array3_mae[(r, i, j)].is_nan() == false {
                                (array3_mae[(r, i, j)] - mae).powf(2.0)
                            } else {
                                0.0
                            };
                    }
                    var = var / (n - 1.00);
                    println!("{}\t{}\t{}\t{}", vec_min_corr[i], vec_max_dist[j], mae, var);
                }

                (optim_corr, optim_dist, optim_mae) = if mae < optim_mae {
                    (vec_min_corr[i], vec_max_dist[j], mae)
                } else {
                    (optim_corr, optim_dist, optim_mae)
                };
                vec_corr.push(vec_min_corr[i]);
                vec_dist.push(vec_max_dist[j]);
                vec_mae.push(mae);
            }
        }
        if *optimise_n_reps > 1 {
            println!("-----------------------------------------");
        }
        // let n = vec_mae.iter().fold(0.0, |n, &x| if !x.is_nan() {n + 1.00} else {n});
        // let mu = vec_mae.iter().fold(0.0, |mu, &x| if !x.is_nan() {mu + x} else {mu}) / n;
        // let sd = (vec_mae.iter().fold(0.0, |ss, &x| if !x.is_nan() {ss + (x - mu).powf(2.0)} else {ss}) / (n - 1.0)).sqrt();
        // println!("n={:?}; mu={:?}; sd={:?}", n, mu, sd);
        (optim_corr, optim_dist, optim_mae)
    } else {
        // If we are not optimising for l or corr and k or dist
        // Estimate predicted imputation accuracy
        let mut genotype_data_for_expected_mae = genotypes_and_phenotypes.clone();
        let (
            _,
            max_sparsity,
            loci_idx,
            vec_masked_loci_idx,
            _vec_masked_coverages,
            vec_vec_masked_alleles_freqs,
        ) = genotype_data_for_expected_mae
            .simulate_missing(&0.1)
            .unwrap();
        assert_eq!(
            max_sparsity, 0.1,
            "Ooops! Missing unexpected simulated sparsity!"
        );
        let mae = genotype_data_for_expected_mae
            .estimate_expected_mae_in_aldknni(
                window_size_bp,
                window_slide_size_bp,
                min_loci_per_window,
                &*min_loci_corr,
                &max_pool_dist,
                &loci_idx,
                &vec_masked_loci_idx,
                &vec_vec_masked_alleles_freqs,
            )
            .unwrap();
        (*min_loci_corr, *max_pool_dist, mae)
    };
    Ok((min_loci_corr, max_pool_dist, mae))
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
