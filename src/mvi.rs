use crate::helpers::*;
use crate::structs_and_traits::*;
use ndarray::prelude::*;
use rand::prelude::IteratorRandom;
use rand::rngs::StdRng;
use rand::SeedableRng;
use std::io;

type SparsitySimulationOutput = (f64, Vec<usize>, Vec<f64>, Vec<Vec<f64>>);

impl GenotypesAndPhenotypes {
    pub fn mean_imputation(&mut self) -> io::Result<&mut Self> {
        self.check().expect("Error calling check() method within mean_imputation() method for GenotypesAndPhenotypes struct.");
        // We are assuming that all non-zero alleles across pools are kept, i.e. biallelic loci have 2 columns, triallelic have 3, and so on.
        let (n, p) = self.intercept_and_allele_frequencies.dim();
        let (loci_idx, _loci_chr, _loci_pos) = self.count_loci().expect("Error defining loci indexes and identities via count_loci() method within mean_imputation() method for GenotypesAndPhenotypes struct.");
        let l = loci_idx.len() - 1; // less the final position, this is equivalent to the the number of columns of self.coverages
        for j in 0..l {
            // Use the indexes of each locus
            let idx_ini = loci_idx[j];
            let idx_fin = if j < (l - 1) { loci_idx[j + 1] } else { p };
            let freqs: Array2<f64> = self
                .intercept_and_allele_frequencies
                .slice(s![.., idx_ini..idx_fin])
                .to_owned();
            let mean_freqs = mean_axis_ignore_nan(&freqs, 0).expect("Error calculating axis-wise mean ignoring NANs within mean_imputation() method for GenotypesAndPhenotypes struct.");
            let sum_freqs = mean_freqs.sum();
            // We need to correct for imputations resulting in a sum of allele frequencies greater or less than 1
            let mean_freqs = if sum_freqs != 1.0 {
                mean_freqs.map(|x| x / sum_freqs)
            } else {
                mean_freqs
            };
            for k in 0..freqs.ncols() {
                for i in 0..n {
                    if self.intercept_and_allele_frequencies[(i, idx_ini + k)].is_nan() {
                        self.intercept_and_allele_frequencies[(i, idx_ini + k)] = mean_freqs[k];
                    } else {
                        continue;
                    };
                }
            }
        }
        // Set missing coverages to infinity to mark imputed data
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
        // println!("self.coverages = {:?}", self.coverages);
        Ok(self)
    }

    pub fn simulate_missing(
        &mut self,
        loci_idx: &[usize],
        max_sparsity: &f64,
        rep: &usize,
    ) -> io::Result<(&mut Self, SparsitySimulationOutput)> {
        self.check().expect("Error calling check() method within simulate_missing() method for GenotypesAndPhenotypes struct.");
        let (n, l) = self.coverages.dim();
        // let (loci_idx, _loci_chr, _loci_pos) = self.count_loci().expect("Error defining loci indexes and identities via count_loci() method within simulate_missing() method for GenotypesAndPhenotypes struct.");
        // Limit max_sparsity to 90%
        let max_sparsity = if *max_sparsity > 0.9 {
            0.9
        } else {
            *max_sparsity
        };
        // Define the number of loci to mask
        let n_masked = (((n * l) as f64) * max_sparsity).floor() as usize;
        // Sample random loci
        // let mut rng = rand::thread_rng();
        let mut rng = StdRng::seed_from_u64(*rep as u64);
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
            (
                max_sparsity,
                vec_masked_loci_idx,
                vec_masked_coverages,
                vec_vec_masked_alleles_freqs,
            ),
        ))
    }

    pub fn extract_imputed_mask(
        &self,
        loci_idx: &[usize],
        vec_masked_loci_idx: &[usize],
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

    pub fn estimate_expected_mae_in_mvi(&self) -> io::Result<f64> {
        let mut genotype_data_for_optimisation = self.clone();
        let (n, p) = genotype_data_for_optimisation
            .intercept_and_allele_frequencies
            .dim();
        let missing_rate_sim = if (n * p) < 20_000 {
            0.01
        } else {
            10_000.0 / ((n * p) as f64)
        };
        // Extract loci indices
        let (loci_idx, _loci_chr, _loci_pos) = self.count_loci().expect("Error calling count_loci() method within adaptive_ld_knn_imputation() method for GenotypesAndPhenotypes struct.");
        let (
            _,
            (max_sparsity,
            vec_masked_loci_idx,
            _vec_masked_coverages,
            vec_vec_masked_alleles_freqs),
        ) = genotype_data_for_optimisation
            .simulate_missing(&loci_idx, &missing_rate_sim, &0)
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
                        || (vec_vec_imputed_mask[i][j].is_nan())
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

/// # `mvi`: mean value imputation
///
/// This imputation uses the arithmetic mean of the observed allele frequencies across all samples where the locus was genotyped:
///
/// $$
/// \hat q_{r,j} = { {1 \over (n-m)} { \sum_{i \ne r}^{n} q_{i,j} } }
/// $$
///
/// where:
/// - $\hat q_{r,j}$ is the imputed allele frequency of sample $r$ at the $j^{\text {th}}$ locus,
/// - $n$ is the total number of samples,
/// - $m$ is the number of samples which are missing data at the $j^{\text {th}}$ locus, and
/// - $q_{i,j}$ is the known allele frequency of the $i^{\text {th}}$ sample at the $j^{\text {th}}$ locus.
pub fn impute_mean(
    mut genotypes_and_phenotypes: GenotypesAndPhenotypes,
    filter_stats: &FilterStats,
    n_threads: &usize,
    out: &str,
) -> io::Result<String> {
    // Estimate predicted imputation accuracy
    let mae = genotypes_and_phenotypes
        .estimate_expected_mae_in_mvi()
        .expect("Error calling estimate_expected_mae_in_mvi() method within impute_mean().");
    println!("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    println!(
        "Expected imputation accuracy in terms of mean absolute error: {}",
        mae
    );
    println!("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .mean_imputation()
        .expect("Error calling mean_imputation() method within impute_mean().");
    let end = std::time::SystemTime::now();
    let duration = end
        .duration_since(start)
        .expect("Error measuring the duration of mean value imputation within impute_mean()");
    println!(
        "Mean value imputation: {} pools x {} loci | Missingness: {}% | Duration: {} seconds",
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        genotypes_and_phenotypes
            .missing_rate()
            .expect("Error measuring sparsity after mean value imputation within impute_mean()."),
        duration.as_secs_f64()
    );

    // println!("genotypes_and_phenotypes={:?}", genotypes_and_phenotypes);

    // Remove 100% of the loci with missing data
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .filter_out_top_missing_loci(&1.00)
        .expect("Error calling filter_out_top_missing_loci() method within impute_mean().");
    let end = std::time::SystemTime::now();
    let duration = end.duration_since(start).expect(
        "Error measuring the duration of filter_out_top_missing_loci() within impute_mean()",
    );
    println!(
        "Missing data removed, i.e. loci which cannot be imputed because of extreme sparsity: {} pools x {} loci | Missingness: {}% | Duration: {} seconds",
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        genotypes_and_phenotypes.missing_rate().expect("Error measuring sparsity after filter_out_top_missing_loci(0) within impute_mean()."),
        duration.as_secs_f64()
    );
    // Output
    let out = genotypes_and_phenotypes
        .write_tsv(filter_stats, false, out, n_threads)
        .expect("Error writing the output of mean value imputation within impute_mean().");

    Ok(out)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_mvi() {
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
        let max_depth_above_which_are_missing = 10.0;
        let _ = frequencies_and_phenotypes
            .set_missing_by_depth(
                &min_depth_below_which_are_missing,
                &max_depth_above_which_are_missing,
            )
            .unwrap();
        println!(
            "Before imputation:\n{:?}",
            frequencies_and_phenotypes.intercept_and_allele_frequencies
        );
        let _ = frequencies_and_phenotypes.mean_imputation().unwrap();
        println!(
            "After imputation:\n{:?}",
            frequencies_and_phenotypes.intercept_and_allele_frequencies
        );

        let keep_p_minus_1 = false;
        let _start = std::time::SystemTime::now();
        let genotypes_and_phenotypes = file_sync_phen
            .convert_into_genotypes_and_phenotypes(&filter_stats, keep_p_minus_1, &n_threads)
            .unwrap();

        let outname = impute_mean(
            genotypes_and_phenotypes,
            &filter_stats,
            &n_threads,
            &"test-impute_mean.tsv".to_owned(),
        )
        .unwrap();
        assert_eq!(outname, "test-impute_mean.tsv".to_owned()); // Do better!!! Load data - thus working on improving load_table()

        assert_eq!(
            frequencies_and_phenotypes.intercept_and_allele_frequencies[(0, 1)],
            Array1::from_vec(vec![0.3333333333333333, 0.2, 0.14285714285714285])
                .mean()
                .unwrap()
        );
        assert_eq!(
            frequencies_and_phenotypes.intercept_and_allele_frequencies[(1, 1)],
            Array1::from_vec(vec![0.3333333333333333, 0.2, 0.14285714285714285])
                .mean()
                .unwrap()
        );
        assert_eq!(
            frequencies_and_phenotypes.intercept_and_allele_frequencies[(0, 2)],
            Array1::from_vec(vec![0.6666666666666666, 0.8, 0.8571428571428571])
                .mean()
                .unwrap()
        );
        assert_eq!(
            frequencies_and_phenotypes.intercept_and_allele_frequencies[(1, 2)],
            Array1::from_vec(vec![0.6666666666666666, 0.8, 0.8571428571428571])
                .mean()
                .unwrap()
        );
        // assert_eq!(0, 1);
    }
}
