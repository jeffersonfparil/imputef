use ndarray::prelude::*;
use std::io;

use crate::helpers::*;
use crate::structs_and_traits::*;

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
}

// Impute using mean allele frequencies across pools
pub fn impute_mean(
    mut genotypes_and_phenotypes: GenotypesAndPhenotypes,
    filter_stats: &FilterStats,
    _min_depth_below_which_are_missing: &f64,
    _max_depth_above_which_are_missing: &f64,
    _frac_top_missing_pools: &f64,
    _frac_top_missing_loci: &f64,
    n_threads: &usize,
    out: &String,
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
        duration.as_secs()
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
        duration.as_secs()
    );
    // Output
    let out = genotypes_and_phenotypes
        .write_csv(filter_stats, false, out, n_threads)
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
            .into_genotypes_and_phenotypes(&filter_stats, false, &n_threads)
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
            .into_genotypes_and_phenotypes(&filter_stats, keep_p_minus_1, &n_threads)
            .unwrap();

        let outname = impute_mean(
            genotypes_and_phenotypes,
            &filter_stats,
            &min_depth_below_which_are_missing,
            &max_depth_above_which_are_missing,
            &0.2,
            &0.5,
            &n_threads,
            &"test-impute_mean.csv".to_owned(),
        )
        .unwrap();
        assert_eq!(outname, "test-impute_mean.csv".to_owned()); // Do better!!! Load data - thus working on improving load_table()

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

// SYNC=/home/jeff/poolgen/tests/test_REMOVE_ME_BEFORE_PUSHING.sync
// PHEN=/home/jeff/poolgen/tests/test_REMOVE_ME_BEFORE_PUSHING.csv
// NCORES=7
// OUT=/home/jeff/poolgen/tests/test-MEAN_IMPUTE-REMOVE_ME_BEFORE_PUSHING.csv
// time cargo run -- impute \
//     --imputation-method "mean" \
//     -f ${SYNC} \
//     -p ${PHEN} \
//     --phen-delim , \
//     --phen-name-col 0 \
//     --phen-pool-size-col 1 \
//     --phen-value-col 2 \
//     --min-allele-frequency 0.0001 \
//     --min-coverage 0 \
//     --min-quality 0.01 \
//     --max-missingness-rate 0.75 \
//     --min-depth-set-to-missing 5 \
//     --window-size-bp 1000000 \
//     --window-slide-size-bp 1000000 \
//     --min-loci-per-window 1 \
//     --min-correlation 0.5 \
//     --k-neighbours 5 \
//     --n-threads ${NCORES} \
//     -o ${OUT}
