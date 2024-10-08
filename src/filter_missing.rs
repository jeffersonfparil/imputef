use ndarray::prelude::*;

use crate::helpers::*;
use crate::structs_and_traits::*;

impl GenotypesAndPhenotypes {
    pub fn set_missing_by_depth(
        &mut self,
        min_depth_below_which_are_missing: &f64,
        max_depth_above_which_are_missing: &f64,
    ) -> Result<&mut Self, ImputefError> {
        match self.check() {
            Ok(x) => x,
            Err(e) => return Err(ImputefError {
                code: 201,
                message:
                    "Error checking GenotypesAndPhenotypes in the method set_missing_by_depth() | "
                        .to_owned()
                        + &e.message,
            }),
        };
        let (n, _p) = self.intercept_and_allele_frequencies.dim();
        let (_n, l) = self.coverages.dim();
        let (loci_idx, _loci_chr, _loci_pos) = self.count_loci().expect("Error defining loci indexes and identities via count_loci() method within set_missing_by_depth() method for GenotypesAndPhenotypes struct.");
        for i in 0..n {
            for j in 0..l {
                if (self.coverages[(i, j)] < *min_depth_below_which_are_missing)
                    || (self.coverages[(i, j)] > *max_depth_above_which_are_missing)
                {
                    self.coverages[(i, j)] = f64::NAN;
                    // Use the indexes of the locus to set missing values to all alleles in the locus
                    let idx_ini = loci_idx[j];
                    let idx_fin = loci_idx[j + 1];
                    for k in idx_ini..idx_fin {
                        self.intercept_and_allele_frequencies[(i, k)] = f64::NAN;
                    }
                }
            }
        }
        match self.check() {
            Ok(x) => x,
            Err(e) => return Err(ImputefError {
                code: 202,
                message:
                    "Error checking GenotypesAndPhenotypes in the method set_missing_by_depth() | "
                        .to_owned()
                        + &e.message,
            }),
        };
        Ok(self)
    }

    pub fn missing_rate(&mut self) -> Result<f64, ImputefError> {
        let (n, _p) = self.intercept_and_allele_frequencies.dim();
        let (_n, l) = self.coverages.dim();
        let sum: f64 = self.coverages.fold(0.00, |sum, &x| {
            if (x.is_nan()) || (x == 0.0) {
                sum + 1.00
            } else {
                sum
            }
        });
        sensible_round(sum * 100.0 / ((n * l) as f64), 5)
    }

    pub fn filter_out_top_missing_pools(
        &mut self,
        frac_top_missing_pools: &f64,
    ) -> Result<&mut Self, ImputefError> {
        match self.check() {
            Ok(x) => x,
            Err(e) => return Err(ImputefError{
                code: 203,
                message: "Error checking GenotypesAndPhenotypes in the method filter_out_top_missing_pools() | ".to_owned() +
                &e.message
            })
        };
        let n = self.intercept_and_allele_frequencies.nrows();
        let p = self.intercept_and_allele_frequencies.ncols() - 1;
        let missingness_per_pool: Array1<f64> = self
            .intercept_and_allele_frequencies
            .map_axis(Axis(1), |x| {
                x.fold(
                    0_u64,
                    |n_missing, &x| if x.is_nan() { n_missing + 1 } else { n_missing },
                )
            })
            .map(|&x| x as f64 / p as f64);
        // println!("missingness_per_pool={:?}", missingness_per_pool);
        // Define the total number of pools that will be retained
        let n_missing =
            missingness_per_pool.fold(0.0, |sum, &x| if x > 0.0 { sum + 1.0 } else { sum });
        let n_after_filtering = n - (n_missing * frac_top_missing_pools).ceil() as usize;
        if n_after_filtering == 0 {
            return Err(ImputefError {
                code: 204,
                message: "No pools left after filtering, please reduce 'frac_top_missing_pools'"
                    .to_owned(),
            });
        }
        // Sort by increasing missingness
        let mut idx = (0..n).collect::<Vec<usize>>();
        idx.sort_by(|&a, &b| {
            missingness_per_pool[a]
                .partial_cmp(&missingness_per_pool[b])
                .expect("Error sorting indexes by sparsity  within filter_out_top_missing_pools() method for GenotypesAndPhenotypes struct.")
        });
        // println!("idx={:?}", idx);
        // Omit the pools with high missingness rates
        idx = idx[0..n_after_filtering].to_vec();
        // println!("idx={:?}", idx);
        // Sort the indexes so that we maintain the order of the pools
        idx.sort();
        // println!("idx={:?}", idx);
        let mut new_intercept_and_allele_frequencies =
            Array2::from_elem((n_after_filtering, p + 1), f64::NAN);
        let mut new_phenotypes: Array2<f64> =
            Array2::from_elem((n_after_filtering, self.phenotypes.ncols()), f64::NAN);
        let mut new_pool_names: Vec<String> = vec![];
        let mut new_coverages: Array2<f64> =
            Array2::from_elem((n_after_filtering, self.coverages.ncols()), f64::NAN);
        for i in 0..n_after_filtering {
            // Pool index in the old data
            let i_ = idx[i];
            // Intercept
            new_intercept_and_allele_frequencies[(i, 0)] = 1.0;
            for j in 1..(p + 1) {
                new_intercept_and_allele_frequencies[(i, j)] =
                    self.intercept_and_allele_frequencies[(i_, j)];
            }
            for k in 0..new_phenotypes.ncols() {
                new_phenotypes[(i, k)] = self.phenotypes[(i_, k)];
            }
            new_pool_names.push(self.pool_names[i_].to_owned());
            for l in 0..new_coverages.ncols() {
                new_coverages[(i, l)] = self.coverages[(i_, l)];
            }
        }
        self.intercept_and_allele_frequencies = new_intercept_and_allele_frequencies;
        self.phenotypes = new_phenotypes;
        self.pool_names = new_pool_names;
        self.coverages = new_coverages;
        // println!("self.chromosome.len()={:?}", self.chromosome.len());
        // println!("self.position.len()={:?}", self.position.len());
        // println!("self.allele.len()={:?}", self.allele.len());
        // println!("self.intercept_and_allele_frequencies.dim()={:?}", self.intercept_and_allele_frequencies.dim());
        // println!("self.coverages.len()={:?}", self.coverages.dim());
        match self.check() {
            Ok(x) => x,
            Err(e) => return Err(ImputefError{
                code: 205,
                message: "Error checking GenotypesAndPhenotypes in the method filter_out_top_missing_pools() | ".to_owned() +
                &e.message
            })
        };
        Ok(self)
    }

    // For thinning the data before imputation and also removing completely missing loci after imputation
    pub fn filter_out_top_missing_loci(
        &mut self,
        frac_top_missing_loci: &f64,
    ) -> Result<&mut Self, ImputefError> {
        match self.check() {
            Ok(x) => x,
            Err(e) => return Err(ImputefError{
                code: 206,
                message: "Error checking GenotypesAndPhenotypes in the method filter_out_top_missing_loci() | ".to_owned() +
                &e.message
            })
        };
        let n = self.intercept_and_allele_frequencies.nrows();
        let (loci_idx, _loci_chr, _loci_pos) = self.count_loci().expect("Error calling count_loci() method within filter_out_top_missing_loci() method for GenotypesAndPhenotypes struct.");
        let l = loci_idx.len() - 1; // Less one for the trailing locus
                                    // Define loci start and end indexes
        let loci_idx_ini = loci_idx[0..l].to_vec();
        let loci_idx_fin = loci_idx[1..(l + 1)].to_vec();
        let missingness_per_locus: Array1<f64> = self
            .coverages
            .map_axis(Axis(0), |x| {
                x.fold(
                    0_u64,
                    |n_missing, &x| if x.is_nan() { n_missing + 1 } else { n_missing },
                )
            })
            .map(|&x| x as f64 / n as f64);
        // println!("missingness_per_locus={:?}", missingness_per_locus);
        // Define the total number of pools that will be retained
        let l_missing =
            missingness_per_locus.fold(0.0, |sum, &x| if x > 0.0 { sum + 1.0 } else { sum });
        // Define the number of loci kept after filtering
        let l_after_filtering = l - (l_missing * frac_top_missing_loci).ceil() as usize;
        if l_after_filtering == 0 {
            return Err(ImputefError {
                code: 207,
                message: "No loci left after filtering, please reduce 'frac_top_missing_loci'"
                    .to_owned(),
            });
        }
        // Sort by increasing missingness
        let mut idx = (0..l).collect::<Vec<usize>>();
        idx.sort_by(|&a, &b| {
            missingness_per_locus[a]
                .partial_cmp(&missingness_per_locus[b])
                .expect("Error sorting indexes by sparsity  within filter_out_top_missing_loci() method for GenotypesAndPhenotypes struct.")
        });
        // println!("idx={:?}", idx);
        // Omit the loci with high missingness rates
        idx = idx[0..l_after_filtering].to_vec();
        // println!("idx={:?}", idx);
        // Sort the indexes so that we maintain the order of the loci
        idx.sort();
        // println!("idx={:?}", idx);
        // Find the number of columns first using the first pool and also set the loci coordinates
        let mut vec_intercept_and_allele_frequencies_pool_0: Vec<f64> =
            vec![self.intercept_and_allele_frequencies[(0, 0)]];
        let mut new_chromosome: Vec<String> = vec![self.chromosome[0].to_owned()];
        let mut new_position: Vec<u64> = vec![self.position[0]];
        let mut new_allele: Vec<String> = vec![self.allele[0].to_owned()];
        let i = 0;
        for j in idx.clone().into_iter() {
            // Use the indexes of the locus to set missing values to all alleles in the locus
            let idx_ini = loci_idx_ini[j];
            let idx_fin = loci_idx_fin[j];
            for k in idx_ini..idx_fin {
                vec_intercept_and_allele_frequencies_pool_0
                    .push(self.intercept_and_allele_frequencies[(i, k)]);
                new_chromosome.push(self.chromosome[k].to_owned());
                new_position.push(self.position[k]);
                new_allele.push(self.allele[k].to_owned());
            }
        }
        // Populate the intercept and allele frequencies matrix
        let p = vec_intercept_and_allele_frequencies_pool_0.len(); // includes the intercept
        let mut new_intercept_and_allele_frequencies = Array2::from_elem((n, p), f64::NAN);
        let mut new_coverages = Array2::from_elem((n, l_after_filtering), f64::NAN);
        for i in 0..n {
            let mut k_new = 1; // Start at the second index after the intercept
            for j in 0..l_after_filtering {
                // Intercept
                new_intercept_and_allele_frequencies[(i, 0)] = 1.0;
                let j_ = idx[j];
                let idx_ini = loci_idx_ini[j_];
                let idx_fin = loci_idx_fin[j_];
                for k in 0..(idx_fin - idx_ini) {
                    let k_old = idx_ini + k;
                    new_intercept_and_allele_frequencies[(i, k_new)] =
                        self.intercept_and_allele_frequencies[(i, k_old)];
                    k_new += 1;
                }
                // Coverages
                new_coverages[(i, j)] = self.coverages[(i, j_)];
            }
        }
        self.chromosome = new_chromosome;
        self.position = new_position;
        self.allele = new_allele;
        self.intercept_and_allele_frequencies = new_intercept_and_allele_frequencies;
        self.coverages = new_coverages;
        match self.check() {
            Ok(x) => x,
            Err(e) => return Err(ImputefError{
                code: 208,
                message: "Error checking GenotypesAndPhenotypes in the method filter_out_top_missing_loci() | ".to_owned() +
                &e.message
            })
        };
        Ok(self)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_filtering_missing() {
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
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // NOTE: MAKE SURE THAT FilterStats is set to no filtering except by zero frequency alleles as in below:
        //       AND THAT WE KEEP ALL NON-ZERO ALLELES ACROSS ALL ENTRIES!
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        let filter_stats = FilterStats {
            remove_ns: false,
            min_quality: 1.0,
            min_coverage: 0,
            min_allele_frequency: 0.000001,
            max_missingness_rate: 0.0,
            pool_sizes: vec![20., 20., 20., 20., 20.],
        };
        let n_threads = 2;
        let mut frequencies_and_phenotypes = file_sync_phen
            .convert_into_genotypes_and_phenotypes(&filter_stats, false, &n_threads)
            .unwrap();
        let min_depth_below_which_are_missing = 5.0;
        let max_depth_above_which_are_missing = 100.0;
        frequencies_and_phenotypes
            .set_missing_by_depth(
                &min_depth_below_which_are_missing,
                &max_depth_above_which_are_missing,
            )
            .unwrap();
        // println!("frequencies_and_phenotypes={:?}", frequencies_and_phenotypes);
        println!(
            "frequencies_and_phenotypes.chromosome[0..5]={:?}",
            &frequencies_and_phenotypes.chromosome[0..5]
        );
        println!(
            "frequencies_and_phenotypes.position[0..5]={:?}",
            &frequencies_and_phenotypes.position[0..5]
        );
        println!(
            "frequencies_and_phenotypes.allele[0..5]={:?}",
            &frequencies_and_phenotypes.allele[0..5]
        );
        println!(
            "frequencies_and_phenotypes.intercept_and_allele_frequencies={:?}",
            frequencies_and_phenotypes.intercept_and_allele_frequencies
        );
        println!(
            "frequencies_and_phenotypes.coverages={:?}",
            frequencies_and_phenotypes.coverages
        );
        assert!(frequencies_and_phenotypes.intercept_and_allele_frequencies[(0, 1)].is_nan());
        assert_eq!(
            frequencies_and_phenotypes.intercept_and_allele_frequencies[(0, 4)],
            1.0 / 21.0
        );
        assert_eq!(
            frequencies_and_phenotypes.intercept_and_allele_frequencies[(3, 13)],
            2.0 / 10.0
        );
        assert!(frequencies_and_phenotypes.intercept_and_allele_frequencies[(3, 18)].is_nan());
        assert!(frequencies_and_phenotypes.intercept_and_allele_frequencies[(3, 18)].is_nan());

        println!("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
        println!("Before filtering out pools");
        println!(
            "frequencies_and_phenotypes.intercept_and_allele_frequencies:\n{:?}",
            frequencies_and_phenotypes.intercept_and_allele_frequencies
        );
        println!(
            "frequencies_and_phenotypes.coverages:\n{:?}",
            frequencies_and_phenotypes.coverages
        );
        frequencies_and_phenotypes
            .filter_out_top_missing_pools(&0.20)
            .unwrap();
        assert_eq!(
            frequencies_and_phenotypes
                .intercept_and_allele_frequencies
                .nrows(),
            4
        );
        println!("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
        println!("After filtering out pools but before filtering out loci");
        println!(
            "frequencies_and_phenotypes.intercept_and_allele_frequencies:\n{:?}",
            frequencies_and_phenotypes.intercept_and_allele_frequencies
        );
        println!(
            "frequencies_and_phenotypes.coverages:\n{:?}",
            frequencies_and_phenotypes.coverages
        );
        frequencies_and_phenotypes
            .filter_out_top_missing_loci(&1.00)
            .unwrap();
        assert_eq!(
            frequencies_and_phenotypes
                .intercept_and_allele_frequencies
                .ncols(),
            7637
        );
        println!("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
        println!("After filtering out pools and loci");
        println!(
            "frequencies_and_phenotypes.intercept_and_allele_frequencies:\n{:?}",
            frequencies_and_phenotypes.intercept_and_allele_frequencies
        );
        println!(
            "frequencies_and_phenotypes.coverages:\n{:?}",
            frequencies_and_phenotypes.coverages
        );
        // assert_eq!(0, 1);
    }
}
