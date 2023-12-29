use ndarray::prelude::*;
use std::fs::{File, OpenOptions};
use std::io::{self, prelude::*, BufReader, Error, ErrorKind, SeekFrom};
use std::str;
use std::sync::{Arc, Mutex};
use std::time::{SystemTime, UNIX_EPOCH};

use crate::helpers::*;
use crate::structs_and_traits::*;

impl CheckStruct for LocusCounts {
    fn check(&self) -> io::Result<()> {
        let (_n, p) = self.matrix.dim();
        let a = self.alleles_vector.len();
        match p == a {
            true => Ok(()),
            false => {
                let locus = [self.chromosome.clone(), self.position.to_string()].join("-");
                Err(Error::new(ErrorKind::Other, locus + " ::: LocusCounts error: the counts matrix does not have the same the number columns as the number of alleles."))
            }
        }
    }
}

impl CheckStruct for LocusFrequencies {
    fn check(&self) -> io::Result<()> {
        let (_n, p) = self.matrix.dim();
        let a = self.alleles_vector.len();
        match p == a {
            true => Ok(()),
            false => {
                let locus = [self.chromosome.clone(), self.position.to_string()].join("-");
                Err(Error::new(ErrorKind::Other, locus + " ::: LocusFrequencies error: the frequencies matrix does not have the same the number columns as the number of alleles."))
            }
        }
    }
}

impl CheckStruct for LocusCountsAndPhenotypes {
    fn check(&self) -> io::Result<()> {
        self.locus_counts.check().unwrap();
        let (n, _p) = self.locus_counts.matrix.dim();
        let (n_, _k) = self.phenotypes.dim();
        let n__ = self.pool_names.len();
        match (n == n_) & (n_ == n__) {
            true => Ok(()),
            false => {
                let locus = [
                    self.locus_counts.chromosome.clone(),
                    self.locus_counts.position.to_string(),
                ]
                .join("-");
                Err(Error::new(ErrorKind::Other, locus + " ::: LocusCountsAndPhenotypes error: the number of pools are inconsistent in the locus counts, and/or phenotypes matrix and/or pool names."))
            }
        }
    }
}

impl CheckStruct for GenotypesAndPhenotypes {
    fn check(&self) -> io::Result<()> {
        let p = self.chromosome.len();
        let p_ = self.position.len();
        let (n, p__) = self.intercept_and_allele_frequencies.dim();
        let (n_, _k) = self.phenotypes.dim();
        let n__ = self.pool_names.len();
        let (n___, l) = self.coverages.dim();
        match (p == p_) & (p_ == p__) & (n == n_) & (n_ == n__) & (n__ == n___) & (l <= p) {
            true => Ok(()),
            false => Err(Error::new(ErrorKind::Other, "GenotypesAndPhenotypes: there are at least 1 mismatch in the number of pools and loci."))
        }
    }
}

impl Count for GenotypesAndPhenotypes {
    fn count_loci(&self) -> io::Result<(Vec<usize>, Vec<String>, Vec<u64>)> {
        let (_, p) = self.intercept_and_allele_frequencies.dim();
        assert_eq!(p, self.chromosome.len(), "The number of entries in the 'chromosome' field and the total number of loci are incompatible. Please check the 'intercept_and_allele_frequencies' and 'chromosome' fields of 'GenotypesAndPhenotypes' struct.");
        assert_eq!(p, self.position.len(), "The number of entries in the 'position' field and the total number of loci are incompatible. Please check the 'intercept_and_allele_frequencies' and 'chromosome' fields of 'GenotypesAndPhenotypes' struct.");
        // Count the number of loci (Note: assumes the loci are sorted) and extract the loci coordinates
        let mut loci_idx: Vec<usize> = vec![];
        let mut loci_chr: Vec<String> = vec![];
        let mut loci_pos: Vec<u64> = vec![];
        for i in 1..p {
            // excludes the intercept
            if (self.chromosome[i - 1] != self.chromosome[i])
                | ((self.chromosome[i - 1] == self.chromosome[i])
                    & (self.position[i - 1] != self.position[i]))
            {
                loci_idx.push(i);
                loci_chr.push(self.chromosome[i].to_owned());
                loci_pos.push(self.position[i]);
            }
        }
        loci_idx.push(p); // last allele of the last locus
        loci_chr.push(self.chromosome.last().unwrap().to_owned()); // last allele of the last locus
        loci_pos.push(self.position.last().unwrap().to_owned()); // last allele of the last locus
        let l = loci_idx.len();
        assert_eq!(l-1, self.coverages.ncols(), "The number of loci with coverage information and the total number of loci are incompatible. You may have duplicate loci in the input genotype file (vcf, sync, or txt). Please check the 'intercept_and_allele_frequencies' and 'coverages' fields of 'GenotypesAndPhenotypes' struct.");
        Ok((loci_idx, loci_chr, loci_pos))
    }
}

impl Parse<LocusCounts> for String {
    // Parse a line of pileup into PileupLine struct
    fn lparse(&self) -> io::Result<Box<LocusCounts>> {
        // Remove trailing newline character in Unix-like (\n) and Windows (\r)
        let mut line = self.clone();
        if line.ends_with('\n') {
            line.pop();
            if line.ends_with('\r') {
                line.pop();
            }
        }
        // Ignore commented-out lines (i.e. '#' => 35)
        if line.as_bytes()[0] == 35_u8 {
            return Err(Error::new(ErrorKind::Other, "Commented out line"));
        }
        // Parse the sync line
        let vec_line = line
            .split('\t')
            .collect::<Vec<&str>>()
            .into_iter()
            .map(|x| x.to_owned())
            .collect::<Vec<String>>();
        let n: usize = vec_line.len() - 3;
        let p: usize = 6;
        let chromosome = vec_line[0].to_owned();
        let position = match vec_line[1].parse::<u64>() {
            Ok(x) => x,
            Err(_) => {
                return Err(Error::new(
                    ErrorKind::Other,
                    "Please check format of the file: position is not and integer.",
                ))
            }
        };
        let alleles_vector = vec!["A", "T", "C", "G", "N", "D"]
            .into_iter()
            .map(|x| x.to_owned())
            .collect::<Vec<String>>();
        // Read the allele counts into the counts matrix
        let mut matrix: Array2<u64> = Array2::from_elem((n, p), 0);
        let mut counts: Vec<u64>;
        for i in 0..n {
            counts = vec_line[i+3].split(':')
                                .map(|x| x.to_string().parse::<u64>().expect("Please check the input sync file as the allele counts are not valid integers."))
                                .collect::<Vec<u64>>();
            for j in 0..p {
                matrix[(i, j)] = counts[j];
            }
        }
        Ok(Box::new(LocusCounts {
            chromosome,
            position,
            alleles_vector,
            matrix, // n pools x 6 alleles
        }))
    }
}

impl Filter for LocusCounts {
    // PileupLine to AlleleCounts
    fn to_counts(&self) -> io::Result<Box<LocusCounts>> {
        let out = self.clone();
        Ok(Box::new(out))
    }

    // PileupLine to AlleleFrequencies
    fn to_frequencies(&self) -> io::Result<Box<LocusFrequencies>> {
        let n = self.matrix.nrows();
        let p = self.matrix.ncols();
        // let row_sums = self.matrix.sum_axis(Axis(1)); // summation across the columns which means sum of all elements per row
        let row_sums = self.matrix.map_axis(Axis(1), |r| {
            r.map(|&x| x as f64)
                .iter()
                .filter(|&&x| !x.is_nan())
                .fold(0.0, |sum, &x| sum + x)
        }); // summation across the columns which means sum of all elements per row while ignoring NANs!
        let mut matrix: Array2<f64> = Array2::from_elem((n, p), 0.0_f64);
        for i in 0..n {
            for j in 0..p {
                matrix[(i, j)] = if row_sums[i] == 0.0 {
                    f64::NAN
                } else {
                    self.matrix[(i, j)] as f64 / row_sums[i]
                };
            }
        }
        Ok(Box::new(LocusFrequencies {
            chromosome: self.chromosome.clone(),
            position: self.position,
            alleles_vector: self.alleles_vector.clone(),
            matrix,
        }))
    }

    // Filter PileupLine by minimum coverage, minimum quality
    fn filter(&mut self, filter_stats: &FilterStats) -> io::Result<&mut Self> {
        // Cannot filter by base qualities as this information is lost and we are assuming this has been performed during pileup to sync conversion
        // Preliminary check of the structure format
        self.check().unwrap();
        // Remove Ns
        if filter_stats.remove_ns {
            let i = match self
                .alleles_vector
                .iter()
                .position(|x| (x == &"N".to_owned()) | (x == &"n".to_owned()))
            {
                Some(x) => x as i32,
                None => -1,
            };
            if i != -1 {
                self.alleles_vector.remove(i as usize);
                self.matrix.remove_index(Axis(1), i as usize);
            }
        }
        // println!("self={:?}", self);
        // Filter by minimum coverage
        // Summation across the columns which means sum of all elements per row while ignoring NANs!
        let sum_coverage = self.matrix.map_axis(Axis(1), |r| {
            r.map(|&x| x as f64)
                .iter()
                .filter(|&&x| !x.is_nan())
                .fold(0.0, |sum, &x| sum + x)
        });
        let min_sum_coverage =
            sum_coverage
                .iter()
                .fold(sum_coverage[0], |min, &x| if x < min { x } else { min });
        if min_sum_coverage < filter_stats.min_coverage as f64 {
            return Err(Error::new(ErrorKind::Other, "Filtered out."));
        }
        // // TODO: convert loci failing the minimum coverage threshold into missing instead of omitting the entire locus
        // for i in 0..self.matrix.nrows() {
        //     if sum_coverage[i] < filter_stats.min_coverage as f64 {
        //         for j in 0..self.matrix.ncols() {
        //             self.matrix[(i, j)] = f64::NAN as u64;
        //         }
        //     }
        // }
        // Filter by minimum allele frequency
        // Before anything else, we clone matrix of allele counts
        let mut matrix = self.matrix.clone();
        //// First convert allele counts into frequencies
        let mut allele_frequencies = match self.to_frequencies() {
            Ok(x) => x,
            Err(_) => {
                return Err(Error::new(
                    ErrorKind::Other,
                    "Cannot convert locus counts to locus frequencies.",
                ))
            }
        };
        //// Next account for pool sizes to get the proper minmum allele frequency across all pools
        let n = allele_frequencies.matrix.nrows();
        let mut p = allele_frequencies.matrix.ncols();
        assert!(
            n == filter_stats.pool_sizes.len(),
            "Please make that the number of pools and the pool sizes in FilterStats match."
        );
        let mut q: f64;
        let mut j: usize = 0;
        while j < p {
            q = 0.0;
            for i in 0..n {
                q += match allele_frequencies.matrix[(i, j)].is_nan() {
                    true => 0.0,
                    false => {
                        allele_frequencies.matrix[(i, j)]
                            * (filter_stats.pool_sizes[i]
                                / filter_stats.pool_sizes.iter().sum::<f64>())
                    }
                };
            }
            if (q < filter_stats.min_allele_frequency)
                | (q > (1.00 - filter_stats.min_allele_frequency))
            {
                allele_frequencies.matrix.remove_index(Axis(1), j);
                matrix.remove_index(Axis(1), j);
                self.alleles_vector.remove(j);
                p -= 1;
            } else {
                j += 1;
            }
        }
        // Check if all alleles have failed the minimum allele frequency, i.e. the locus has been filtered out
        if p < 2 {
            return Err(Error::new(ErrorKind::Other, "Filtered out."));
        }
        // Filter out if the locus is missing across all pools using the first allele where if the locus is missing then all
        let (n, _p) = allele_frequencies.matrix.dim();
        let n_missing_across_pools = allele_frequencies
            .matrix
            .slice(s![.., 0])
            .fold(0, |sum, &x| if x.is_nan() { sum + 1 } else { sum });
        if n_missing_across_pools == n {
            return Err(Error::new(ErrorKind::Other, "Filtered out."));
        }
        // Filter-out the locus if the rate of missingness, i.e. the fraction of the pools missing coverage of the current locus is below the minimum threshold
        if (n_missing_across_pools as f64 / n as f64) > filter_stats.max_missingness_rate {
            return Err(Error::new(ErrorKind::Other, "Filtered out."));
        }
        // Return the locus if it passed all the filtering steps
        self.matrix = matrix;
        Ok(self)
    }
}

impl Filter for LocusFrequencies {
    // PileupLine to AlleleCounts
    fn to_counts(&self) -> io::Result<Box<LocusCounts>> {
        let n = self.matrix.nrows();
        let p = self.matrix.ncols();
        let mut matrix: Array2<u64> = Array2::from_elem((n, p), 0);
        let mut max_n: f64;
        for i in 0..n {
            let row = self.matrix.row(i);
            let min = row
                .iter()
                .filter(|&&x| x != 0.0)
                .fold(1.0, |min, &x| if x < min { x } else { min });
            if min == 0.0 {
                return Err(Error::new(
                    ErrorKind::Other,
                    "At least one of the pools have no coverage.",
                ));
            }
            max_n = 1.00 / min;
            for j in 0..p {
                matrix[(i, j)] = (max_n * self.matrix[(i, j)]).round() as u64;
            }
        }
        Ok(Box::new(LocusCounts {
            chromosome: self.chromosome.clone(),
            position: self.position,
            alleles_vector: self.alleles_vector.clone(),
            matrix,
        }))
    }

    // PileupLine to AlleleFrequencies
    fn to_frequencies(&self) -> io::Result<Box<LocusFrequencies>> {
        // Recompute the frequencies using frequencies when the number of colulmns or one or more alleles have been filtered out/removed
        let n = self.matrix.nrows();
        let p = self.matrix.ncols();
        // let row_sums = self.matrix.sum_axis(Axis(1)); // summation across the columns which means sum of all elements per row
        let row_sums = self.matrix.map_axis(Axis(1), |r| {
            r.iter()
                .filter(|&&x| !x.is_nan())
                .fold(0.0, |sum, &x| sum + x)
        }); // summation across the columns which means sum of all elements per row while ignoring NANs!
            // // Make sure all pools have been convered
            // if row_sums
            //     .iter()
            //     .fold(row_sums[0], |min, &x| if x < min { x } else { min })
            //     == 0.0
            // {
            //     return Err(Error::new(
            //         ErrorKind::Other,
            //         "At least one pool did not have coverage.",
            //     ));
            // }
        let mut matrix: Array2<f64> = Array2::from_elem((n, p), 0.0_f64);
        for i in 0..n {
            for j in 0..p {
                matrix[(i, j)] = if row_sums[i] == 0.0 {
                    f64::NAN
                } else {
                    self.matrix[(i, j)] / row_sums[i]
                };
            }
        }
        Ok(Box::new(LocusFrequencies {
            chromosome: self.chromosome.clone(),
            position: self.position,
            alleles_vector: self.alleles_vector.clone(),
            matrix,
        }))
    }

    // Filter PileupLine by minimum coverage, minimum quality
    fn filter(&mut self, filter_stats: &FilterStats) -> io::Result<&mut Self> {
        // Cannot filter by base qualities as this information is lost and we are assuming this has been performed during pileup to sync conversion
        // Also, cannot filter by minimum coverage as that data is lost from counts to frequencies conversion
        // Preliminary check of the structure format
        self.check().unwrap();
        // Remove Ns
        if filter_stats.remove_ns {
            let i = match self
                .alleles_vector
                .iter()
                .position(|x| (x == &"N".to_owned()) | (x == &"n".to_owned()))
            {
                Some(x) => x as i32,
                None => -1,
            };
            if i != -1 {
                self.alleles_vector.remove(i as usize);
                self.matrix.remove_index(Axis(1), i as usize);
            }
        }
        // println!("self={:?}", self);
        // Recompute frequencies after removing Ns
        let recomputed_self = match self.to_frequencies() {
            Ok(x) => x,
            Err(_) => {
                return Err(Error::new(
                    ErrorKind::Other,
                    "T_T Cannot convert locus counts to locus frequencies.",
                ))
            }
        };
        self.alleles_vector = recomputed_self.alleles_vector;
        self.matrix = recomputed_self.matrix;
        // Filter by minimum allele frequency
        // Before anything else, we clone matrix of allele frequencies
        let mut matrix = self.matrix.clone();
        //// First convert allele counts into frequencies
        let mut allele_frequencies = match self.to_frequencies() {
            Ok(x) => x,
            Err(_) => {
                return Err(Error::new(
                    ErrorKind::Other,
                    "Cannot convert locus counts to locus frequencies.",
                ))
            }
        };
        //// Next account for pool sizes to get the proper minmum allele frequency across all pools
        let n = allele_frequencies.matrix.nrows();
        let p_ = allele_frequencies.matrix.ncols();
        let mut p = p_;
        let mut q: f64;
        let mut j: usize = 0;
        while j < p {
            q = 0.0;
            for i in 0..n {
                q += match allele_frequencies.matrix[(i, j)].is_nan() {
                    true => 0.0,
                    false => {
                        allele_frequencies.matrix[(i, j)]
                            * (filter_stats.pool_sizes[i]
                                / filter_stats.pool_sizes.iter().sum::<f64>())
                    }
                };
            }
            if (q < filter_stats.min_allele_frequency)
                | (q > (1.00 - filter_stats.min_allele_frequency))
            {
                allele_frequencies.matrix.remove_index(Axis(1), j);
                matrix.remove_index(Axis(1), j);
                self.alleles_vector.remove(j);
                p -= 1;
            } else {
                j += 1;
            }
        }
        // Check if all alleles have failed the minimum allele frequency, i.e. the locus has been filtered out
        if p < 2 {
            return Err(Error::new(ErrorKind::Other, "Filtered out."));
        }
        // Filter out if the locus is missing across all pools using the first allele where if the locus is missing then all
        let (n, _p) = allele_frequencies.matrix.dim();
        let n_missing_across_pools = allele_frequencies
            .matrix
            .slice(s![.., 0])
            .fold(0, |sum, &x| if x.is_nan() { sum + 1 } else { sum });
        if n_missing_across_pools == n {
            return Err(Error::new(ErrorKind::Other, "Filtered out."));
        }
        // Filter-out the locus if the rate of missingness, i.e. the fraction of the pools missing coverage of the current locus is below the minimum threshold
        if (n_missing_across_pools as f64 / n as f64) > filter_stats.max_missingness_rate {
            return Err(Error::new(ErrorKind::Other, "Filtered out."));
        }
        // Correct allele frequencies if one or more alleles were filtered out
        for i in 0..matrix.nrows() {
            let row_sum = matrix.row(i).sum();
            for j in 0..matrix.ncols() {
                matrix[(i, j)] = matrix[(i, j)] / row_sum;
            }
        }
        // Return the locus if it passed all the filtering steps
        self.matrix = matrix;
        Ok(self)
    }
}

impl Sort for LocusFrequencies {
    fn sort_by_allele_freq(&mut self, decreasing: bool) -> io::Result<&mut Self> {
        let n = self.matrix.nrows();
        let p = self.matrix.ncols();
        let mut sorted_matrix: Array2<f64> = Array2::from_elem((n, p), f64::NAN);
        let mut sorted_alleles_vector: Vec<String> = vec![];
        let mut idx = (0..self.matrix.ncols()).collect::<Vec<usize>>();
        // let column_sums = self.matrix.sum_axis(Axis(0));
        let column_sums = self.matrix.map_axis(Axis(0), |r| {
            r.iter()
                .filter(|&&x| !x.is_nan())
                .fold(0.0, |sum, &x| sum + x)
        }); // ignoring NANs!
            // println!("self={:?}", self);
        if decreasing {
            idx.sort_by(|&a, &b| column_sums[b].partial_cmp(&column_sums[a]).unwrap());
        } else {
            idx.sort_by(|&a, &b| column_sums[a].partial_cmp(&column_sums[b]).unwrap());
        }
        for i in 0..p {
            sorted_matrix
                .column_mut(i)
                .assign(&self.matrix.column(idx[i]));
            sorted_alleles_vector.push(self.alleles_vector[idx[i]].clone());
        }
        self.matrix = sorted_matrix;
        self.alleles_vector = sorted_alleles_vector;
        Ok(self)
    }
}

impl RemoveMissing for LocusCountsAndPhenotypes {
    // Remove pools with missing data in the phenotype file
    fn remove_missing(&mut self) -> io::Result<&mut Self> {
        let (n, k) = self.phenotypes.dim();
        let n_ = self.pool_names.len();
        let (n__, p) = self.locus_counts.matrix.dim();
        assert_eq!(n, n_);
        assert_eq!(n, n__);
        let pool_means: Array1<f64> = self.phenotypes.mean_axis(Axis(1)).unwrap();
        let mut idx: Vec<usize> = vec![];
        for i in 0..n {
            if !pool_means[i].is_nan() {
                idx.push(i);
            }
        }
        if !idx.is_empty() {
            let mut new_phenotypes: Array2<f64> = Array2::from_elem((idx.len(), k), f64::NAN);
            let mut new_pool_names: Vec<String> = vec![];
            let mut new_locus_counts_matrix: Array2<u64> = Array2::from_elem((idx.len(), p), 0);
            let mut i_new: usize = 0;
            for i in idx {
                for j in 0..k {
                    new_phenotypes[(i_new, j)] = self.phenotypes[(i, j)];
                }
                new_pool_names.push(self.pool_names[i].clone());
                for j in 0..p {
                    new_locus_counts_matrix[(i_new, j)] = self.locus_counts.matrix[(i, j)];
                }
                i_new += 1;
            }
            self.phenotypes = new_phenotypes;
            self.pool_names = new_pool_names;
            self.locus_counts.matrix = new_locus_counts_matrix;
        } else {
            return Err(Error::new(
                ErrorKind::Other,
                "All pools have missing data. Please check the phenotype file.",
            ));
        }
        Ok(self)
    }
}

impl RemoveMissing for GenotypesAndPhenotypes {
    // Remove pools with missing data in the phenotype file
    fn remove_missing(&mut self) -> io::Result<&mut Self> {
        let (n, k) = self.phenotypes.dim();
        let n_ = self.pool_names.len();
        let (n__, p) = self.intercept_and_allele_frequencies.dim();
        let (n___, l) = self.coverages.dim();
        assert_eq!(n, n_);
        assert_eq!(n, n__);
        assert_eq!(n, n___);
        let pool_means: Array1<f64> = self.phenotypes.mean_axis(Axis(1)).unwrap();
        let mut idx: Vec<usize> = vec![];
        for i in 0..n {
            if !pool_means[i].is_nan() {
                idx.push(i);
            }
        }
        if !idx.is_empty() {
            let mut new_phenotypes: Array2<f64> = Array2::from_elem((idx.len(), k), f64::NAN);
            let mut new_pool_names: Vec<String> = vec![];
            let mut new_intercept_and_allele_frequencies: Array2<f64> =
                Array2::from_elem((idx.len(), p), f64::NAN);
            let mut new_coverages: Array2<f64> = Array2::from_elem((idx.len(), l), f64::NAN);
            let mut i_new: usize = 0;
            for i in idx.clone() {
                for j in 0..k {
                    new_phenotypes[(i_new, j)] = self.phenotypes[(i, j)];
                }
                new_pool_names.push(self.pool_names[i].clone());
                for j in 0..p {
                    new_intercept_and_allele_frequencies[(i_new, j)] =
                        self.intercept_and_allele_frequencies[(i, j)];
                }
                for j in 0..l {
                    new_coverages[(i_new, j)] = self.coverages[(i, j)];
                }
                i_new += 1;
            }
            self.phenotypes = new_phenotypes;
            self.pool_names = new_pool_names;
            self.intercept_and_allele_frequencies = new_intercept_and_allele_frequencies;
            self.coverages = new_coverages;
        } else {
            return Err(Error::new(
                ErrorKind::Other,
                "All pools have missing data. Please check the phenotype file.",
            ));
        }
        // println!("self={:?}", self);
        // println!("self.phenotypes={:?}", self.phenotypes);
        // println!("idx={:?}", idx);
        Ok(self)
    }
}

impl LoadAll for FileSyncPhen {
    fn per_chunk_load(
        &self,
        start: &u64,
        end: &u64,
        filter_stats: &FilterStats,
        keep_p_minus_1: bool,
    ) -> io::Result<(Vec<LocusFrequencies>, Vec<LocusCounts>)> {
        // Input syn file
        let fname = self.filename_sync.clone();

        // Prepare output vectors
        let mut freq: Vec<LocusFrequencies> = Vec::new();
        let mut cnts: Vec<LocusCounts> = Vec::new();
        // Input file chunk
        let file = File::open(fname.clone()).unwrap();
        let mut reader = BufReader::new(file);
        // Navigate to the start of the chunk
        let mut i: u64 = *start;
        reader.seek(SeekFrom::Start(*start)).unwrap();
        // Read and parse until the end of the chunk
        while i < *end {
            // Instantiate the line
            let mut line = String::new();
            // Read the line which automatically movesthe cursor position to the next line
            let _ = reader.read_line(&mut line).unwrap();
            // Find the new cursor position
            i = reader.stream_position().unwrap();
            // Remove trailing newline character in Unix-like (\n) and Windows (\r)
            if line.ends_with('\n') {
                line.pop();
                if line.ends_with('\r') {
                    line.pop();
                }
            }
            // Parse the pileup line
            let mut locus_counts: LocusCounts = match line.lparse() {
                Ok(x) => *x,
                Err(x) => match x.kind() {
                    ErrorKind::Other => continue,
                    _ => {
                        return Err(Error::new(
                            ErrorKind::Other,
                            "T_T Input sync file error, i.e. '".to_owned()
                                + &fname
                                + "' at line with the first 20 characters as: "
                                + &line[0..20]
                                + ".",
                        ))
                    }
                },
            };
            match locus_counts.filter(filter_stats) {
                Ok(x) => x,
                Err(_) => continue,
            };
            let mut locus_frequencies = match locus_counts.to_frequencies() {
                Ok(x) => *x,
                Err(_) => continue,
            };
            // Remove minor allele
            if keep_p_minus_1 {
                locus_frequencies.sort_by_allele_freq(true).unwrap();
                locus_frequencies.matrix.remove_index(Axis(1), 0);
                locus_frequencies.alleles_vector.remove(0);
            }
            freq.push(locus_frequencies);
            cnts.push(locus_counts);
        }
        Ok((freq, cnts))
    }

    fn load(
        &self,
        filter_stats: &FilterStats,
        keep_p_minus_1: bool,
        n_threads: &usize,
    ) -> io::Result<(Vec<LocusFrequencies>, Vec<LocusCounts>)> {
        let fname = self.filename_sync.clone();
        // Find the positions whereto split the file into n_threads pieces
        let chunks = find_file_splits(&fname, n_threads).unwrap();
        let n_threads = chunks.len() - 1;
        println!("Chunks: {:?}", chunks);
        // Tuple arguments of pileup2sync_chunks
        // Instantiate thread object for parallel execution
        let mut thread_objects = Vec::new();
        // Vector holding all returns from pileup2sync_chunk()
        let thread_ouputs_freq: Arc<Mutex<Vec<LocusFrequencies>>> =
            Arc::new(Mutex::new(Vec::new())); // Mutated within each thread worker
        // Mutated within each thread worker
        let thread_ouputs_cnts: Arc<Mutex<Vec<LocusCounts>>> = Arc::new(Mutex::new(Vec::new()));
        for i in 0..n_threads {
            // Clone pileup2sync_chunk parameters
            let self_clone = self.clone();
            let start = chunks[i];
            let end = chunks[i + 1];
            let filter_stats = filter_stats.clone();
            let thread_ouputs_freq_clone = thread_ouputs_freq.clone(); // Mutated within the current thread worker
            let thread_ouputs_cnts_clone = thread_ouputs_cnts.clone(); // Mutated within the current thread worker
            let thread = std::thread::spawn(move || {
                let (mut freq, mut cnts) = self_clone
                    .per_chunk_load(&start, &end, &filter_stats, keep_p_minus_1)
                    .unwrap();
                thread_ouputs_freq_clone.lock().unwrap().append(&mut freq);
                thread_ouputs_cnts_clone.lock().unwrap().append(&mut cnts);
            });
            thread_objects.push(thread);
        }
        // Waiting for all threads to finish
        for thread in thread_objects {
            thread.join().expect("Unknown thread error occured.");
        }
        // Extract output filenames from each thread into a vector and sort them
        let mut freq: Vec<LocusFrequencies> = Vec::new();
        let mut cnts: Vec<LocusCounts> = Vec::new();
        for x in thread_ouputs_freq.lock().unwrap().iter() {
            freq.push(x.clone());
        }
        for x in thread_ouputs_cnts.lock().unwrap().iter() {
            cnts.push(x.clone());
        }
        freq.sort_by(|a, b| {
            a.chromosome
                .cmp(&b.chromosome)
                .then(a.position.cmp(&b.position))
        });
        cnts.sort_by(|a, b| {
            a.chromosome
                .cmp(&b.chromosome)
                .then(a.position.cmp(&b.position))
        });

        Ok((freq, cnts))
    }

    fn into_genotypes_and_phenotypes(
        &self,
        filter_stats: &FilterStats,
        keep_p_minus_1: bool,
        n_threads: &usize,
    ) -> io::Result<GenotypesAndPhenotypes> {
        let (freqs, cnts) = self.load(filter_stats, keep_p_minus_1, n_threads).unwrap();
        let n = self.pool_names.len();
        let m = freqs.len(); // total number of loci
                             // Find the total number of alleles across all loci
        let mut p = 1; // start with the intercept
        for f in freqs.iter() {
            p += f.matrix.ncols();
        }
        // println!("p={}", p);
        let mut chromosome: Vec<String> = Vec::with_capacity(p);
        chromosome.push("intercept".to_owned());
        let mut position: Vec<u64> = Vec::with_capacity(p);
        position.push(0);
        let mut allele: Vec<String> = Vec::with_capacity(p);
        allele.push("intercept".to_owned());
        let mut coverages: Array2<f64> = Array2::from_elem((n, m), f64::NAN);
        let mut l: usize = 0; // locus index
        let mut mat: Array2<f64> = Array2::from_elem((n, p), 1.0);
        let mut j: usize = 1; // SNP index across loci, start after the intercept
        assert_eq!(
            freqs.len(),
            cnts.len(),
            "Frequencies and counts not the same length."
        );
        for idx in 0..freqs.len() {
            // Allele frequencies
            let f = &freqs[idx];
            for j_ in 0..f.matrix.ncols() {
                chromosome.push(f.chromosome.clone());
                position.push(f.position);
                allele.push(f.alleles_vector[j_].clone());
                for i in 0..f.matrix.nrows() {
                    mat[(i, j)] = f.matrix[(i, j_)];
                }
                j += 1; // next allele
            }
            // Coverages
            let c = &cnts[idx];
            let cov: Array1<f64> = c.matrix.map_axis(Axis(1), |r| {
                r.map(|&x| x as f64)
                    .iter()
                    .filter(|&&x| !x.is_nan())
                    .fold(0.0, |sum, &x| sum + x)
            }); // summation across the columns which means sum of all elements per row while ignoring NANs!
            for i in 0..cov.len() {
                coverages[(i, l)] = cov[i];
            }
            l += 1; // next locus
        }
        // println!("coverages={:?}", coverages);
        // println!("mat={:?}", mat);
        // println!("mat={:?}", mat.slice(s![0..5, 0..4]));
        // println!("chromosome[0]={:?}", chromosome[0]);
        // println!("chromosome[1]={:?}", chromosome[1]);
        // println!("chromosome[2]={:?}", chromosome[2]);
        // println!("chromosome[3]={:?}", chromosome[3]);
        // println!("position[0]={:?}", position[0]);
        // println!("position[1]={:?}", position[1]);
        // println!("position[2]={:?}", position[2]);
        // println!("position[3]={:?}", position[3]);
        Ok(GenotypesAndPhenotypes {
            chromosome,
            position,
            allele,
            intercept_and_allele_frequencies: mat,
            phenotypes: self.phen_matrix.clone(),
            pool_names: self.pool_names.clone(),
            coverages,
        })
    }
}

impl SaveCsv for FileSyncPhen {
    fn write_csv(
        &self,
        filter_stats: &FilterStats,
        keep_p_minus_1: bool,
        out: &String,
        n_threads: &usize,
    ) -> io::Result<String> {
        // Output filename
        let out = if *out == "".to_owned() {
            let time = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_secs_f64();
            let bname = self
                .filename_sync
                .split('.')
                .collect::<Vec<&str>>()
                .into_iter()
                .map(|a| a.to_owned())
                .collect::<Vec<String>>()
                .into_iter()
                .rev()
                .collect::<Vec<String>>()[1..]
                .to_owned()
                .into_iter()
                .rev()
                .collect::<Vec<String>>()
                .join(".");
            bname.to_owned() + "-" + &time.to_string() + "-allele_frequencies.csv"
        } else {
            out.clone()
        };
        // Instantiate output file
        let error_writing_file = "Unable to create file: ".to_owned() + &out;
        let mut file_out = OpenOptions::new()
            .create_new(true)
            .write(true)
            .append(false)
            .open(&out)
            .expect(&error_writing_file);
        // Load the full sync file in parallel and sort
        let (freqs, _cnts) = self.load(filter_stats, keep_p_minus_1, n_threads).unwrap();
        // Make sure that we have the same number of pools in the genotype and phenotype files
        assert!(!freqs.is_empty(), "No data passed the filtering variables. Please decrease minimum depth, and/or minimum allele frequency.");
        assert!(
            freqs[0].matrix.nrows() == self.pool_names.len(),
            "Please check that the pools are consistent across the genotype and phenotype files."
        );
        // Write the header
        file_out
            .write_all(
                ("#chr,pos,allele,".to_owned() + &self.pool_names.join(",") + "\n").as_bytes(),
            )
            .unwrap();
        // Write allele frequencies line by line
        for f in freqs.iter() {
            for i in 0..f.alleles_vector.len() {
                let freqs_per_pool = f
                    .matrix
                    .column(i)
                    .iter()
                    .map(|x| parse_f64_roundup_and_own(*x, 6))
                    .collect::<Vec<String>>()
                    .join(",");
                let line = [
                    f.chromosome.to_owned(),
                    f.position.to_string(),
                    f.alleles_vector[i].to_owned(),
                    freqs_per_pool,
                ]
                .join(",")
                    + "\n";
                file_out.write_all(line.as_bytes()).unwrap();
            }
        }
        Ok(out)
    }
}

impl SaveCsv for GenotypesAndPhenotypes {
    fn write_csv(
        &self,
        _filter_stats: &FilterStats,
        _keep_p_minus_1: bool,
        out: &String,
        _n_threads: &usize,
    ) -> io::Result<String> {
        // Note: All input parameters are not used except for one - out, the rest are for other implementations of this trait i.e. filter_stats, keep_p_minus_1, and n_threads
        // Sanity checks
        let (n, p) = self.intercept_and_allele_frequencies.dim();
        let (n_, _l) = self.coverages.dim();
        let (n__, _k) = self.phenotypes.dim();
        let n___ = self.pool_names.len();
        let p_ = self.chromosome.len();
        let p__ = self.position.len();
        let p___ = self.allele.len();
        assert_eq!(p_, p__, "Please check the genotypes and phenotypes data: the number of elements in the chromosome names vector is not equal to the number of elements in the positions vector.");
        assert_eq!(p_, p___, "Please check the genotypes and phenotypes data: the number of elements in the chromosome names vector is not equal to the number of elements in the alleles vector.");
        assert_eq!(p_, p, "Please check the genotypes and phenotypes data: the number of elements in the chromosome names vector is not equal to the number of columns in the allele frequencies matrix.");
        assert_eq!(n___, n__, "Please check the genotypes and phenotypes data: the number of elements in the pool names vector is not equal to the number of rows in the phenotypes matrix.");
        assert_eq!(n___, n_, "Please check the genotypes and phenotypes data: the number of elements in the pool names vector is not equal to the number of rows in the coverages matrix.");
        assert_eq!(n___, n, "Please check the genotypes and phenotypes data: the number of elements in the pool names vector is not equal to the number of rows in the allele frequencies matrix.");
        // Output filename
        let out = if *out == "".to_owned() {
            let time = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_secs_f64();
            "genotypes_and_phenotypes".to_owned()
                + "-"
                + &time.to_string()
                + "-allele_frequencies.csv"
        } else {
            out.clone()
        };
        // Instantiate output file
        let error_writing_file = "Unable to create file: ".to_owned() + &out;
        let mut file_out = OpenOptions::new()
            .create_new(true)
            .write(true)
            .append(false)
            .open(&out)
            .expect(&error_writing_file);
        // Write the header
        file_out
            .write_all(
                ("#chr,pos,allele,".to_owned() + &self.pool_names.join(",") + "\n").as_bytes(),
            )
            .unwrap();
        // Write allele frequencies line by line (skip the intercept)
        for i in 1..p {
            let freqs_per_pool = self
                .intercept_and_allele_frequencies
                .column(i)
                .iter()
                .map(|&x| parse_f64_roundup_and_own(x, 6))
                .collect::<Vec<String>>()
                .join(",");
            let line = [
                self.chromosome[i].to_owned(),
                self.position[i].to_string(),
                self.allele[i].to_owned(),
                freqs_per_pool,
            ]
            .join(",")
                + "\n";
            file_out.write_all(line.as_bytes()).unwrap();
        }
        Ok(out)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use rand::prelude::Distribution;

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    #[test]
    fn test_sync_check() {
        // (01) CheckStruct: LocusCounts
        let mut locus_counts = LocusCounts {
            chromosome: "chr1".to_owned(),
            position: 12345,
            alleles_vector: vec!["A".to_owned()], // Missing an allele in alleles_vector
            matrix: Array2::from_shape_vec((5, 2), vec![10, 90, 20, 80, 50, 50, 80, 20, 90, 10])
                .unwrap(),
        };
        // Missing an allele in alleles_vector
        assert!(locus_counts.check().is_err());
        // Append the missing allele
        locus_counts.alleles_vector.push("D".to_owned());
        // No error
        assert_eq!((), locus_counts.check().unwrap());
        // (02) CheckStruct: LocusFrequencies
        let mut locus_frequncies = LocusFrequencies {
            chromosome: "chr1".to_owned(),
            position: 12345,
            alleles_vector: vec!["A".to_owned()], // Missing an allele in alleles_vector
            matrix: Array2::from_shape_vec(
                (5, 2),
                vec![0.1, 0.9, 0.2, 0.8, 0.5, 0.5, 0.8, 0.2, 0.9, 0.1],
            )
            .unwrap(),
        };
        // Missing an allele in alleles_vector
        assert!(locus_frequncies.check().is_err());
        // Append the missing allele
        locus_frequncies.alleles_vector.push("D".to_owned());
        // No error
        assert_eq!((), locus_frequncies.check().unwrap());
        // (03) CheckStruct: LocusCountsAndPhenotypes
        let mut locus_counts_and_phenotypes = LocusCountsAndPhenotypes {
            locus_counts,
            phenotypes: Array2::from_shape_vec((4, 1), vec![0.1, 0.2, 0.3, 0.4]).unwrap(),
            pool_names: (1..5)
                .map(|x| "Pop".to_owned() + &x.to_string()[..])
                .collect(),
        };
        // Missing phenotype, and pool name
        assert!(locus_counts_and_phenotypes.check().is_err());
        // Replace phenotype array with the complete data
        locus_counts_and_phenotypes.phenotypes =
            Array2::from_shape_vec((5, 1), vec![0.1, 0.2, 0.3, 0.4, 0.5]).unwrap();
        // Missing pool name
        assert!(locus_counts_and_phenotypes.check().is_err());
        // Append missing pool name
        locus_counts_and_phenotypes
            .pool_names
            .push("Pop5".to_owned());
        // No error
        assert_eq!((), locus_counts_and_phenotypes.check().unwrap());
        // (04) CheckStruct: GenotypesAndPhenotypes
        let rng = rand::thread_rng();
        let dist_gaus = statrs::distribution::Normal::new(0.0, 1.0).unwrap();
        let mut genotypes_and_phenotypes = GenotypesAndPhenotypes {
            chromosome: [
                "intercept",
                "chr1",
                "chr1",
                "chr1",
                "chr2",
                "chr2",
                "chr2",
                "chrX",
                "chrX",
            ]
            .iter()
            .map(|&x| x.to_owned())
            .collect(),
            position: vec![0, 123, 123, 123, 50005, 50005, 50005, 701, 701],
            allele: ["intercept", "A", "T", "C", "A", "T", "C", "A", "T"]
                .iter()
                .map(|&x| x.to_owned())
                .collect(),
            intercept_and_allele_frequencies: Array2::from_shape_vec(
                (5, 9),
                vec![
                    1.0, 0.1, 0.0, 0.9, 0.1, 0.1, 0.8, 0.6, 0.4, 1.0, 0.0, 0.1, 0.9, 0.2, 0.0, 0.8,
                    0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 0.1, 0.9, 0.0, 0.1, 0.9, 1.0, 0.0, 0.0, 1.0, 0.0,
                    0.0, 1.0, 0.3, 0.7, 1.0, 0.3, 0.5, 0.2, 0.0, 0.5, 0.5, 0.5, 0.5,
                ],
            )
            .unwrap(),
            phenotypes: Array2::from_shape_vec(
                (5, 2),
                dist_gaus.sample_iter(rng.clone()).take(5 * 2).collect(),
            )
            .unwrap(),
            pool_names: (0..4)
                .map(|x| "Pop".to_owned() + &x.to_string()[..])
                .collect(), // missing the 5th pool
            coverages: Array2::from_shape_vec(
                (5, 3),
                std::iter::repeat(100.0).take(5 * 3).collect(),
            )
            .unwrap(),
        };
        println!("genotypes_and_phenotypes={:?}", genotypes_and_phenotypes);
        assert!(genotypes_and_phenotypes.check().is_err());
        // Add missing poool
        genotypes_and_phenotypes.pool_names.push("Pop5".to_owned());
        // No error
        assert_eq!((), genotypes_and_phenotypes.check().unwrap());
    }

    #[test]
    fn test_sync_count() {
        let rng = rand::thread_rng();
        let dist_gaus = statrs::distribution::Normal::new(0.0, 1.0).unwrap();
        let genotypes_and_phenotypes = GenotypesAndPhenotypes {
            chromosome: [
                "intercept",
                "chr1",
                "chr1",
                "chr1",
                "chr2",
                "chr2",
                "chr2",
                "chrX",
                "chrX",
            ]
            .iter()
            .map(|&x| x.to_owned())
            .collect(),
            position: vec![0, 123, 123, 123, 50005, 50005, 50005, 701, 701],
            allele: ["intercept", "A", "T", "C", "A", "T", "C", "A", "T"]
                .iter()
                .map(|&x| x.to_owned())
                .collect(),
            intercept_and_allele_frequencies: Array2::from_shape_vec(
                (5, 9),
                vec![
                    1.0, 0.1, 0.0, 0.9, 0.1, 0.1, 0.8, 0.6, 0.4, 1.0, 0.0, 0.1, 0.9, 0.2, 0.0, 0.8,
                    0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 0.1, 0.9, 0.0, 0.1, 0.9, 1.0, 0.0, 0.0, 1.0, 0.0,
                    0.0, 1.0, 0.3, 0.7, 1.0, 0.3, 0.5, 0.2, 0.0, 0.5, 0.5, 0.5, 0.5,
                ],
            )
            .unwrap(),
            phenotypes: Array2::from_shape_vec(
                (5, 2),
                dist_gaus.sample_iter(rng.clone()).take(5 * 2).collect(),
            )
            .unwrap(),
            pool_names: (0..4)
                .map(|x| "Pop".to_owned() + &x.to_string()[..])
                .collect(), // missing the 5th pool
            coverages: Array2::from_shape_vec(
                (5, 3),
                std::iter::repeat(100.0).take(5 * 3).collect(),
            )
            .unwrap(),
        };
        println!("genotypes_and_phenotypes={:?}", genotypes_and_phenotypes);
        println!(
            "genotypes_and_phenotypes.count_loci().unwrap()={:?}",
            genotypes_and_phenotypes.count_loci().unwrap()
        );
        let loci = genotypes_and_phenotypes.count_loci().unwrap();
        assert_eq!(vec![1, 4, 7, 9], loci.0); // includes the trailing locus position + 1 for index slicing since the end position is exclusive
        assert_eq!(
            vec![
                "chr1".to_owned(),
                "chr2".to_owned(),
                "chrX".to_owned(),
                "chrX".to_owned()
            ],
            loci.1
        ); // includes the trailing locus
        assert_eq!(vec![123, 50005, 701, 701], loci.2); // includes the trailing locus
                                                        // assert_eq!((), genotypes_and_phenotypes.count_loci().unwrap());
    }

    #[test]
    fn test_sync_parse() {
        let line = "Chromosome1\t456527\tC\t1:0:999:0:4:0\t0:1:2:0:0:0\t0:2:4:0:0:0\t0:1:4:0:0:0\t0:1:6:0:0:0".to_owned();
        let counts: Box<LocusCounts> = line.lparse().unwrap();
        println!("counts={:?}", counts);
        assert!(
            Array2::from_shape_vec(
                (5, 6),
                vec![
                    1, 0, 999, 0, 4, 0, 0, 1, 2, 0, 0, 0, 0, 2, 4, 0, 0, 0, 0, 1, 4, 0, 0, 0, 0, 1,
                    6, 0, 0, 0,
                ]
            )
            .unwrap()
                == (*counts).matrix
        );
    }

    #[test]
    fn test_sync_filter() {
        let line = "Chromosome1\t456527\tC\t1:0:999:0:4:0\t0:1:2:0:0:0\t0:2:4:0:0:0\t0:1:4:0:0:0\t0:1:6:0:0:0".to_owned();
        let mut counts: Box<LocusCounts> = line.lparse().unwrap();
        let mut frequencies: Box<LocusFrequencies> = counts.to_frequencies().unwrap();
        let filter_stats = FilterStats {
            remove_ns: true,
            min_quality: 0.005,
            min_coverage: 1,
            min_allele_frequency: 0.005,
            max_missingness_rate: 0.0,
            pool_sizes: vec![20., 20., 20., 20., 20.],
        };
        assert!(counts == counts.to_counts().unwrap());
        assert!(
            Array2::from_shape_vec(
                (5, 6),
                vec![
                    1. / 1004.,
                    0.,
                    999. / 1004.,
                    0.,
                    4. / 1004.,
                    0.,
                    0.,
                    1. / 3.,
                    2. / 3.,
                    0.,
                    0.,
                    0.,
                    0.,
                    2. / 6.,
                    4. / 6.,
                    0.,
                    0.,
                    0.,
                    0.,
                    1. / 5.,
                    4. / 5.,
                    0.,
                    0.,
                    0.,
                    0.,
                    1. / 7.,
                    6. / 7.,
                    0.,
                    0.,
                    0.,
                ]
            )
            .unwrap()
                == frequencies.matrix
        );
        counts.filter(&filter_stats).unwrap();
        assert!(
            Array2::from_shape_vec((5, 2), vec![0, 999, 1, 2, 2, 4, 1, 4, 1, 6]).unwrap()
                == counts.matrix
        );
        assert!(
            Array2::from_shape_vec(
                (5, 6),
                vec![
                    1, 0, 999, 0, 4, 0, 0, 1, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 1, 4, 0, 0, 0, 0, 1,
                    6, 0, 0, 0,
                ]
            )
            .unwrap()
                == frequencies.to_counts().unwrap().matrix
        );
        assert!(frequencies == frequencies.to_frequencies().unwrap());
        frequencies.filter(&filter_stats).unwrap();
        assert!(
            Array2::from_shape_vec(
                (5, 2),
                vec![
                    0. / 999.,
                    999. / 999.,
                    1. / 3.,
                    2. / 3.,
                    2. / 6.,
                    4. / 6.,
                    1. / 5.,
                    4. / 5.,
                    1. / 7.,
                    6. / 7.,
                ]
            )
            .unwrap()
                == frequencies.matrix
        );
    }

    #[test]
    fn test_sync_sort() {
        let line = "Chromosome1\t456527\tC\t1:0:999:0:4:0\t0:1:2:0:0:0\t0:2:4:0:0:0\t0:1:4:0:0:0\t0:1:6:0:0:0".to_owned();
        let counts: Box<LocusCounts> = line.lparse().unwrap();
        let mut frequencies: Box<LocusFrequencies> = counts.to_frequencies().unwrap();
        frequencies.sort_by_allele_freq(true).unwrap();
        assert!(vec!["C".to_owned(), "T".to_owned()] == frequencies.alleles_vector[0..2]);
        assert!(Array1::from_vec(vec![0.0; 5]).view() == frequencies.matrix.column(5));
        frequencies.sort_by_allele_freq(false).unwrap();
        assert!(vec!["T".to_owned(), "C".to_owned()] == frequencies.alleles_vector[4..6]);
        assert!(Array1::from_vec(vec![0.0; 5]).view() == frequencies.matrix.column(0));
    }

    #[test]
    fn test_sync_remove_missing() {
        let lines = vec!["Chromosome1\t12345\tC\t1:0:999:0:4:0\t0:1:2:0:0:0\t0:2:4:0:0:0\t0:1:4:0:0:0\t0:1:6:0:0:0".to_owned(),
        "Chromosome2\t6789\tA\t1:0:0:3:0:0\t1:0:0:5:0:0\t7:0:0:2:0:0\t2:0:0:10:0:0\t10:0:0:50:0:0".to_owned(),
        "Chromosome3\t42069\tT\t10:10:0:0:0:0\t5:12:0:0:0:0\t5:24:0:0:0:0\t5:14:0:0:0:0\t3:16:0:0:0:0".to_owned(),
        ];
        let mut counts: LocusCounts = *(lines[0].lparse().unwrap());
        let filter_stats = FilterStats {
            remove_ns: true,
            min_quality: 0.005,
            min_coverage: 1,
            min_allele_frequency: 0.005,
            max_missingness_rate: 0.0,
            pool_sizes: vec![20., 20., 20., 20., 20.],
        };
        counts.filter(&filter_stats).unwrap();
        let mut locus_counts_and_phenotypes = LocusCountsAndPhenotypes {
            locus_counts: counts,
            phenotypes: Array2::from_shape_vec((5, 1), vec![0.1, 0.2, 0.3, 0.4, 0.5]).unwrap(), // n pools x k traits
            pool_names: (0..5)
                .map(|x| ("pool_".to_owned() + &x.to_string()).to_owned())
                .collect::<Vec<String>>(),
        };

        // Extract loci counts and frequencies
        let mut loci_counts: Vec<LocusCounts> = vec![];
        let mut loci_frequencies: Vec<LocusFrequencies> = vec![];
        for i in 0..lines.len() {
            let mut counts: LocusCounts = *(lines[i].lparse().unwrap());
            let mut frequencies = counts.to_frequencies().unwrap();
            counts.filter(&filter_stats).unwrap();
            frequencies.filter(&filter_stats).unwrap();
            frequencies.filter(&filter_stats).unwrap();
            loci_counts.push(counts.clone());
            loci_frequencies.push(*frequencies);
        }
        // Extract coverages n x l loci
        let n = loci_counts[0].matrix.nrows();
        let p = loci_counts.len();
        let mut coverages: Array2<f64> = Array2::from_elem((n, p), f64::NAN);
        for j in 0..p {
            for i in 0..n {
                coverages[(i, j)] = loci_counts[j].matrix.row(i).sum() as f64;
            }
        }
        // println!("coverages={:?}", coverages);
        // Build the GenotypesAndPhenotypes struct
        let mut chromosome: Vec<String> = vec![];
        let mut position: Vec<u64> = vec![];
        let mut allele: Vec<String> = vec![];
        let mut vec_intercept_and_allele_frequencies: Vec<f64> = vec![];
        let phenotypes: Array2<f64> =
            Array2::from_shape_vec((5, 1), vec![0.1, 0.2, 0.3, 0.4, 0.5]).unwrap();
        let pool_names: Vec<String> = (0..5)
            .map(|x| ("pool_".to_owned() + &x.to_string()).to_owned())
            .collect();
        let coverages: Array2<f64> = coverages;
        for j in 0..loci_frequencies.len() {
            let l = loci_frequencies[j].matrix.ncols();
            chromosome.append(&mut vec![loci_frequencies[j].chromosome.clone(); l]);
            position.append(&mut vec![loci_frequencies[j].position; l]);
            allele.append(&mut loci_frequencies[j].alleles_vector.clone());
        }
        for i in 0..n {
            vec_intercept_and_allele_frequencies.push(1.0); // append the intercept
            for j in 0..loci_frequencies.len() {
                vec_intercept_and_allele_frequencies.append(
                    &mut loci_frequencies[j]
                        .matrix
                        .row(i)
                        .map(|&x| x.to_owned())
                        .into_iter()
                        .collect::<Vec<f64>>(),
                );
            }
        }
        // println!("vec_intercept_and_allele_frequencies={:?}", vec_intercept_and_allele_frequencies);
        let p_plus_one = vec_intercept_and_allele_frequencies.len() / n;
        let intercept_and_allele_frequencies: Array2<f64> =
            Array2::from_shape_vec((n, p_plus_one), vec_intercept_and_allele_frequencies).unwrap();
        let mut genotypes_and_phenotypes = GenotypesAndPhenotypes {
            chromosome: chromosome,
            position: position,
            allele: allele,
            intercept_and_allele_frequencies: intercept_and_allele_frequencies,
            phenotypes: phenotypes,
            pool_names: pool_names,
            coverages: coverages,
        };
        // println!("loci_frequencies={:?}", loci_frequencies);
        // println!("loci_counts={:?}", loci_counts);
        // println!("genotypes_and_phenotypes={:?}", genotypes_and_phenotypes);
        assert_eq!(
            vec![1, 1, 2, 2, 3, 3]
                .iter()
                .map(|&x| "Chromosome".to_owned() + &x.to_string())
                .collect::<Vec<String>>(),
            genotypes_and_phenotypes.chromosome
        );
        assert_eq!(
            vec![12345, 12345, 6789, 6789, 42069, 42069],
            genotypes_and_phenotypes.position
        );
        assert_eq!(
            Array1::from_vec(vec![0.0, 1.0, 2.0, 1.0, 1.0]),
            &genotypes_and_phenotypes
                .intercept_and_allele_frequencies
                .column(1)
                * &genotypes_and_phenotypes.coverages.column(0)
        );
        assert_eq!(
            Array1::from_vec(vec![1.0, 1.0, 7.0, 2.0, 10.]),
            &genotypes_and_phenotypes
                .intercept_and_allele_frequencies
                .column(3)
                * &genotypes_and_phenotypes.coverages.column(1)
        );
        assert_eq!(
            Array1::from_vec(vec![10., 5.0, 5.0, 5.0, 3.0]),
            &genotypes_and_phenotypes
                .intercept_and_allele_frequencies
                .column(5)
                * &genotypes_and_phenotypes.coverages.column(2)
        );

        locus_counts_and_phenotypes.phenotypes[(0, 0)] = f64::NAN;
        genotypes_and_phenotypes.phenotypes[(1, 0)] = f64::NAN;
        locus_counts_and_phenotypes.remove_missing().unwrap();
        genotypes_and_phenotypes.remove_missing().unwrap();
        println!(
            "locus_counts_and_phenotypes={:?}",
            locus_counts_and_phenotypes
        );
        println!("genotypes_and_phenotypes={:?}", genotypes_and_phenotypes);

        assert_eq!(4, locus_counts_and_phenotypes.phenotypes.nrows());
        assert_eq!(4, locus_counts_and_phenotypes.locus_counts.matrix.nrows());
        assert_eq!(
            4,
            genotypes_and_phenotypes
                .intercept_and_allele_frequencies
                .nrows()
        );
        assert_eq!(4, genotypes_and_phenotypes.coverages.nrows());
    }

    #[test]
    fn test_sync_load() {
        let file_sync_phen = FileSyncPhen {
            filename_sync: "tests/test.sync".to_owned(),
            pool_names: (0..5)
                .map(|x| "Pop".to_owned() + &x.to_string()[..])
                .collect(),
            pool_sizes: vec![100.0; 5],
            phen_matrix: Array2::from_shape_vec((5, 1), vec![0.1, 0.2, 0.3, 0.4, 0.5]).unwrap(),
            test: "test".to_owned(),
        };
        let filter_stats = FilterStats {
            remove_ns: true,
            min_quality: 1.00,
            min_coverage: 0,
            min_allele_frequency: 0.01,
            max_missingness_rate: 1.0,
            pool_sizes: vec![20., 20., 20., 20., 20.],
        };
        // Per chunk
        let (locus_frequencies, locus_counts) = file_sync_phen
            .per_chunk_load(&0, &300, &filter_stats, false)
            .unwrap();
        assert!(3 == locus_frequencies.len());
        assert!(3 == locus_counts.len());
        assert_eq!(
            LocusFrequencies {
                chromosome: "Chromosome1".to_owned(),
                position: 1131021,
                alleles_vector: vec!["T".to_owned(), "C".to_owned()],
                matrix: Array2::from_shape_vec(
                    (5, 2),
                    vec![0.5, 0.5, 0.5, 0.5, 0.0, 1.0, 1. / 6., 5. / 6., 1.0, 0.0]
                )
                .unwrap()
            },
            locus_frequencies[2]
        );
        assert_eq!(
            LocusCounts {
                chromosome: "Chromosome1".to_owned(),
                position: 1041321,
                alleles_vector: vec!["T".to_owned(), "G".to_owned()],
                matrix: Array2::from_shape_vec((5, 2), vec![20, 1, 12, 0, 15, 1, 22, 1, 26, 1])
                    .unwrap()
            },
            locus_counts[1]
        );
        // Load
        let (locus_frequencies_load, locus_counts_load) =
            file_sync_phen.load(&filter_stats, false, &2).unwrap();
        assert_eq!(5745, locus_frequencies_load.len());
        assert_eq!(5745, locus_counts_load.len());
        assert_eq!(locus_frequencies[0], locus_frequencies_load[0]);
        assert_eq!(locus_frequencies[1], locus_frequencies_load[1]);
        assert_eq!(locus_frequencies[2], locus_frequencies_load[2]);
        assert_eq!(locus_counts[0], locus_counts_load[0]);
        assert_eq!(locus_counts[1], locus_counts_load[1]);
        assert_eq!(locus_counts[2], locus_counts_load[2]);
        // Convert into genotypes and phenotypes
        let genotypes_and_phenotypes: GenotypesAndPhenotypes = file_sync_phen
            .into_genotypes_and_phenotypes(&filter_stats, false, &2)
            .unwrap();
        assert_eq!(12669, genotypes_and_phenotypes.chromosome.len());
        assert_eq!(12669, genotypes_and_phenotypes.allele.len());
        assert_eq!(12669, genotypes_and_phenotypes.position.len());
        assert_eq!(
            (5, 12669),
            genotypes_and_phenotypes
                .intercept_and_allele_frequencies
                .dim()
        );
    }

    #[test]
    fn test_sync_save() {
        let file_sync_phen = FileSyncPhen {
            filename_sync: "tests/test.sync".to_owned(),
            pool_names: (0..5)
                .map(|x| "Pop".to_owned() + &x.to_string()[..])
                .collect(),
            pool_sizes: vec![100.0; 5],
            phen_matrix: Array2::from_shape_vec((5, 1), vec![0.1, 0.2, 0.3, 0.4, 0.5]).unwrap(),
            test: "test".to_owned(),
        };
        let filter_stats = FilterStats {
            remove_ns: true,
            min_quality: 1.00,
            min_coverage: 0,
            min_allele_frequency: 0.01,
            max_missingness_rate: 1.0,
            pool_sizes: vec![20., 20., 20., 20., 20.],
        };
        // File sync phen struct
        let fname_out_file_sync_phen_save = "tests/test-file_sync_phen-save.csv".to_owned();
        let _ = file_sync_phen.write_csv(&filter_stats, false, &fname_out_file_sync_phen_save, &2);
        let file = std::fs::File::open(fname_out_file_sync_phen_save).unwrap();
        let reader = std::io::BufReader::new(file);
        let mut vec_lines: Vec<String> = vec![];
        for line in reader.lines() {
            let line = line.unwrap();
            vec_lines.push(line);
        }
        assert_eq!(
            "#chr,pos,allele,Pop0,Pop1,Pop2,Pop3,Pop4".to_owned(),
            vec_lines[0]
        );
        assert_eq!(
            "Chromosome1,456527,T,0,0.333333,0.333333,0.2,0.142857".to_owned(),
            vec_lines[1]
        );
        assert_eq!(
            "Chromosome1,1131021,T,0.5,0.5,0,0.166667,1".to_owned(),
            vec_lines[5]
        );
        assert_eq!(
            "Chromosome1,25785190,C,0.615385,0.833333,0.777778,1,0.9375".to_owned(),
            vec_lines[100]
        );
        // Genotypes and phenotypes struct
        let genotypes_and_phenotypes: GenotypesAndPhenotypes = file_sync_phen
            .into_genotypes_and_phenotypes(&filter_stats, false, &2)
            .unwrap();
        let fname_out_genotypes_and_phenotypes_save =
            "tests/test-genotypes_and_phenotypes-save.csv".to_owned();
        let _ = genotypes_and_phenotypes.write_csv(
            &filter_stats,
            false,
            &fname_out_genotypes_and_phenotypes_save,
            &2,
        );
        let file = std::fs::File::open(fname_out_genotypes_and_phenotypes_save).unwrap();
        let reader = std::io::BufReader::new(file);
        let mut vec_lines: Vec<String> = vec![];
        for line in reader.lines() {
            let line = line.unwrap();
            vec_lines.push(line);
        }
        assert_eq!(
            "#chr,pos,allele,Pop0,Pop1,Pop2,Pop3,Pop4".to_owned(),
            vec_lines[0]
        );
        assert_eq!(
            "Chromosome1,456527,T,0,0.333333,0.333333,0.2,0.142857".to_owned(),
            vec_lines[1]
        );
        assert_eq!(
            "Chromosome1,1131021,T,0.5,0.5,0,0.166667,1".to_owned(),
            vec_lines[5]
        );
        assert_eq!(
            "Chromosome1,25785190,C,0.615385,0.833333,0.777778,1,0.9375".to_owned(),
            vec_lines[100]
        );
    }
}
