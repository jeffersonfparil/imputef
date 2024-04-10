//! Structs and Traits

use ndarray::prelude::*;
use std::io;

////////////////////////////////////////////////////////////////////////////////
/// # STRUCTS
////////////////////////////////////////////////////////////////////////////////

/// The alternative entry point genotype file struct, i.e. describing the genotype data in vcf format
/// Note: Make sure that the annotate the vcf file with allele depths, e.g. `bcftools mpileup -a AD ...`
/// - `filename` - filename of the vcf file (`*.vcf` or `*.vcf.gz`)
#[derive(Debug, Clone)]
pub struct FileVcf {
    pub filename: String,
}

/// The main genotype file struct used for most of the analyses
/// - `filename` - filename of the synchronised pileup file (`*.sync`)
/// - `test` - name of statistical test, i.e. "sync2csv", "fisher_exact_test", "chisq_test", "fst", "heterozygosity, "pearson_corr", "ols_iter", "ols_iter_with_kinship", "mle_iter", "mle_iter_with_kinship", "gwalpha", "genomic_prediction_cross_validation""
#[derive(Debug, Clone)]
pub struct FileSync {
    pub filename: String,
    pub test: String,
}

/// The simple genotype file which is the only output format at the moment
/// - `filename` - filename of the tab-delimited genotype file (*.txt)
#[derive(Debug, Clone)]
pub struct FileGeno {
    pub filename: String,
}

/// Non-critical in this project as we do not need phenotype data for imputation.
/// This is reserved for quantitative and population genetics analyses in poolgen.
/// Filename of the phenotype file which can be a simple delimited file (e.g. csv,  tsv or even ssv; semi-colon-delimited) or a specialised GWAlpha phenotype information file in a python file.
/// - `filename` - filename of the phenotype file (e.g. `*.csv`, `*.txt`, `*.tsv`, or `*.ssv`)
/// - `delim` - string delimiter of the phenotype file (e.g. `","` or `"\t"`)
/// - `names_column_id` - index of the column containing the names of the pools or populations
/// - `sizes_column_id` - index of the column containing the sizes of each pool or population
/// - `trait_values_column_ids` - vector of indexes corresponding to the column containing the trait values to be included in the analyses (Note that multi-trait analyses may not be available to all analyses types)
/// - `format` - string defining the format of the phenotype file as `default` for simple delimited file or `gwalpha_fmt` for GWAlpha-required format (for back-compatibility with github.com/aflevel/GWAlpha/GWAlpha.py)
#[derive(Debug, Clone)]
pub struct FilePhen {
    pub filename: String,
    pub delim: String,
    pub names_column_id: usize,
    pub sizes_column_id: usize,
    pub trait_values_column_ids: Vec<usize>,
    pub format: String,
}

/// Non-critical in this project as we do not need phenotype data for imputation.
/// This is reserved for quantitative and population genetics analyses in poolgen.
/// Phenotype data including the names of the pools, the size of each pool, and the trait values
#[derive(Debug, Clone, PartialEq)]
pub struct Phen {
    pub pool_names: Vec<String>,
    pub pool_sizes: Vec<f64>,
    pub phen_matrix: Array2<f64>,
}

/// Filename of the synchronised pileup file and its corresponding phenotype data
/// Note that in this project, the `phen_matrix` field is not used because imputation does not require phenotype data. This is reserved for quantitative and population genetics analyses in poolgen.
#[derive(Debug, Clone, PartialEq)]
pub struct FileSyncPhen {
    pub filename_sync: String,
    pub pool_names: Vec<String>,
    pub pool_sizes: Vec<f64>,
    pub phen_matrix: Array2<f64>,
    pub test: String,
}

/// Locus and allele filtering settings
#[derive(Debug, Clone)]
pub struct FilterStats {
    pub remove_ns: bool,
    pub min_quality: f64,
    pub min_coverage: u64,
    pub min_allele_frequency: f64,
    pub max_missingness_rate: f64,
    pub pool_sizes: Vec<f64>,
}

/// A line of a vcf file corresponding to a single locus across all the pools
/// We are interested in extracting allele counts.
/// We are not interested in the genotype calls and their corresponding likelihoods.
#[derive(Debug, Clone, PartialEq)]
pub struct VcfLine {
    pub chromosome: String,             // chromosome or scaffold name
    pub position: u64,                  // position in number of bases
    pub reference_allele: char,         // reference allele
    pub alternative_alleles: Vec<char>, // vector of alternative alleles
    pub allele_depths: Vec<Vec<u64>>, // across samples average utf8 base quality codes which can be transformed into bases error rate as 10^(-(u8 - 33)/10)
}

/// Allele counts at a locus across pools
#[derive(Debug, Clone, PartialEq)]
pub struct LocusCounts {
    pub chromosome: String,
    pub position: u64,
    pub alleles_vector: Vec<String>,
    pub matrix: Array2<u64>, // n pools x p alleles
}

// Struct of allele frequencies to convert reads into sync
// #[derive(Debug, Clone, PartialEq, PartialOrd)]
#[derive(Debug, Clone, PartialEq)]
pub struct LocusFrequencies {
    pub chromosome: String,
    pub position: u64,
    pub alleles_vector: Vec<String>,
    pub matrix: Array2<f64>, // n pools x p alleles
}

// Struct of allele counts and phenotypes per pool
/// Note that in this project, the `phenotypes` field is not used because imputation does not require phenotype data. This is reserved for quantitative and population genetics analyses in poolgen.
#[derive(Debug, Clone)]
pub struct LocusCountsAndPhenotypes {
    pub locus_counts: LocusCounts,
    pub phenotypes: Array2<f64>, // n pools x k traits
    pub pool_names: Vec<String>,
}

// Struct of allele frequencies and phenotypes
/// Note that in this project, the `phenotypes` field is not used because imputation does not require phenotype data. This is reserved for quantitative and population genetics analyses in poolgen.#[derive(Debug, Clone)]
pub struct GenotypesAndPhenotypes {
    pub chromosome: Vec<String>,                       // 1 + p
    pub position: Vec<u64>,                            // 1 + p
    pub allele: Vec<String>,                           // 1 + p
    pub intercept_and_allele_frequencies: Array2<f64>, // n pools x 1 + p alleles across loci
    pub phenotypes: Array2<f64>,                       // n pools x k traits
    pub pool_names: Vec<String>,                       // n
    pub coverages: Array2<f64>,                        // n pools x m loci
}

// Struct of parameters chosen to identify linked loci and nearest neighbours
#[derive(Debug, Clone)]
pub struct ParametersOfLinkedLociAndNearestNeighbours {
    pub min_loci_corr: f64,
    pub max_pool_dist: f64,
    pub l_linked_loci: f64,
    pub k_nearest_neighbours: f64,
}

////////////////////////////////////////////////////////////////////////////////
/// # TRAITS
////////////////////////////////////////////////////////////////////////////////

pub trait CheckStruct {
    fn check(&self) -> io::Result<()>;
}

pub trait Count {
    fn count_loci(&self) -> io::Result<(Vec<usize>, Vec<String>, Vec<u64>)>;
}

pub trait Parse<T> {
    fn lparse(&self) -> io::Result<Box<T>>;
}

pub trait Filter {
    fn to_counts(&self) -> io::Result<Box<LocusCounts>>;
    fn to_frequencies(&self) -> io::Result<Box<LocusFrequencies>>;
    fn filter(&mut self, filter_stats: &FilterStats) -> io::Result<&mut Self>;
}

pub trait Sort {
    fn sort_by_allele_freq(&mut self, decreasing: bool) -> io::Result<&mut Self>;
}

pub trait RemoveMissing {
    fn remove_missing(&mut self) -> io::Result<&mut Self>;
}

pub trait ChunkyReadAnalyseWrite<T, F> {
    fn per_chunk(
        &self,
        start: &u64,
        end: &u64,
        outname_ndigits: &usize,
        filter_stats: &FilterStats,
        function: F,
    ) -> io::Result<String>
    where
        F: Fn(&mut T, &FilterStats) -> Option<String>;
    fn read_analyse_write(
        &self,
        filter_stats: &FilterStats,
        out: &str,
        n_threads: &usize,
        function: F,
    ) -> io::Result<String>
    where
        F: Fn(&mut T, &FilterStats) -> Option<String>;
}

pub trait LoadAll {
    fn per_chunk_load(
        &self,
        start: &u64,
        end: &u64,
        filter_stats: &FilterStats,
        keep_n_minus_1: bool,
    ) -> io::Result<(Vec<LocusFrequencies>, Vec<LocusCounts>)>; // Allele frequencies and counts across pools and alleles per locus
    fn load(
        &self,
        filter_stats: &FilterStats,
        keep_n_minus_1: bool,
        n_threads: &usize,
    ) -> io::Result<(Vec<LocusFrequencies>, Vec<LocusCounts>)>; // Allele frequencies and counts across pools and alleles per locus
    fn convert_into_genotypes_and_phenotypes(
        &self,
        filter_stats: &FilterStats,
        keep_n_minus_1: bool,
        n_threads: &usize,
    ) -> io::Result<GenotypesAndPhenotypes>;
}

pub trait SaveCsv {
    fn write_tsv(
        &self,
        filter_stats: &FilterStats,
        keep_p_minus_1: bool,
        out: &str,
        n_threads: &usize,
    ) -> io::Result<String>;
}
