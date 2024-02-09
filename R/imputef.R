#' @title
#' mvi
#' @description
#' Mean value imputation of allele frequencies
#' @usage
#' mvi(fname, 
#'        min_coverage=0,
#'        min_allele_frequency=0.0001,
#'        max_missingness_rate_per_locus=1.00,
#'        pool_sizes=c(100),
#'        min_depth_below_which_are_missing=1,
#'        max_depth_above_which_are_missing=1000000,
#'        frac_top_missing_pools=0.0,
#'        frac_top_missing_loci=0.0,
#'        n_threads=2,
#'        fname_out_prefix="")
#' @param fname
#' name of the genotype file to be imputed in uncompressed vcf, sync or allele frequency table format. See genotype format details below.
#' @param min_coverage
#' minimum coverage per locus, i.e. if at a locus, a pool falls below this value (does not skip missing data, i.e. missing locus has a depth of zero), then the whole locus is omitted. Set this to zero if the vcf has been filtered and contains missing values, i.e. `./.` or `.|.`. [Default=0]
#' @param min_allele_frequency
#' minimum allele frequency per locus, i.e. if at a locus, a pool has all its alleles below this value and/or above the additive complement of this value (skipping missing data), then the entire locus is omitted. [Default=0.0001]
#' @param max_missingness_rate_per_locus
#' maximum fraction of pools missing per locus, i.e. if at a locus, there were more pools missing than the coverage dictated by this threshold, then the locus is omitted. [Default=1.00]
#' @param pool_sizes
#' vector of pool sizes, i.e. the number of individuals included in each pool, or can be set to an arbitrarily large value like 100 for individual polyploids or if allele frequency estimates are expected to be accurate. [Default=100]
#' @param min_depth_below_which_are_missing
#' minimum depth at which loci with depth below this threshold are set to missing. Set to one if the input vcf has already been filtered and the loci beyond the depth thresholds have been set to missing, otherwise set to an integer above zero. [Default=1]
#' @param max_depth_above_which_are_missing
#' maximum depth at which loci with depth above this threshold are set to missing. Set to some large arbitrarily large value (e.g. 1000000) if the input vcf has already been filtered and the loci beyond the depth thresholds have been set to missing, otherwise set to an integer above zero. [Default=1000000]
#' @param frac_top_missing_pools
#' fraction of pools with the highest number of missing loci to be omitted. Set to zero if the input vcf has already been filtered and the loci beyond the depth thresholds have been set to missing, otherwise set to a decimal number between zero and one. [Default=0.0]
#' @param frac_top_missing_loci
#' fraction of loci with the highest number of pools with missing data to be omitted. Set to zero if the input vcf has already been filtered and the loci beyond the depth thresholds have been set to missing, otherwise set to an decimal number between zero and one. [Default=0.0]
#' @param n_threads
#' number of computing threads or processor cores to use in the computations. [Default=2]
#' @param fname_out_prefix
#' prefix of the output files including the [imputed allele frequency table](#allele-frequency-table-csv) (`<fname_out_prefix>-<time>-<random_id>-IMPUTED.csv`). [Default="" which will use the name of the input file as the prefix including the path]
#' @details
#' This mean value imputation method simply uses the arithmetic mean of the allele frequencies across all non-missing samples to impute missing data.
#' This function prints out the expected mean absolute error (MAE) of the imputation using 10% simulated missing data. Repeat the imputation manually to get a range of these MAEs for a better estimate of the expected MAE.
#' The allele depth information (AD), i.e. the unfiltered allele depth which includes the reads which did not pass the variant caller filters are used to calculate allele frequencies. 
#' Optional filtering steps based on minimum depth, minimum allele frequency, and maximum sparsity are available. 
#' If the GT field is present but the AD field is absent, then each sample is assumed to be an individual diploid, i.e., neither a polyploid nor a pool. 
#' Genotype data are not imported into R, LD estimation and imputation per se are multi-threaded, and imputation output is written into disk as an allele frequency table. 
#' The structs, traits, methods, and functions defined in this library are subsets of [poolgen](https://github.com/jeffersonfparil/poolgen), and will eventually be merged. 
#' Genotype file formats:
#' - vcf: canonical variant calling or genotype data format for individual samples. This should include the `AD` field (allele depth), and may or may not have genotypes called (e.g. generated via bctools mpileup -a AD,DP ...). If the `GT` field is present but the `AD` field is absent, then each sample is assumed to be an individual diploid, i.e., neither a polyploid nor a pool. The [`vcf2sync`](#vcf2sync) utility is expected to work with vcf versions 4.2 and 4.3. See [VCFv4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf) and [VCFv4.3](https://samtools.github.io/hts-specs/VCFv4.3.pdf) for details in the format specifications.
#' - sync: an extension of [popoolation2's](https://academic.oup.com/bioinformatics/article/27/24/3435/306737) sync or synchronised pileup file format, which includes a header line prepended by '#' showing the names of each column including the names of each pool. Additional header line/s and comments prepended with '#' may be added anywhere within the file.
#'    + tab-delimited
#'    + *Header line/s*:  optional header line/s including the names of the pools, e.g. `# chr pos ref pool1 pool2 pool3 pool4 pool5`
#'    + *Column 1*:       chromosome or scaffold name
#'    + *Column 2*:       locus position 
#'    + *Column 3*:       reference allele, e.g. A, T, C, G 
#'    + *Column/s 4 to n*:  colon-delimited allele counts: A:T:C:G:DEL:N, where "DEL" refers to insertion/deletion, and "N" is unclassified. A pool or population or polyploid individual is represented by a single column of this colon-delimited allele counts.
#' - allele frequency table
#'    + comma-delimited
#'    + *Header line*: ` #chr,pos,allele,<pool_name_1>,...,<pool_name_n>`
#'    + each locus is represented by 2 or more rows, i.e. 2 for biallelic loci, and >2 for multi-allelic loci
#' @returns 
#' imputed allele frequency table with the following filename: `<fname_out_prefix>-<time>-<random_id>-IMPUTED.csv`. Additionally, a sync file is generated if the input was an uncompressed vcf.
#' @export
mvi = function(fname, 
                min_coverage=0,
                min_allele_frequency=0.0001,
                max_missingness_rate_per_locus=1.00,
                pool_sizes=c(100),
                min_depth_below_which_are_missing=1,
                max_depth_above_which_are_missing=1000000,
                frac_top_missing_pools=0.0,
                frac_top_missing_loci=0.0,
                n_threads=2,
                fname_out_prefix="") {
    out = impute(fname=fname,
        imputation_method="mean",
        min_coverage=min_coverage,
        min_allele_frequency=min_allele_frequency,
        max_missingness_rate_per_locus=max_missingness_rate_per_locus,
        pool_sizes=pool_sizes,
        min_depth_below_which_are_missing=min_depth_below_which_are_missing,
        max_depth_above_which_are_missing=max_depth_above_which_are_missing,
        frac_top_missing_pools=frac_top_missing_pools,
        frac_top_missing_loci=frac_top_missing_loci,
        min_loci_corr=0,
        max_pool_dist=0,
        min_l_loci=0,
        min_k_neighbours=0,
        restrict_linked_loci_per_chromosome=FALSE,
        n_reps=0,
        n_threads=n_threads,
        fname_out_prefix=fname_out_prefix)
    return(out)
}

#' @title
#' aldknni
#' @description
#' Adaptive linkage-informed k-nearest neighbour imputation of allele frequencies
#' @usage
#' aldknni(fname, 
#'         min_coverage=0,
#'         min_allele_frequency=0.0001,
#'         max_missingness_rate_per_locus=1.00,
#'         pool_sizes=c(100),
#'         min_depth_below_which_are_missing=1,
#'         max_depth_above_which_are_missing=1000000,
#'         frac_top_missing_pools=0.0,
#'         frac_top_missing_loci=0.0,
#'         min_loci_corr=0.9,
#'         max_pool_dist=0.1,
#'         min_l_loci=1,
#'         min_k_neighbours=1,
#'         restrict_linked_loci_per_chromosome=TRUE,
#'         n_reps=20,
#'         n_threads=2,
#'         fname_out_prefix="")
#' @param fname
#' name of the genotype file to be imputed in uncompressed vcf, sync or allele frequency table format. See genotype format details below.
#' @param min_coverage
#' minimum coverage per locus, i.e. if at a locus, a pool falls below this value (does not skip missing data, i.e. missing locus has a depth of zero), then the whole locus is omitted. Set this to zero if the vcf has been filtered and contains missing values, i.e. `./.` or `.|.`. [Default=0]
#' @param min_allele_frequency
#' minimum allele frequency per locus, i.e. if at a locus, a pool has all its alleles below this value and/or above the additive complement of this value (skipping missing data), then the entire locus is omitted. [Default=0.0001]
#' @param max_missingness_rate_per_locus
#' maximum fraction of pools missing per locus, i.e. if at a locus, there were more pools missing than the coverage dictated by this threshold, then the locus is omitted. [Default=1.00]
#' @param pool_sizes
#' vector of pool sizes, i.e. the number of individuals included in each pool, or can be set to an arbitrarily large value like 100 for individual polyploids or if allele frequency estimates are expected to be accurate. [Default=100]
#' @param min_depth_below_which_are_missing
#' minimum depth at which loci with depth below this threshold are set to missing. Set to one if the input vcf has already been filtered and the loci beyond the depth thresholds have been set to missing, otherwise set to an integer above zero. [Default=1]
#' @param max_depth_above_which_are_missing
#' maximum depth at which loci with depth above this threshold are set to missing. Set to some large arbitrarily large value (e.g. 1000000) if the input vcf has already been filtered and the loci beyond the depth thresholds have been set to missing, otherwise set to an integer above zero. [Default=1000000]
#' @param frac_top_missing_pools
#' fraction of pools with the highest number of missing loci to be omitted. Set to zero if the input vcf has already been filtered and the loci beyond the depth thresholds have been set to missing, otherwise set to a decimal number between zero and one. [Default=0.0]
#' @param frac_top_missing_loci
#' fraction of loci with the highest number of pools with missing data to be omitted. Set to zero if the input vcf has already been filtered and the loci beyond the depth thresholds have been set to missing, otherwise set to an decimal number between zero and one. [Default=0.0]
#' @param min_loci_corr
#' Minimum correlation (Pearson's correlation) between the locus requiring imputation and other loci deemed to be in linkage with it. Ranges from 0.0 to 1.0. If using the default value with is NA, then this threshold will be optimised to find the best value minimising imputation error. [Default=NA]
#' @param max_pool_dist
#' Maximum genetic distance (mean absolute difference in allele frequencies) between the pool or sample requiring imputation and pools or samples deemed to be the closest neighbours. Ranges from 0.0 to 1.0. If using the default value with is NA, then this threshold will be optimised to find the best value minimising imputation error. [Default=NA]
#' @param min_l_loci
#' Minimum number of linked loci to be used in estimating genetic distances between the pool or sample requiring imputation and other pools or samples. Minimum value of 1. [Default=1]
#' @param min_k_neighbours
#' Minimum number of k-nearest neighbours of the pool or sample requiring imputation. Minimum value of 1. [Default=1]
#' @param restrict_linked_loci_per_chromosome
#' Restrict the choice of linked loci to within the chromosome the locus requiring imputation belong to? [Default=TRUE]
#' @param n_reps
#' Number of replications for the optimisation for the minimum loci correlation, and/or maximum genetic distance. Minimum value of 1. [Default=20]
#' @param n_threads
#' number of computing threads or processor cores to use in the computations. [Default=2]
#' @param fname_out_prefix
#' prefix of the output files including the [imputed allele frequency table](#allele-frequency-table-csv) (`<fname_out_prefix>-<time>-<random_id>-IMPUTED.csv`). [Default="" which will use the name of the input file as the prefix including the path]
#' @details
#' This is an attempt to extend the [LD-kNNi method of Money et al, 2015, i.e. LinkImpute](https://doi.org/10.1534/g3.115.021667), which was an extension of the [kNN imputation of Troyanskaya et al, 2001](https://doi.org/10.1093/bioinformatics/17.6.520). 
#' Similar to LD-kNNi, linkage disequilibrium (LD) is estimated using Pearson's product moment correlation per pair of loci, which is computed per chromosome by default, but can be computed across the entire genome. 
#' We use the mean absolute difference/error (MAE) between allele frequencies among linked loci as an estimate of genetic distance between samples. 
#' Fixed values for the minimum correlation to identify loci used in distance estimation, and maximum genetic distance to select the k-nearest neighbours can be defined. 
#' Additionally, minimum number of loci to include in distance estimation, and minimum number of nearest neighbours can be set. 
#' Moreover, all four parameters can be optimised, i.e. the minimum correlation and/or maximum distance and/or minimum number of loci and/or minimum number of nearest neighbours which minimises the MAE between predicted and expected allele frequencies after simulating 10% missing data are identified.
#' The allele depth information (AD), i.e. the unfiltered allele depth which includes the reads which did not pass the variant caller filters are used to calculate allele frequencies. 
#' If the GT field is present but the AD field is absent, then each sample is assumed to be an individual diploid, i.e., neither a polyploid nor a pool. 
#' Optional filtering steps based on minimum depth, minimum allele frequency, and maximum sparsity are available. 
#' Genotype data are not imported into R, LD estimation and imputation per se are multi-threaded, and imputation output is written into disk as an allele frequency table. 
#' The structs, traits, methods, and functions defined in this library are subsets of [poolgen](https://github.com/jeffersonfparil/poolgen), and will eventually be merged. 
#' Genotype file formats:
#' - vcf: canonical variant calling or genotype data format for individual samples. This should include the `AD` field (allele depth), and may or may not have genotypes called (e.g. generated via bctools mpileup -a AD,DP ...). If the `GT` field is present but the `AD` field is absent, then each sample is assumed to be an individual diploid, i.e., neither a polyploid nor a pool. The [`vcf2sync`](#vcf2sync) utility is expected to work with vcf versions 4.2 and 4.3. See [VCFv4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf) and [VCFv4.3](https://samtools.github.io/hts-specs/VCFv4.3.pdf) for details in the format specifications.
#' - sync: an extension of [popoolation2's](https://academic.oup.com/bioinformatics/article/27/24/3435/306737) sync or synchronised pileup file format, which includes a header line prepended by '#' showing the names of each column including the names of each pool. Additional header line/s and comments prepended with '#' may be added anywhere within the file.
#'    + tab-delimited
#'    + *Header line/s*:  optional header line/s including the names of the pools, e.g. `# chr pos ref pool1 pool2 pool3 pool4 pool5`
#'    + *Column 1*:       chromosome or scaffold name
#'    + *Column 2*:       locus position 
#'    + *Column 3*:       reference allele, e.g. A, T, C, G 
#'    + *Column/s 4 to n*:  colon-delimited allele counts: A:T:C:G:DEL:N, where "DEL" refers to insertion/deletion, and "N" is unclassified. A pool or population or polyploid individual is represented by a single column of this colon-delimited allele counts.
#' - allele frequency table
#'    + comma-delimited
#'    + *Header line*: ` #chr,pos,allele,<pool_name_1>,...,<pool_name_n>`
#'    + each locus is represented by 2 or more rows, i.e. 2 for biallelic loci, and >2 for multi-allelic loci
#' @returns 
#' imputed allele frequency table with the following filename: `<fname_out_prefix>-<time>-<random_id>-IMPUTED.csv`. Additionally, a sync file is generated if the input was an uncompressed vcf.


#' @export
aldknni = function(fname, 
                    min_coverage=0,
                    min_allele_frequency=0.0001,
                    max_missingness_rate_per_locus=1.00,
                    pool_sizes=c(100),
                    min_depth_below_which_are_missing=1,
                    max_depth_above_which_are_missing=1000000,
                    frac_top_missing_pools=0.0,
                    frac_top_missing_loci=0.0,
                    min_loci_corr=NA,
                    max_pool_dist=NA,
                    min_l_loci=1,
                    min_k_neighbours=1,
                    restrict_linked_loci_per_chromosome=TRUE,
                    n_reps=20,
                    n_threads=2,
                    fname_out_prefix="") {
    ### Handling NA conversion into Rust's f64:NAN
    if (is.na(min_loci_corr)) {
        min_loci_corr = -1.0
    }
    if (is.na(max_pool_dist)) {
        max_pool_dist = -1.0
    }
    out = impute(fname=fname,
        imputation_method="aLDkNNi",
        min_coverage=min_coverage,
        min_allele_frequency=min_allele_frequency,
        max_missingness_rate_per_locus=max_missingness_rate_per_locus,
        pool_sizes=pool_sizes,
        min_depth_below_which_are_missing=min_depth_below_which_are_missing,
        max_depth_above_which_are_missing=max_depth_above_which_are_missing,
        frac_top_missing_pools=frac_top_missing_pools,
        frac_top_missing_loci=frac_top_missing_loci,
        min_loci_corr=min_loci_corr,
        max_pool_dist=max_pool_dist,
        min_l_loci=min_l_loci,
        min_k_neighbours=min_k_neighbours,
        restrict_linked_loci_per_chromosome=restrict_linked_loci_per_chromosome,
        n_reps=n_reps,
        n_threads=n_threads,
        fname_out_prefix=fname_out_prefix)
    return(out)
}
