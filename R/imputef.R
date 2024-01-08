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
#' filename of the uncompressed vcf file which has the AD (allele depth) field and may or may not have genotypes called (e.g. generated via `bctools mpileup -a AD,DP ...`). If the GT field is present but the AD field is absent, then each sample is assumed to be an individual diploid, i.e., neither a polyploid nor a pool.
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
#' prefix of the output files (sync and csv files). [Default="" which will use the entire name of the input vcf file as the prefix]
#' @details
#' This mean value imputation method simply uses the arithmetic mean of the allele frequencies across all non-missing samples to impute missing data.
#' This function prints out the expected mean absolute error (MAE) of the imputation using 10% simulated missing data. Repeat the imputation manually to get a range of these MAEs for a better estimate of the expected MAE.
#' This function does not import any genotype data into R. Most processes are multi-threaded and outputs are written into disk as text files, i.e. sync file, and csv file of imputed allele frequencies. 
#' It converts vcf into sync with locus filtering based on minimum depth, minimum allele frequency, and maximum missingness rate with minimal memory footprint as large vcf files are split into chunks equivalent to the number of threads and processed line-by-line.
#' The allele depth information (AD), i.e. the unfiltered allele depth which includes the reads which did not pass the variant caller filters are used to calculate allele frequencies. 
#' The genotype calls (GT), if present are not used, hence the variant caller filtering is unimportant as only the allele frequencies are extracted from the the vcf file. 
#' However, if the AD information is absent then the GT information will be used assuming that each sample is an individual diploid, i.e., neither a polyploid nor a pool.
#' The entire sync file is then loaded into memory and imputed in parallel across windows. 
#' The structs, traits, methods, and functions defined in this library are subsets of [poolgen](https://github.com/jeffersonfparil/poolgen), and will eventually be merged with it. 
#' @returns sync file:
#' - filename: <fname_out_prefix>-<time>.sync
#' - an extension of [popoolation2's](https://academic.oup.com/bioinformatics/article/27/24/3435/306737) sync or synchronised pileup file format, which includes a header line prepended by '#' showing the names of each column including the names of each pool. Additional header line/s and comments prepended with '#' may be added anywhere within the file.
#'  + *Header line/s*:  optional header line/s including the names of the pools, e.g. `# chr pos ref pool1 pool2 pool3 pool4 pool5`
#'  + *Column 1*:       chromosome or scaffold name
#'  + *Column 2*:       locus position 
#'  + *Column 3*:       reference allele, e.g. A, T, C, G 
#'  + *Column/s 4 to n*:  colon-delimited allele counts: A:T:C:G:DEL:N, where "DEL" refers to insertion/deletion, and "N" is unclassified. A pool or population or polyploid individual is represented by a single column of this colon-delimited allele counts.
#' @returns imputed allele frequencies:
#' - filename: <fname_out_prefix>-<time>-<random_id>-IMPUTED.csv
#' - comma-separated file
#' - header: ` #chr,pos,allele,<pool_names>,...`
#' - each locus is represented by 2 or more rows, i.e. 2 for biallelic loci, and >2 for multi-allelic loci
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
        window_size_bp=0,
        min_loci_per_window=0,
        min_loci_corr=0,
        max_pool_dist=0,
        optimise_for_thresholds=FALSE,
        optimise_n_steps_corr=0,
        optimise_n_steps_dist=0,
        optimise_n_reps=0,
        misc_min_l=0,
        misc_min_k=0,
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
#'         window_size_bp=0,
#'         min_loci_per_window=1,
#'         min_loci_corr=0.9,
#'         max_pool_dist=0.1,
#'         optimise_for_thresholds=TRUE,
#'         optimise_n_steps_corr=10,
#'         optimise_n_steps_dist=10,
#'         optimise_n_reps=1,
#'         misc_min_l=0,
#'         misc_min_k=0,
#'         n_threads=2,
#'         fname_out_prefix="")
#' @param fname
#' filename of the uncompressed vcf file which has the AD (allele depth) field and may or may not have genotypes called (e.g. generated via `bctools mpileup -a AD,DP ...`). If the GT field is present but the AD field is absent, then each sample is assumed to be an individual diploid, i.e., neither a polyploid nor a pool.
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
#' @param window_size_bp
#' non-overlapping window size in bases to be used in the adaptive linkage-informed k-nearest neighbour imputation (aldknni) of allele frequencies. By default, set to the length of the longest chromosome or scaffold. Alternatively, set to the expected linkage block size (e.g. 10,000,000 or 10 megabases) to maximise the use available computing threads and decrease computation time. [Default=0]
#' @param min_loci_per_window
#' minimum number of loci per window to be used in the imputation of allele frequencies. Windows which fail this threshold are omitted. If the default value of one is used, then at very sparse windows the mean value will likely be used for imputation. [Default=1]
#' @param min_loci_corr
#' minimum correlation between loci or number of linked loci. If the former is intended, then this value ranges between 0 and 1, while greater than or equal to 1 if otherwise. Definition 1: Minimum correlation between the locus requiring imputation and other loci deemed to be in linkage with it. Definition 2: The number of linked loci. The resulting linked loci will be used to estimate the distances between pools. If this is set to 0, then this value will be optimised - see additional parameters below. [Default=0.9]
#' @param max_pool_dist
#' maximum distance between the pool or number of nearest neighbours. If the former is intended, then this value ranges between 0 and 1, while greater than or equal to 1 if otherwise. Definition 1: Maximum distance between the pool requiring imputation and other pools deemed to be the closest neighbours. Definition 2: Number of nearest neighbours. The resulting close neighbours will be used to impute. The distance metric is the mean absolute difference between pools across the linked loci. If this is set to 0, then this value will be optimised - see additional parameters below. [Default=0.1]
#' @param optimise_for_thresholds
#' optimise for minimum correlation and maximum distance thresholds if TRUE, else optimise for the number of linked loci and nearest neighbours. [Default=TRUE]
#' @param optimise_n_steps_corr
#' number levels for the optimisation of the minimum loci correlation or number of linked loci. [Default=10]
#' @param optimise_n_steps_dist
#' number levels for the optimisation of the maximum pool distance or number of nearest neighbours. [Default=10]
#' @param optimise_n_reps
#' number of replications for the optimisation of the minimum loci correlation or number of linked loci and maximum pool distance or number of nearest neighbours. [Default=1]
#' @param misc_min_l
#' Minimum number of linked loci to be included in imputation if using minimum loci correlation threshold. If the default value of zero is used, then mean value imputation will be used if no loci passed the minimum correlation threshold [Default=0].
#' @param misc_min_k
#' Minimum number of nearest neighbours to be included in imputation if using maximum distance threshold. If the default value of zero is used, then mean value imputation will be used if no neighbours passed the maximum distance threshold [Default=0].
#' @param n_threads
#' number of computing threads or processor cores to use in the computations. [Default=2]
#' @param fname_out_prefix
#' prefix of the output files (sync and csv files). [Default="" which will use the entire name of the input vcf file as the prefix]
#' @details
#' This is an attempt to extend the [LD-kNNi method of Money et al, 2015, i.e. LinkImpute](https://doi.org/10.1534/g3.115.021667), which was an extension of the [kNN imputation of Troyanskaya et al, 2001](https://doi.org/10.1093/bioinformatics/17.6.520). 
#' Similar to LD-kNNi, LD is estimated using Pearson's product moment correlation across loci per pair of samples, but instead of computing this across all the loci, we divide the genome into windows which respect chromosomal/scaffold boundaries. 
#' We use the mean absolute difference (MAD or MAE where E stands for error) between allele frequencies as an estimate of distance between samples. 
#' Instead of optimising for the number of loci to include in the distance estimation and the number neighbours to include in the weighted allele frequency mean, we use a minimum correlation coefficient for the former, and a maximum distance for the latter.
#' Both of these parameters can range from 0 to 1 and can be separately optimised, but a single pair of reasonable values is expected to result in good imputation accuracy, e.g. the default values of 0.9 minimum correlation, and 0.1 maximum distance.
#' The adaptive behaviour of our algorithm can be described in cases where:
#'  - sparsity in the data is too high, or 
#'  - loci are too uncorrelated because the breadth of coverage is too sparse, or
#'  - pools are too unrelated.
#'
#' These cases can mean that the data may not be informative enough to yield LD-kNN imputations better than mean value imputation. 
#' Hence, under these cases mean value imputation will be used instead of LD-kNN imputation.
#' This function prints out the expected mean absolute error (MAE) of the imputation using 10% simulated missing data. Repeat the imputation manually to get a range of these MAEs for a better estimate of the expected MAE.
#' This function does not import any genotype data into R. Most processes are multi-threaded and outputs are written into disk as text files, i.e. sync file, and csv file of imputed allele frequencies. 
#' It converts vcf into sync with locus filtering based on minimum depth, minimum allele frequency, and maximum missingness rate with minimal memory footprint as large vcf files are split into chunks equivalent to the number of threads and processed line-by-line.
#' The allele depth information (AD), i.e. the unfiltered allele depth which includes the reads which did not pass the variant caller filters are used to calculate allele frequencies. 
#' The genotype calls, if present are not used, hence the variant caller filtering is unimportant as only the allele frequencies are extracted from the the vcf file. 
#' However, if the AD information is absent then the GT information will be used assuming that each sample is an individual diploid, i.e., neither a polyploid nor a pool.
#' The entire sync file is then loaded into memory and imputed in parallel across windows. 
#' The structs, traits, methods, and functions defined in this library are subsets of [poolgen](https://github.com/jeffersonfparil/poolgen), and will eventually be merged with it. 
#' @returns sync file:
#' - filename: <fname_out_prefix>-<time>.sync
#' - an extension of [popoolation2's](https://academic.oup.com/bioinformatics/article/27/24/3435/306737) sync or synchronised pileup file format, which includes a header line prepended by '#' showing the names of each column including the names of each pool. Additional header line/s and comments prepended with '#' may be added anywhere within the file.
#'  + *Header line/s*:  optional header line/s including the names of the pools, e.g. `# chr pos ref pool1 pool2 pool3 pool4 pool5`
#'  + *Column 1*:       chromosome or scaffold name
#'  + *Column 2*:       locus position 
#'  + *Column 3*:       reference allele, e.g. A, T, C, G 
#'  + *Column/s 4 to n*:  colon-delimited allele counts: A:T:C:G:DEL:N, where "DEL" refers to insertion/deletion, and "N" is unclassified. A pool or population or polyploid individual is represented by a single column of this colon-delimited allele counts.
#' @returns imputed allele frequencies:
#' - filename: <fname_out_prefix>-<time>-<random_id>-IMPUTED.csv
#' - comma-separated file
#' - header: ` #chr,pos,allele,<pool_names>,...`
#' - each locus is represented by 2 or more rows, i.e. 2 for biallelic loci, and >2 for multi-allelic loci
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
                    window_size_bp=0,
                    min_loci_per_window=1,
                    min_loci_corr=0.9,
                    max_pool_dist=0.1,
                    optimise_for_thresholds=TRUE,
                    optimise_n_steps_corr=10,
                    optimise_n_steps_dist=10,
                    optimise_n_reps=1,
                    misc_min_l=0,
                    misc_min_k=0,
                    n_threads=2,
                    fname_out_prefix="") {
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
        window_size_bp=window_size_bp,
        min_loci_per_window=min_loci_per_window,
        min_loci_corr=min_loci_corr,
        max_pool_dist=max_pool_dist,
        optimise_for_thresholds=optimise_for_thresholds,
        optimise_n_steps_corr=optimise_n_steps_corr,
        optimise_n_steps_dist=optimise_n_steps_dist,
        optimise_n_reps=optimise_n_reps,
        misc_min_l=misc_min_l,
        misc_min_k=misc_min_k,
        n_threads=n_threads,
        fname_out_prefix=fname_out_prefix)
    return(out)
}
