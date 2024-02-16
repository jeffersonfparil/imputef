### Load imputef library
# dir_src = "/group/pasture/Jeff/imputef"
# dir_src = "/home/jeff/imputef"

### Extract allele frequencies into a pxn matrix where we have p loci and n entries
### Assumes all loci have a maximum of 2 alleles
fn_extract_allele_frequencies = function(vcf, min_depth=10, max_depth=200) {
    vec_loci_names = paste(vcfR::getCHROM(vcf), vcfR::getPOS(vcf), vcfR::getREF(vcf), sep="_")
    vec_pool_names = colnames(vcf@gt)[-1]
    mat_allele_counts = vcfR::extract.gt(vcf, element="AD")
    ### Extract biallelic diploid allele frequencies if the AD field is missing
    if (sum(is.na(mat_allele_counts)) == prod(dim(mat_allele_counts))) {
        GT = vcfR::extract.gt(vcf, element="GT")
        mat_genotypes = matrix(NA, nrow=nrow(mat_allele_counts), ncol=ncol(mat_allele_counts))
        mat_genotypes[(GT == "0/0") | (GT == "0|0")] = 1.0
        mat_genotypes[(GT == "1/1") | (GT == "1|1")] = 0.0
        mat_genotypes[(GT == "0/1") | (GT == "0|1") | (GT == "1|0")] = 0.5
        mat_coverage = matrix(min_depth + (max_depth - min_depth)/2, nrow=nrow(mat_genotypes), ncol=ncol(mat_genotypes)) 
    } else {
        mat_ref_counts = vcfR::masplit(mat_allele_counts, delim=',', record=1, sort=0)
        mat_alt_counts = vcfR::masplit(mat_allele_counts, delim=',', record=2, sort=0)
        ### Set missing allele counts to 0, if the other allele is non-missing and non-zero
        idx_for_ref = which(!is.na(mat_alt_counts) & (mat_alt_counts != 0) & is.na(mat_ref_counts))
        idx_for_alt = which(!is.na(mat_ref_counts) & (mat_ref_counts != 0) & is.na(mat_alt_counts))
        mat_ref_counts[idx_for_ref] = 0
        mat_alt_counts[idx_for_alt] = 0
        ### Calculate reference allele frequencies
        mat_coverage = mat_ref_counts + mat_alt_counts
        mat_genotypes = mat_ref_counts / mat_coverage
    }
    ### Identify high-confidence loci for assessing accuracy (miscellaneous metrics where we will only compute accuracy metrics on these high-confidence loci)
    mat_idx_high_conf_data = (mat_coverage >= min_depth) & (mat_coverage <= max_depth)
    ### Label the loci and pools
    rownames(mat_genotypes) = vec_loci_names
    colnames(mat_genotypes) = vec_pool_names
    rownames(mat_idx_high_conf_data) = vec_loci_names
    colnames(mat_idx_high_conf_data) = vec_pool_names
    ### Output
    return(list(mat_genotypes=mat_genotypes, mat_idx_high_conf_data=mat_idx_high_conf_data))
}

### Identity genotypes and generate a genotype (allele frequency or genotype class) matrix
### This uses the ploidy level of the species to define these genotype classes, 
### e.g. for diploids we expect 3 genotype classes - AA, AB/BA, and BB, while for tetraploids we expect 5 genotype classes - AAAA, AAAB, AABB, ABBB, and BBBB.
### The default behaviour is to define strict boundaries for the extreme genotypes, i.e. we only consider genotypes to be homozygotes if the allele depth is fixed for one allele.
### Non-strict boundaries is reserved for imputation where we use weighted means and hence the probability of the imputed genotype to belong to one class or another is not strictly bounded at the extremes.
fn_classify_allele_frequencies = function(mat_genotypes, ploidy, strict_boundaries=TRUE) {
    if (strict_boundaries) {
        vec_expected_frequencies = c(0, c(0:ploidy))/ploidy
        ### We are only setting the unfixed frequencies
        for (i in 3:(ploidy+1)) {
            # i = 3
            dfreq0 = (vec_expected_frequencies[i-1] - vec_expected_frequencies[i-2])/2
            dfreq1 = (vec_expected_frequencies[i-0] - vec_expected_frequencies[i-1])/2
            # print(vec_expected_frequencies[i-1] + dfreq0)
            # print(vec_expected_frequencies[i+1] - dfreq1)
            idx = (!is.na(mat_genotypes) & 
                (mat_genotypes >  (vec_expected_frequencies[i-1] + dfreq0)) & 
                (mat_genotypes <= (vec_expected_frequencies[i-0] + dfreq1)))
            if (i == ploidy+1) {
                idx = (!is.na(mat_genotypes) & 
                    (mat_genotypes > (vec_expected_frequencies[i-1] + dfreq0)) & 
                    (mat_genotypes < 1.00))
            }
            mat_genotypes[idx] = vec_expected_frequencies[i]
        }
    } else {
        # dfreq = (1/ploidy)/2
        # vec_expected_frequencies = c(-dfreq, c(0:ploidy)/ploidy)
        # for (i in 2:(ploidy+2)) {
        #     idx = (!is.na(mat_genotypes)) &
        #           (mat_genotypes >= (vec_expected_frequencies[i-1] + dfreq)) &
        #           (mat_genotypes <  (vec_expected_frequencies[i-0] + dfreq))
        #     # print(vec_expected_frequencies[i-1] + dfreq)
        #     # print(vec_expected_frequencies[i-0] + dfreq)
        #     mat_genotypes[idx] = vec_expected_frequencies[i]
        # }
        mat_genotypes = round(mat_genotypes * ploidy) / ploidy
    }
    return(mat_genotypes)
}

### Simulate missing data (we are setting maximum possible missing rate at 90% of entries x loci massing the maf threshold)
fn_simulate_missing_data = function(vcf, mat_genotypes, maf=0.25, missing_rate=0.5) {
    # maf = 0.25
    # missing_rate = 0.5
    ### Filter by maf
    vec_mean_freqs = rowMeans(mat_genotypes)
    idx = which((vec_mean_freqs >= maf) & (vec_mean_freqs <= (1-maf)))
    vcf_filtered = vcf[idx, , ]
    mat_genotypes_filtered = mat_genotypes[idx, ]
    ### Genotype data dimensions
    n = ncol(mat_genotypes_filtered)
    p = nrow(mat_genotypes_filtered)
    ### Randomly sample unique missing loci
    idx = sort(sample.int(n*p, size=n*p*missing_rate, replace=FALSE))
    idx_row = ceiling(idx / n)
    idx_col = idx %% n
    idx_col[idx_col==0] = n
    idx_mat = cbind(idx_row, idx_col)
    ### Extract loci and pools names corresponding to the simulated missing data
    vec_loci_names = rownames(mat_genotypes_filtered)
    vec_pool_names = colnames(mat_genotypes_filtered)
    vec_missing_loci = vec_loci_names[idx_mat[,1]]
    vec_missing_pools = vec_pool_names[idx_mat[,2]]
    ### Extract the expected allele frequencies
    expected_allele_frequencies = mat_genotypes_filtered[idx_mat]
    ### Statistics
    print(paste0("Total number of pools = ", n))
    print(paste0("Total number of loci = ", p))
    print(paste0("Simulated missing data rate = ", missing_rate))
    print("Distribution of the expected allele frequencies masked/simulted to be missing:")
    txtplot::txtdensity(expected_allele_frequencies)
    ### Mask the simulated missing data
    idx_mat[,2] = 1 + idx_mat[,2] ### Skip the first column of GT, i.e. the format column
    vcf_filtered@gt[idx_mat] = NA
    print(vcf_filtered)
    ### Remove fixed loci
    idx_mat[,2] = idx_mat[,2] - 1 ### Un-skip the first column of GT, i.e. the format column
    mat_genotypes_filtered[idx_mat] = NA
    vec_var_freqs = apply(mat_genotypes_filtered, MARGIN=1, FUN=var, na.rm=TRUE)
    idx_non_fixed_loci_after_simulating_missing_data = which(vec_var_freqs > 0.0)
    vcf_filtered = vcf_filtered[idx_non_fixed_loci_after_simulating_missing_data, ]
    print(vcf_filtered)
    ### Save the vcf and unzip
    print("Saving and unzipping vcf with simulated missing data...")
    fname_vcf_gz = paste0("SIMULATED_MISSING-", missing_rate, "-", as.numeric(Sys.time()), ".vcf.gz")
    fname_vcf = gsub(".gz$", "", fname_vcf_gz)
    vcfR::write.vcf(vcf_filtered, file=fname_vcf_gz)
    system(paste0("gunzip -f ", fname_vcf_gz))
    ### Output
    return(list(
        vec_missing_loci=vec_missing_loci,
        vec_missing_pools=vec_missing_pools,
        expected_allele_frequencies=expected_allele_frequencies,
        fname_vcf=fname_vcf
    ))
}

### Subsample loci to artificially reduce marker density
fn_simulate_marker_density_reduction = function(vcf, mat_genotypes, reduction_rate=0.5) {
    p = dim(vcf)[1]
    if (p != nrow(mat_genotypes)) {
        print("The vcfR object and the matrix of genotypes have incompatible dimensions.")
        return(0)
    }
    idx = sort(sample(c(1:p), size=min(c(p, round(p*reduction_rate))), replace=FALSE))
    return(list(vcf=vcf[idx, ], mat_genotypes=mat_genotypes[idx, ]))
}

### Imputation accuracy metrics
fn_metrics = function(q_predicted, q_expected) {
    ## Overall metrics
    deviation = q_predicted - q_expected
    mae = mean(abs(deviation), na.rm=TRUE)
    mse = mean(deviation^2, na.rm=TRUE)
    rmse = sqrt(mse)
    r2 = 1.00 - (mean(mse, na.rm=TRUE) / mean((q_expected-mean(q_expected))^2, na.rm=TRUE))
    concordance = mean(q_predicted == q_expected, na.rm=TRUE)
    ### Metrics across the range of expected allele frequencies
    # vec_q_max = c(0.00, 0.01, 0.05, c(1:9)/10, 0.95, 0.99, 1.00)
    vec_q_max = seq(0, 1, by=0.1)
    vec_n = rep(0, each=length(vec_q_max))
    vec_mae = rep(NA, each=length(vec_q_max))
    vec_mse = rep(NA, each=length(vec_q_max))
    vec_rmse = rep(NA, each=length(vec_q_max))
    vec_r2 = rep(NA, each=length(vec_q_max))
    vec_concordance = rep(NA, each=length(vec_q_max))
    for (i in 1:length(vec_q_max)) {
        # i = 1
        if (i==1) {
            q_min = -1e-20
        } else {
            q_min = vec_q_max[i-1]
        }
        q_max = vec_q_max[i]
        idx = which((q_expected > q_min) & (q_expected <= q_max) & (is.na(q_expected)==FALSE) & (is.na(q_predicted)==FALSE))
        if (length(idx) > 0) {
            vec_n[i] = length(idx)
            vec_deviations = q_predicted[idx] - q_expected[idx]
            vec_mae[i] = mean(abs(vec_deviations), na.rm=TRUE)
            vec_mse[i] = mean(vec_deviations^2, na.rm=TRUE)
            vec_rmse[i] = sqrt(vec_mse[i])
            vec_r2[i] = 1.00 - (mean(vec_mse[i], na.rm=TRUE) / mean((q_expected[idx]-mean(q_expected[idx]))^2, na.rm=TRUE))
            vec_concordance[i] = mean(q_predicted[idx] == q_expected[idx], na.rm=TRUE)
        }
    }
    df_metrics_across_allele_freqs = data.frame(q=vec_q_max, n=vec_n, mae=vec_mae, mse=vec_mse, rmse=vec_rmse, r2=vec_r2, concordance=vec_concordance)
    return(list(
        mae=mae,
        mse=mse,
        rmse=rmse,
        r2=r2,
        concordance=concordance,
        df_metrics_across_allele_freqs=df_metrics_across_allele_freqs
    ))
}

### Assess imputation accuracies
fn_imputation_accuracy = function(fname_imputed, list_sim_missing, mat_idx_high_conf_data, ploidy=4, strict_boundaries=FALSE, lukes=FALSE, mat_genotypes=NULL, n_threads=10) {
    # fname_imputed = fname_out_mvi
    # # fname_imputed = fname_out_aldknni
    # ploidy = 2
    # lukes = FALSE
    # mat_genotypes = NULL
    # n_threads = 32
    ### Extract imputed allele frequencies corresponding to the expected allele frequencies
    n_missing = length(list_sim_missing$vec_missing_loci)
    ### Load imputed (and lightly filtered ;-P) genotype data
    df = read.csv(fname_imputed)
    imputed_loci_names = paste(df$X.chr, df$pos, df$allele, sep="_")
    imputed_pool_names = colnames(df)
    imputed_pool_names = gsub("-", ".", imputed_pool_names)
    print(paste0("Identifying the imputed data points in: ", fname_imputed, "..."))
    vec_imputed = parallel::mclapply(X=c(1:n_missing),
        FUN=function(i, list_sim_missing, imputed_loci_names){
            # i = 1 # for (i in 1:n_missing) {
            locus = list_sim_missing$vec_missing_loci[i]
            pool = list_sim_missing$vec_missing_pools[i]
            idx_locus = which(imputed_loci_names %in% locus)
            idx_pool = which(imputed_pool_names %in% gsub("-", ".", pool))
            if (length(idx_pool) == 0) {
                idx_pool = which(imputed_pool_names %in% gsub(":", ".", pool))
            }
            if ((length(idx_locus)==0) | (length(idx_pool)==0)) {
                out = NA
            } else {
                out = df[idx_locus, idx_pool]
            }
            return(out)
        }, list_sim_missing=list_sim_missing, imputed_loci_names=imputed_loci_names,
    mc.cores=n_threads)
    vec_imputed = unlist(vec_imputed)
    mat_high_conf_true_and_imputed = parallel::mclapply(X=c(1:n_missing),
        FUN=function(i, list_sim_missing, imputed_loci_names){
            # i = 1 # for (i in 1:n_missing) {
            locus = list_sim_missing$vec_missing_loci[i]
            pool = list_sim_missing$vec_missing_pools[i]
            idx_locus = which(imputed_loci_names %in% locus)
            idx_pool = which(imputed_pool_names %in% gsub("-", ".", pool))
            if (length(idx_pool) == 0) {
                idx_pool = which(imputed_pool_names %in% gsub(":", ".", pool))
            }
            idx_high_conf_data = mat_idx_high_conf_data[which(rownames(mat_idx_high_conf_data) %in% locus), which(colnames(mat_idx_high_conf_data) %in% pool)]
            if ((length(idx_locus)==0) | (length(idx_pool)==0) | (idx_high_conf_data==FALSE)) {
                out_conf = NA
            } else {
                out_conf = df[idx_locus, idx_pool]
            }
            # return(out)
            return(c(list_sim_missing$expected_allele_frequencies[i], out_conf))
        }, list_sim_missing=list_sim_missing, imputed_loci_names=imputed_loci_names,
    mc.cores=n_threads)
    mat_high_conf_true_and_imputed = matrix(unlist(mat_high_conf_true_and_imputed), byrow=TRUE, ncol=2)
    ### Missing stats
    n_imputed = sum(!is.na(vec_imputed))
    # print(paste0("Total number of simulated missing data = ", n_missing))
    # print(paste0("Total number of imputed missing data = ", n_imputed))
    ### Metrics using allele frequencies
    metrics_allele_frequencies = fn_metrics(q_predicted=vec_imputed, q_expected=list_sim_missing$expected_allele_frequencies)
    ### Metrics using genotype classes
    vec_expected_classes = fn_classify_allele_frequencies(mat_genotypes=list_sim_missing$expected_allele_frequencies, ploidy=ploidy, strict_boundaries=strict_boundaries)
    vec_imputed_classes = fn_classify_allele_frequencies(mat_genotypes=vec_imputed, ploidy=ploidy, strict_boundaries=strict_boundaries)
    metrics_genotype_classes = fn_metrics(q_predicted=vec_imputed_classes, q_expected=vec_expected_classes)
    ### Metrics using high-confidence data points
    metrics_allele_frequencies_high_conf = fn_metrics(q_predicted=mat_high_conf_true_and_imputed[,2], q_expected=mat_high_conf_true_and_imputed[,1])
    vec_expected_classes_high_conf = fn_classify_allele_frequencies(mat_genotypes=mat_high_conf_true_and_imputed[,1], ploidy=ploidy, strict_boundaries=strict_boundaries)
    vec_imputed_classes_high_conf = fn_classify_allele_frequencies(mat_genotypes=mat_high_conf_true_and_imputed[,2], ploidy=ploidy, strict_boundaries=strict_boundaries)
    metrics_genotype_classes_high_conf = fn_metrics(q_predicted=vec_imputed_classes_high_conf, q_expected=vec_expected_classes_high_conf)
    return(list(
        frac_imputed = n_imputed / n_missing,
        mae_frequencies = metrics_allele_frequencies$mae,
        rmse_frequencies = metrics_allele_frequencies$rmse,
        r2_frequencies = metrics_allele_frequencies$r2,
        mae_classes = metrics_genotype_classes$mae,
        rmse_classes = metrics_genotype_classes$rmse,
        r2_classes = metrics_genotype_classes$r2,
        concordance_classes = metrics_genotype_classes$concordance,
        df_metrics_across_allele_freqs_frequencies = metrics_allele_frequencies$df_metrics_across_allele_freqs,
        df_metrics_across_allele_freqs_classes = metrics_genotype_classes$df_metrics_across_allele_freqs,
        highConf_mae_frequencies = metrics_allele_frequencies_high_conf$mae,
        highConf_rmse_frequencies = metrics_allele_frequencies_high_conf$rmse,
        highConf_r2_frequencies = metrics_allele_frequencies_high_conf$r2,
        highConf_mae_classes = metrics_genotype_classes_high_conf$mae,
        highConf_rmse_classes = metrics_genotype_classes_high_conf$rmse,
        highConf_r2_classes = metrics_genotype_classes_high_conf$r2,
        highConf_concordance_classes = metrics_genotype_classes_high_conf$concordance,
        highConf_df_metrics_across_allele_freqs_frequencies = metrics_allele_frequencies_high_conf$df_metrics_across_allele_freqs,
        highConf_df_metrics_across_allele_freqs_classes = metrics_genotype_classes_high_conf$df_metrics_across_allele_freqs
    ))
}

### Performance assessment function
fn_test_imputation = function(vcf, mat_genotypes, mat_idx_high_conf_data, ploidy=4, maf=0.25, missing_rate=0.5, strict_boundaries=FALSE, n_threads=10) {
    # vcf = vcfR::read.vcfR(file.path(dir_src, "res/grape.vcf"))
    # list_genotypes = fn_extract_allele_frequencies(vcf)
    # mat_genotypes = list_genotypes$mat_genotypes
    # mat_idx_high_conf_data = list_genotypes$mat_idx_high_conf_data
    # maf = 0.01
    # missing_rate = 0.1
    # strict_boundaries = FALSE
    # n_threads = 32
    ## Simulate missing data
    list_sim_missing = fn_simulate_missing_data(vcf=vcf,
                                                mat_genotypes=mat_genotypes,
                                                maf=maf,
                                                missing_rate=missing_rate)
    ### Define the actual number of missing loci after simulating missing data to account for cases when missingness (sparsity) is above 90% which is the maximum sparsity we have artificially set for computationally efficiency and to avoid errors due to too much sparsity
    n_missing = length(list_sim_missing$vec_missing_loci)
    rand_number_id = sample.int(1e9, 1)
    ### Use all loci if the number of loci is less than 10k lco
    if (ncol(mat_genotypes) < 1e4) {
        restrict_linked_loci_per_chromosome = FALSE
    } else {
        restrict_linked_loci_per_chromosome = TRUE
    }
    print(paste0("restrict_linked_loci_per_chromosome=", restrict_linked_loci_per_chromosome))
    ### (1) Mean value imputation
    time_ini = Sys.time()
    fname_out_mvi = tail(system(paste0(dir_src, "/target/release/imputef -f ", list_sim_missing$fname_vcf, " -m='mean' --fname-out-prefix='MVI-maf", maf, "-missing_rate", missing_rate, "-", rand_number_id, "' --n-threads=", n_threads), intern=TRUE), n=1)
    duration_mvi = difftime(Sys.time(), time_ini, units="mins")
    ### (2) Adaptive LD-kNN imputation using default fixed min_loci_corr and max_pool_dist at 0.9 and 0.1, respectively (where min_l_loci=20 and min_k_neighbours=5)
    time_ini = Sys.time()
    if (restrict_linked_loci_per_chromosome) {
        fname_out_aldknni_fixed = tail(system(paste0(dir_src, "/target/release/imputef -f ", list_sim_missing$fname_vcf, " --fname-out-prefix='AFIXED-maf", maf, "-missing_rate", missing_rate, "-", rand_number_id, "' --restrict-linked-loci-per-chromosome --n-threads=", n_threads), intern=TRUE), n=1)
    } else {
        fname_out_aldknni_fixed = tail(system(paste0(dir_src, "/target/release/imputef -f ", list_sim_missing$fname_vcf, " --fname-out-prefix='AFIXED-maf", maf, "-missing_rate", missing_rate, "-", rand_number_id, "' --n-threads=", n_threads), intern=TRUE), n=1)
    }
    duration_aldknni_fixed = difftime(Sys.time(), time_ini, units="mins")
    ### (3) Adaptive LD-kNN imputation using optimised min_loci_corr and max_pool_dist (where min_l_loci=20 and min_k_neighbours=5)
    time_ini = Sys.time()
    if (restrict_linked_loci_per_chromosome) {
        tmp_fname_out_aldknni_optim = system(paste0(dir_src, "/target/release/imputef -f ", list_sim_missing$fname_vcf, " --fname-out-prefix='AOPTIM-maf", maf, "-missing_rate", missing_rate, "-", rand_number_id, "' --restrict-linked-loci-per-chromosome --n-threads=", n_threads, " --min-loci-corr='-1.0' --max-pool-dist='-1.0'"), intern=TRUE)
        print(tmp_fname_out_aldknni_optim)
        fname_out_aldknni_optim = tail(tmp_fname_out_aldknni_optim, n=1)
    } else {
        tmp_fname_out_aldknni_optim = system(paste0(dir_src, "/target/release/imputef -f ", list_sim_missing$fname_vcf, " --fname-out-prefix='AOPTIM-maf", maf, "-missing_rate", missing_rate, "-", rand_number_id, "' --min-loci-corr='-1.0' --max-pool-dist='-1.0'"), intern=TRUE)
        print(tmp_fname_out_aldknni_optim)
        fname_out_aldknni_optim = tail(tmp_fname_out_aldknni_optim, n=1)
    }
    duration_aldknni_optim = difftime(Sys.time(), time_ini, units="mins")
    ### LinkImpute's LD-kNN imputation algorithm for unordered genotype data (forcing all data to be diploids)
    fname_for_linkimpute = paste0("LINKIMPUTE_INPUT-maf", maf, "-missing_rate", missing_rate, "-", rand_number_id,".tsv")
    output_for_linkimpute = paste0("LINKIMPUTE_INPUT-maf", maf, "-missing_rate", missing_rate, "-", rand_number_id,"-IMPUTED.tsv")
    fname_out_linkimpute = paste0("LINKIMPUTE_INPUT-maf", maf, "-missing_rate", missing_rate, "-", rand_number_id,"-IMPUTED.csv")
    vcf_for_linkimpute = vcfR::read.vcfR(list_sim_missing$fname_vcf)
    list_genotypes_for_linkimpute = fn_extract_allele_frequencies(vcf_for_linkimpute)
    mat_genotypes_for_linkimpute = t(fn_classify_allele_frequencies(list_genotypes_for_linkimpute$mat_genotypes, ploidy=2)) * 2
    bool_enough_data_to_simulate_10k_missing = sum(!is.na(mat_genotypes_for_linkimpute)) >= 11000
    if (bool_enough_data_to_simulate_10k_missing == TRUE) {
        ### LinkImpute stalls if it cannot mask 10,000 data points for optimising l and k, because the number of non-missing data points is not enough to reach the fixed 10,000 random data points.
        mat_genotypes_for_linkimpute[is.na(mat_genotypes_for_linkimpute)] = -1
        write.table(mat_genotypes_for_linkimpute, file=fname_for_linkimpute, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
        time_ini = Sys.time()
        system(paste0("java -jar ", dir_src, "/res/linkimpute/LinkImpute.jar --verbose -a ", fname_for_linkimpute, " ", output_for_linkimpute)) ### dir_src is defined in perf.R which sources this Rsccipt
        duration_linkimpute = difftime(Sys.time(), time_ini, units="mins")
        mat_linkimputed = read.delim(output_for_linkimpute, header=FALSE)
        rownames(mat_linkimputed) = rownames(mat_genotypes_for_linkimpute)
        colnames(mat_linkimputed) = colnames(mat_genotypes_for_linkimpute)
        mat_linkimputed = t(mat_linkimputed / 2)
        list_loci_names = strsplit(rownames(mat_linkimputed), "_")
        chr = unlist(lapply(list_loci_names, FUN=function(x){(paste(x[1:(length(x)-2)], collapse="_"))}))
        pos = unlist(lapply(list_loci_names, FUN=function(x){as.numeric(x[(length(x)-1)])}))
        allele = unlist(lapply(list_loci_names, FUN=function(x){x[length(x)]}))
        df_linkimputed = data.frame(chr, pos, allele)
        df_linkimputed = cbind(df_linkimputed, mat_linkimputed)
        colnames(df_linkimputed)[1] = "#chr"
        write.table(df_linkimputed, file=fname_out_linkimpute, sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)
    }
    ### Validating imputation
    metrics_mvi = fn_imputation_accuracy(fname_imputed=fname_out_mvi,
        list_sim_missing=list_sim_missing,
        mat_idx_high_conf_data=mat_idx_high_conf_data,
        ploidy=ploidy,
        strict_boundaries=strict_boundaries,
        n_threads=n_threads)
    metrics_aldknni_fixed = fn_imputation_accuracy(fname_imputed=fname_out_aldknni_fixed,
        list_sim_missing=list_sim_missing,
        mat_idx_high_conf_data=mat_idx_high_conf_data,
        ploidy=ploidy,
        strict_boundaries=strict_boundaries,
        n_threads=n_threads)
    metrics_aldknni_optim = fn_imputation_accuracy(fname_imputed=fname_out_aldknni_optim,
        list_sim_missing=list_sim_missing,
        mat_idx_high_conf_data=mat_idx_high_conf_data,
        ploidy=ploidy,
        strict_boundaries=strict_boundaries,
        n_threads=n_threads)
    if (bool_enough_data_to_simulate_10k_missing == TRUE) {
        metrics_linkimpute = fn_imputation_accuracy(fname_imputed=fname_out_linkimpute,
            list_sim_missing=list_sim_missing,
            mat_idx_high_conf_data=mat_idx_high_conf_data,
            ploidy=2,
            strict_boundaries=strict_boundaries,
            n_threads=n_threads)
    } else {
        df_metrics_across_allele_freqs_frequencies = metrics_aldknni_optim$df_metrics_across_allele_freqs_frequencies
        df_metrics_across_allele_freqs_classes = metrics_aldknni_optim$df_metrics_across_allele_freqs_classes
        df_metrics_across_allele_freqs_frequencies[,2:7] = NA
        df_metrics_across_allele_freqs_classes[,2:7] = NA
        metrics_linkimpute = list(
            frac_imputed = 0.0,
            mae_frequencies = NA,
            rmse_frequencies = NA,
            r2_frequencies = NA,
            mae_classes = NA,
            rmse_classes = NA,
            r2_classes = NA,
            concordance_classes = NA,
            df_metrics_across_allele_freqs_frequencies = df_metrics_across_allele_freqs_frequencies,
            df_metrics_across_allele_freqs_classes = df_metrics_across_allele_freqs_classes,
            highConf_mae_frequencies = NA,
            highConf_rmse_frequencies = NA,
            highConf_r2_frequencies = NA,
            highConf_mae_classes = NA,
            highConf_rmse_classes = NA,
            highConf_r2_classes = NA,
            highConf_concordance_classes = NA,
            highConf_df_metrics_across_allele_freqs_frequencies = df_metrics_across_allele_freqs_frequencies,
            highConf_df_metrics_across_allele_freqs_classes = df_metrics_across_allele_freqs_classes
        )
        duration_linkimpute = NA
    }
    ### Merge imputation accuracy metrics into the output data.frame
    string_metric_lists = c("metrics_mvi", "metrics_aldknni_fixed", "metrics_aldknni_optim", "metrics_linkimpute")
    for (m in string_metric_lists) {
        # m = string_metric_lists[1]
        algorithm = gsub("metrics_", "", m)
        vec_basic_metric_names = names(eval(parse(text=m)))
        vec_basic_metric_names = vec_basic_metric_names[grepl("df_", vec_basic_metric_names) == FALSE] ### Excluding the data.frames of MAE across allele frequency bins
        if (algorithm=="linkimpute") {
            ### Forcing all data to be diploids
            df_metrics = data.frame(
                maf,
                n_missing,
                missing_rate,
                n_threads,
                algorithm, 
                2, 
                as.numeric(eval(parse(text=paste0("duration_", algorithm)))),
                matrix(unlist(eval(parse(text=paste0("c(", paste(paste0(m, "$", vec_basic_metric_names), collapse=","), ")")))), nrow=1))
        } else {
            df_metrics = data.frame(
                maf,
                n_missing,
                missing_rate,
                n_threads,
                algorithm, 
                ploidy, 
                as.numeric(eval(parse(text=paste0("duration_", algorithm)))),
                matrix(unlist(eval(parse(text=paste0("c(", paste(paste0(m, "$", vec_basic_metric_names), collapse=","), ")")))), nrow=1))
        }
        colnames(df_metrics) = c(
            "maf",
            "n_missing",
            "missing_rate",
            "n_threads",
            "algorithm",
            "ploidy",
            "duration_mins",
            vec_basic_metric_names)
        ### Insert metrics across allele frequency bins
        df_metrics_across_allele_freqs_frequencies = eval(parse(text=paste0(m, "$df_metrics_across_allele_freqs_frequencies")))
        for (q in df_metrics_across_allele_freqs_frequencies$q) {
            idx = which(df_metrics_across_allele_freqs_frequencies$q == q)
            eval(parse(text=paste0("df_metrics$`mae_", q, "` = df_metrics_across_allele_freqs_frequencies$mae[idx]")))
        }
        ### Bind
        if (m==string_metric_lists[1]) {
            df_out = df_metrics
            colnames(df_out) = names(df_metrics)
        } else {
            df_out = rbind(df_out, df_metrics)
        }
    }
    ### Cleanup
    system(paste0("rm ", list_sim_missing$fname_vcf))
    system(paste0("rm ", gsub("-IMPUTED.csv$", "*", fname_out_mvi)))
    system(paste0("rm ", gsub("-IMPUTED.csv$", "*", fname_out_aldknni_fixed)))
    system(paste0("rm ", gsub("-IMPUTED.csv$", "*", fname_out_aldknni_optim)))
    if (bool_enough_data_to_simulate_10k_missing == TRUE) {
        system(paste0("rm ", gsub("-IMPUTED.csv$", "*", fname_out_linkimpute)))
    }
    ### Output
    return(df_out)
}
