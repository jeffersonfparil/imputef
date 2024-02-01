args = commandArgs(trailingOnly=TRUE)
dir_src = args[1]
fname_vcf = args[2]
n_reps = as.numeric(args[3])
ploidy = as.numeric(args[4])
strict_boundaries = as.logical(args[5])
n_threads = as.numeric(args[6])
fname_csv_out = args[7]
# dir_src = "/group/pasture/Jeff/imputef/res"
# fname_vcf = "/group/pasture/Jeff/imputef/misc/lucerne_chromosome1_population2.vcf"
# n_reps = 3
# ploidy = 4
# strict_boundaries = FALSE
# n_threads = 32
# fname_csv_out = "/group/pasture/Jeff/imputef/res/sensitivity_analysis_output.csv"
### Load functions
source(paste0(dir_src, "/perf_functions.R"))
### Load genotype
vcf = vcfR::read.vcfR(fname_vcf)
list_genotypes = fn_extract_allele_frequencies(vcf)
mat_genotypes = list_genotypes$mat_genotypes
### Filter by 5% minor allele frequency:
maf = 0.05
mean_allele_freqs = rowMeans(mat_genotypes, na.rm=TRUE)
idx = which((mean_allele_freqs>=maf) & ((1-mean_allele_freqs)>=maf))
vcf = vcf[idx, , ]
mat_genotypes = mat_genotypes[idx, ]
### Define variable combinations
df_variables = expand.grid(marker_density = seq(from=0.2, to=1.0, by=0.2),
                           sparsity = c(0.01, seq(from=0.1, to=0.6, by=0.1)),
                           min_loci_corr = seq(from=0.0, to=1.0, by=0.1),
                           max_pool_dist = seq(from=0.0, to=1.0, by=0.1),
                           min_l_loci = c(1:5, round(seq(from=6, to=50, length=5))),
                           min_k_neighbours = c(1:5, round(seq(from=6, to=50, length=5))))
idx_sort_rand = sample(1:nrow(df_variables), size=nrow(df_variables), replace=FALSE)
df_variables = df_variables[idx_sort_rand, ]
# df_variables = df_variables[1:10, ]
### Write header
header = c("rep",
           "marker_density",
           "sparsity",
           "min_loci_corr",
           "max_pool_dist",
           "min_l_loci",
           "min_k_neighbours",
           "frac_imputed",
           "mae_frequencies",
           "rmse_frequencies",
           "r2_frequencies",
           "mae_classes",
           "rmse_classes",
           "r2_classes",
           "concordance_classes",
           paste0("mae_q", seq(0,1,by=0.1))
           )
cat(paste0(paste(header, collapse=","), "\n"), file=fname_csv_out)
### Sensitivity analysis across variable combinations
pb = txtProgressBar(min=0, max=n_reps*nrow(df_variables), style=3)
for (i in 1:nrow(df_variables)) {
    for (r in 1:n_reps) {
        list_subsamp = fn_simulate_marker_density_reduction(vcf=vcf, mat_genotypes=mat_genotypes, reduction_rate=df_variables$marker_density[i])
        list_sim_missing = fn_simulate_missing_data(vcf=list_subsamp$vcf,
                                                    mat_genotypes=list_subsamp$mat_genotypes,
                                                    maf=maf,
                                                    missing_rate=df_variables$sparsity[i])
        fname_out = aldknni(fname=list_sim_missing$fname_vcf,
            min_loci_corr=df_variables$min_loci_corr[i],
            max_pool_dist=df_variables$max_pool_dist[i],
            min_l_loci=df_variables$min_l_loci[i],
            min_k_neighbours=df_variables$min_k_neighbours[i],
            optimise_n_steps_min_loci_corr=1,
            optimise_n_steps_max_pool_dist=1,
            optimise_max_l_loci=1,
            optimise_max_k_neighbours=1,
            n_threads=n_threads)
        metrics_out = fn_imputation_accuracy(fname_imputed=fname_out,
            list_sim_missing=list_sim_missing,
            mat_idx_high_conf_data=!is.na(list_subsamp$mat_genotypes),
            ploidy=ploidy,
            strict_boundaries=strict_boundaries,
            n_threads=n_threads)
        vec_metric_names = names(metrics_out)[(!grepl("df_", names(metrics_out))) & (!grepl("highConf_", names(metrics_out)))]
        vec_variables_and_metrics = c(r,
                                      df_variables$marker_density[i],
                                      df_variables$sparsity[i],
                                      df_variables$min_loci_corr[i],
                                      df_variables$max_pool_dist[i],
                                      df_variables$min_l_loci[i],
                                      df_variables$min_k_neighbours[i],
                                      eval(parse(text=paste0("c(", paste(paste0("metrics_out$", vec_metric_names), collapse=","), ")"))),
                                      metrics_out$df_metrics_across_allele_freqs_frequencies$mae
                                    )
        metrics_string = paste0(paste(vec_variables_and_metrics, collapse=","), "\n")
        cat(metrics_string, file=fname_csv_out, append=TRUE)
        system(paste0("rm ", list_sim_missing$fname_vcf))
        system(paste0("rm ", fname_out))
        system(paste0("rm ", gsub("vcf", "vcf-*sync", list_sim_missing$fname_vcf)))
        setTxtProgressBar(pb, r*i)
    }
}
close(pb)
