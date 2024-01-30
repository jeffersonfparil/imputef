### Load functions
dir = dirname(sys.frame(1)$ofile)
# dir = "/group/pasture/Jeff/imputef/res"
source(paste0(dir, "/perf_functions.R"))

### Define variable combinations
df_variables = expand.grid(dataset=c("grape", "lucerne", "soybean"), 
                           ploidy=2, 
                           maf=c(0.01, 0.05), 
                           missing_rate=c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
df_variables$ploidy[df_variables$dataset=="lucerne"] = 4
df_variables$ploidy[df_variables$dataset=="soybean"] = 2*42
df_variables = df_variables[order(df_variables$maf, decreasing=TRUE), ]
df_variables = df_variables[order(df_variables$dataset), ]

### Load input variables
args = commandArgs(trailingOnly=TRUE)
# args = c("/group/pasture/Jeff/imputef/misc", "19", "3", "32", "FALSE", "1")
dir_data = args[1]
i = as.numeric(args[2])
n_reps = as.numeric(args[3])
n_threads = as.numeric(args[4])
fname_vcf = paste0(dir_data, "/", df_variables$dataset[i], ".vcf")
ploidy = df_variables$ploidy[i]
maf = df_variables$maf[i]
missing_rate = df_variables$missing_rate[i]
strict_boundaries=FALSE

### Load genotype data
vcf = vcfR::read.vcfR(fname_vcf)
list_genotypes = fn_extract_allele_frequencies(vcf)
mat_genotypes = list_genotypes$mat_genotypes
mat_idx_high_conf_data = list_genotypes$mat_idx_high_conf_data

### Filter by maf >= 0.01
mean_allele_freqs = rowMeans(mat_genotypes, na.rm=TRUE)
idx = which((mean_allele_freqs>=0.01) & ((1-mean_allele_freqs)>=0.01))
vcf = vcf[idx, , ]
mat_genotypes = mat_genotypes[idx, ]
mat_idx_high_conf_data = mat_idx_high_conf_data[idx, ]

### Assess imputation accuracies
for (r in c(1:n_reps)) {
    # r = 1
    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    print(paste0("IMPUTATION PERFORMANCE ASSESSMENT: vcf=", fname_vcf, "; ploidy=", ploidy, " maf=", maf, "; missing_rate=", missing_rate, "; rep=", r))
    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    df_perf = fn_test_imputation(vcf=vcf, 
                                 mat_genotypes=mat_genotypes, 
                                 mat_idx_high_conf_data=mat_idx_high_conf_data, 
                                 ploidy=ploidy, 
                                 maf=maf, 
                                 missing_rate=missing_rate, 
                                 strict_boundaries=strict_boundaries, 
                                 n_threads=n_threads)
    df_perf$rep = rep(r, times=nrow(df_perf))
    if (r==1) {
        TESTS_OUTPUT = df_perf
    } else {
        TESTS_OUTPUT = rbind(TESTS_OUTPUT, df_perf)
    }
    write.table(TESTS_OUTPUT, file=paste0(gsub(".gz$", "", gsub(".vcf$", "", basename(fname_vcf))), "-performance_assessment-maf_", maf, "-missing_rate_", missing_rate, ".csv"), sep=",", quote=FALSE, row.names=FALSE)
}
