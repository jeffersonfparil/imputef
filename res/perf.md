# Assessing imputation accuracies

## Variables and datasets

- 3 imputation algorithms: 
    + **ALDKNNI**: adaptive LD-kNN imputation
    + **ALDKNNI_OPTIM_THRESHOLDS**: adaptive LD-kNN imputation with optimisation for the minimum loci correlation and maximum pool distance
    + **ALDKNNI_OPTIM_COUNTS**: adaptive LD-kNN imputation with optimisation for the number of linked loci correlation and k-nearest neighbours
    + **MVI**: mean value imputation
    <!-- + **SAMP**: Luke's LD-kNN imputation algorithm vias sampling from a normal distribution whose parameters depend on the k-nearest neighbours -->
- 2 minor allele frequency thresholds:
    + 0.01
    + 0.05
- 10 sparsity levels (missing rate):
    + 0.01
    + 0.1
    + 0.2
    + 0.3
    + 0.4
    + 0.5
    + 0.6
    + 0.7
    + 0.8
    + 0.9
- 5 datasets (**Note**: Place these data into `imputef/misc`)
    + autotetraploid *Medicago sativa* (2n=4x=32; 2.74 Gb genome; 155 samples x 124,151 biallelic loci; in-house source)
    + pools of diploid *Glycine max* (2n=2x=20; 1.15 Gb genome; 478 pools (each pool comprised of 42 individuals) x 39,636 biallelic loci; source: [http://gong_lab.hzau.edu.cn/Plant_imputeDB/#!/download_soybean](http://gong_lab.hzau.edu.cn/Plant_imputeDB/#!/download_soybean))
    + diploid *Vitis vinifera* (2n=2x=38; 0.5 Gb genome; 77 samples x 8,506 biallelic loci; source: [https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/g3journal/5/11/10.1534_g3.115.021667/5/021667_files1.zip?Expires=1706750617&Signature=yMBQeDumKhnNHIhhUdwsdac~D81t~5RfRi39Bqs4fA8sdE27FVMiyYI7xL8OvLupTqXUim2qC5mgvd5eqby4WCWxCw8x25xnkd6~05gC6puXpHloQSbesTQGrTFios7JeCnXUf306Z~p2vMi0TRgX8qpNTWiwGwwyn2wYAr1tbWIN4EwTQvN8~BgJF31Tj8xJoCVJm2uTpA7~hhsSidJgxVqL4aO20CvwAI1iDcx1gxvienNDS1rYTOruLhwXDif4RGFv8tAb2W5SK3qt4bjgpD6mP8gghv7BWGf0g-arYQywL1fmLCia35qJr7Umxc3LM8iPvWabo5K0sTlRH1oHw__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/g3journal/5/11/10.1534_g3.115.021667/5/021667_files1.zip?Expires=1706750617&Signature=yMBQeDumKhnNHIhhUdwsdac~D81t~5RfRi39Bqs4fA8sdE27FVMiyYI7xL8OvLupTqXUim2qC5mgvd5eqby4WCWxCw8x25xnkd6~05gC6puXpHloQSbesTQGrTFios7JeCnXUf306Z~p2vMi0TRgX8qpNTWiwGwwyn2wYAr1tbWIN4EwTQvN8~BgJF31Tj8xJoCVJm2uTpA7~hhsSidJgxVqL4aO20CvwAI1iDcx1gxvienNDS1rYTOruLhwXDif4RGFv8tAb2W5SK3qt4bjgpD6mP8gghv7BWGf0g-arYQywL1fmLCia35qJr7Umxc3LM8iPvWabo5K0sTlRH1oHw__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA))

- **Genome partitioning:**
    + **By chromosomes or scaffolds**
    + **Rename all chromosomes as chr1**

*prepping_datasets.sh*

```shell
#!/bin/bash
conda activate bcftools
DIR=/group/pasture/Jeff/imputef/misc
cd $DIR
### (2) pools of diploid *Glycine max*
for chr in $(seq 1 20); do wget http://gong_lab.hzau.edu.cn/static/PLimputeDB/download/species/soybean/panel/soybean_impute_Chr${chr}.vcf.gz; done
ls soybean_impute_Chr* > soybean-vcf.txt
for f in $(cat soybean-vcf.txt); do echo $f; gunzip ${f}; bgzip ${f%.gz*}; tabix ${f}; done
bcftools concat $(cat soybean-vcf.txt) -Ov -o soybean-indi_temp.vcf
time Rscript poolify.R \
    soybean-indi_temp.vcf \
    42 \
    100 \
    soybean.vcf
rm soybean_impute_Chr* soybean-vcf.txt soybean-indi_temp.vcf
### (3) diploid grape from the LinkImpute paper with missing data filled in with MVI to be fair with the other datasets
time Rscript \
    ssv2vcf.R \
    LinkImpute.data.grape.num.raw.txt \
    zucchini.vcf
mv LinkImpute.data.grape.num.raw.txt.vcf grape.vcf
```

*poolify.R* - assumes each individual per pool are perfectly equally represented (best-case scenario in the real-world)

```R
args = commandArgs(trailingOnly=TRUE)
# args = c("/group/pasture/Jeff/imputef/misc/soybean-indi_temp.vcf", "42", "100", "/group/pasture/Jeff/imputef/misc/soybean.vcf")
vcf_fname = args[1]
min_pool_size = as.numeric(args[2])
depth = as.numeric(args[3])
out_fname = args[4]
print("###############################################################################################################")
print("Pooling individual diploid genotypes, i.e. each pool is comprised of {min_pool_size} most closely-related samples using k-means clustering using 1000 evenly indexes loci")
print("Note: assumes each individual per pool are perfectly equally represented (best-case scenario in the real-world)")
vcf = vcfR::read.vcfR(vcf_fname)
vec_loci_names = paste(vcfR::getCHROM(vcf), vcfR::getPOS(vcf), vcfR::getREF(vcf), sep="_")
vec_sample_names = colnames(vcf@gt)[-1]
### Extract biallelic diploid allele frequencies
G = matrix(NA, nrow=length(vec_loci_names), ncol=length(vec_sample_names))
GT = vcfR::extract.gt(vcf, element="GT")
G[(GT == "0/0") | (GT == "0|0")] = 1.0
G[(GT == "1/1") | (GT == "1|1")] = 0.0
G[(GT == "0/1") | (GT == "0|1") | (GT == "1|0")] = 0.5
rm(GT)
gc()
rownames(G) = vec_loci_names
colnames(G) = vec_sample_names
### Create n pools of the most related individuals
print("Clustering. This may take a while.")
PCA = prcomp(G, rank=100)
C = PCA$rotation
p = length(vec_loci_names)
# C = t(G[seq(from=1, to=p, length=1000), ])
clustering = kmeans(x=C, centers=200, iter.max=20)
vec_clusters = table(clustering$cluster)
vec_clusters = vec_clusters[vec_clusters >= min_pool_size]
n = length(vec_clusters)
GT = matrix("AD", nrow=p, ncol=n+1)
pb = txtProgressBar(min=0, max=n, initial=0, style=3)
for (i in 1:n) {
    idx = which(clustering$cluster == i)
    q = rowMeans(G[, idx])
    for (j in 1:p) {
        # j = 1
        ref = rbinom(n=1, size=depth, prob=q[j])
        alt = rbinom(n=1, size=depth, prob=(1-q[j]))
        GT[j, (i+1)] = paste0(ref, ",", alt) ### skipping the first column containing the FORMAT field "AD"
    }
    setTxtProgressBar(pb, i)
}
close(pb)
colnames(GT) = c("FORMAT", paste0("Pool-", c(1:n)))
### Create the output vcfR object
vcf_out = vcf
str(vcf_out)
str(vcf_out@gt)
vcf_out@gt = GT
fname_vcf_gz = paste0(out_fname, ".gz")
vcfR::write.vcf(vcf_out, file=fname_vcf_gz)
system(paste0("gunzip -f ", fname_vcf_gz))
print(paste0("Output: ", out_fname))
```

*ssv2vcf.R* - convert space-delimited genotype data from LinkImpute paper

```R
args = commandArgs(trailingOnly = TRUE)
# args = c("/group/pasture/Jeff/imputef/misc/LinkImpute.data.apple.num.raw.txt", "/group/pasture/Jeff/imputef/misc/zucchini.vcf")
# args = c("/group/pasture/Jeff/imputef/misc/LinkImpute.data.grape.num.raw.txt", "/group/pasture/Jeff/imputef/misc/zucchini.vcf")
fname_geno_txt = args[1] 
fname_geno_dummy_vcf = args[2]
dat = read.table(fname_geno_txt, header=TRUE, sep=" ")
vcf = vcfR::read.vcfR(fname_geno_dummy_vcf)
idx_col_start = 7
n = nrow(dat)
p = ncol(dat) - (idx_col_start-1)
vec_loci_names = colnames(dat)[idx_col_start:ncol(dat)]
vec_pool_names = dat$IID
G = dat[, idx_col_start:ncol(dat)]
# ### Missing data assessment
# idx_row_without_missing = apply(G, MARGIN=1, FUN=function(x){sum(is.na(x))==0})
# idx_col_without_missing = apply(G, MARGIN=2, FUN=function(x){sum(is.na(x))==0})
# sum(idx_row_without_missing)
# sum(idx_col_without_missing)
# sum(is.na(G)) / prod(dim(G))
# ### Fill missing with mean (because LinkImpute fails to finish running at MAF=0.25 for reasons I do not know)
# ploidy = 2
# pb = txtProgressBar(min=0, max=p, style=3)
# for (j in 1:p) {
#     # j = 1
#     idx_missing = which(is.na(G[, j]))
#     G[idx_missing, j] = round(mean(G[, j], na.rm=TRUE) * ploidy)
#     setTxtProgressBar(pb, j)
# }
# close(pb)
### Create meta, fit and gt fields of the vcfR object
mat_loci_ids = matrix(gsub("^X", "chr_", unlist(strsplit(unlist(strsplit(vec_loci_names, "[.]")), "_"))), ncol=3, byrow=TRUE)
### Randomly choose the alternative alleles as the genotype data (*.raw file) do not have that information.
vec_alt = c()
vec_alleles = c("A", "T", "C", "G")
for (i in 1:nrow(mat_loci_ids)) {
    vec_alt = c(vec_alt, sample(vec_alleles[vec_alleles != mat_loci_ids[i,3]], size=1))
}
vec_chr = unique(mat_loci_ids[,1])
META = c("##fileformat=VCFv", paste0("##", vec_chr), "##Extracted from text file.")
FIX = cbind(mat_loci_ids[,1], mat_loci_ids[,2], vec_loci_names, mat_loci_ids[,3], vec_alt, rep(NA, each=p), rep("PASS", each=p), rep(NA, each=p))
colnames(FIX) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
GT = matrix("", nrow=p, ncol=(n+1))
GT[,1] = "GT"
colnames(GT) = c("FORMAT", vec_pool_names)
pb = txtProgressBar(min=0, max=n, style=3)
for (i in 1:n) {
    for (j in 1:p) {
        ### We're assuming that the counts in G are for the reference alleles
        if (is.na(G[i, j])) {
            GT[j, (i+1)] = "./."
        } else {
            if (G[i, j] == 0) {
                GT[j, (i+1)] = "1/1"
            } else if (G[i, j] == 1) {
                GT[j, (i+1)] = "0/1"
            } else if (G[i, j] == 2) {
                GT[j, (i+1)] = "0/0"
            } else {
                GT[j, (i+1)] = "./."
            }
        }
    }
    setTxtProgressBar(pb, i)
}
close(pb)
### Create the new vcfR object
vcf_new = vcf
vcf_new@meta = META
vcf_new@fix = FIX
vcf_new@gt = GT
print(vcf)
print(vcf_new)
### Save and unzip the object
vcfR::write.vcf(vcf_new, file=paste0(fname_geno_txt, ".vcf.gz"))
system(paste0("gunzip -f ", fname_geno_txt, ".vcf.gz"))
vcf_new_loaded = vcfR::read.vcfR(paste0(fname_geno_txt, ".vcf"))
print(vcf_new_loaded)
```


## Assess the genetic relationships between samples per dataset

```R
dir = "/group/pasture/Jeff/imputef/res/"
vec_fnames = paste0("/group/pasture/Jeff/imputef/misc/", c("grape.vcf", "lucerne.vcf", "soybean.vcf"))
setwd(dir)
### Allele frequency extraction
fn_extract_allele_frequencies = function(vcf) {
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
    } else {
        mat_ref_counts = vcfR::masplit(mat_allele_counts, delim=',', record=1, sort=0)
        mat_alt_counts = vcfR::masplit(mat_allele_counts, delim=',', record=2, sort=0)
        ### Set missing allele counts to 0, if the other allele is non-missing and non-zero
        idx_for_ref = which(!is.na(mat_alt_counts) & (mat_alt_counts != 0) & is.na(mat_ref_counts))
        idx_for_alt = which(!is.na(mat_ref_counts) & (mat_ref_counts != 0) & is.na(mat_alt_counts))
        mat_ref_counts[idx_for_ref] = 0
        mat_alt_counts[idx_for_alt] = 0
        ### Calculate reference allele frequencies
        mat_genotypes = mat_ref_counts / (mat_ref_counts + mat_alt_counts)
    }
    ### Label the loci and pools
    rownames(mat_genotypes) = vec_loci_names
    colnames(mat_genotypes) = vec_pool_names
    ### Output
    return(mat_genotypes)
}
### Assess genotype relationships per dataset
for (i in 1:length(vec_fnames)) {
    # i = 1
    fname = vec_fnames[i]
    fname_png = gsub(".vcf", ".png", basename(fname))
    print("=====================================")
    print(fname)
    vcf = vcfR::read.vcfR(fname, verbose=FALSE)
    # print(vcf)
    mat_geno = fn_extract_allele_frequencies(vcf)
    # str(mat_geno)
    C = cor(mat_geno)
    png(fname_png, width=2000, height=2000)
    heatmap(C, scale="none", main=gsub(".vcf", "", basename(fname)))
    dev.off()
}
```

![estimated relationships between samples](../res/lucerne.png)

![estimated relationships between samples](../res/soybean.png)

![estimated relationships between samples](../res/grape.png)


## Prepare LinkImpute for testing against diploid imputation

*prepare_linkimpute.sh*

```shells
#!/bin/bash
DIR=/group/pasture/Jeff/imputef/res
cd $DIR
mkdir linkimpute/
cd linkimpute/
wget http://www.cultivatingdiversity.org/uploads/5/8/2/9/58294859/20210706_linkimpute.tar.gz
tar -xvzf 20210706_linkimpute.tar.gz
rm 20210706_linkimpute.tar.gz
java -jar LinkImpute.jar -h
```


## Metrics

- Concordance: $c = {{1 \over n} \Sigma_{i=1}^{n} p}$, where: $p=
\begin{cases}
0 \text{ if } \hat g \ne g_{true}\\
1 \text{ if } \hat g = g_{true}
\end{cases}
$.
This is used for genotype classes, i.e., binned allele frequencies: $g = {{1 \over {ploidy}} round(q*ploidy)}$, here $q = P(allele)$. Note that there is alternative way of defining these genotype classes with strict boundaries, i.e., homozygotes have fixed allele frequencies.
- Mean absolute error: $mae = {{1 \over n} \Sigma_{i=1}^{n}|\hat q - q_{true}|}$.
- Coefficient of determination: $R^2 = { 1 - {{\Sigma_{}^{}(\hat q - q_{true})^2} \over {\Sigma_{}^{}(\hat q_{true} - \bar q_{true})^2}} }$


## Execution

```shell
### Submit array jobs for each dataset using specific memory and time limits
DIR=/group/pasture/Jeff/imputef/res
cd $DIR
for DATASET in "grape" "lucerne" "soybean"
do
    if [ $DATASET == "grape" ]
    then
        INI=1
        FIN=20
        sed 's/--job-name="imputef"/--job-name="grapeImp"/g' perf.slurm | \
            sed 's/--mem=250G/--mem=100G/g' | \
            sed 's/--time=14-0:0:00/--time=0-0:30:00/g' > perf_${DATASET}.slurm
    elif [ $DATASET == "lucerne" ]
    then
        INI=21
        FIN=40
        sed 's/--job-name="imputef"/--job-name="lucerImp"/g' perf.slurm > perf_${DATASET}.slurm
    else
        INI=41
        FIN=60
        sed 's/--job-name="imputef"/--job-name="soyImp"/g' perf.slurm | \
            sed 's/--mem=250G/--mem=200G/g' | \
            sed 's/--time=14-0:0:00/--time=7-0:0:00/g' > perf_${DATASET}.slurm
    fi
    echo ${DATASET}: ${INI}-${FIN}
    sbatch --array=${INI}-${FIN} perf_${DATASET}.slurm
done

### Monitor the jobs
conda activate rustenv
DIR=/group/pasture/Jeff/imputef/res
cd $DIR
squeue -u jp3h | sort
SLURMOUT_GRAPE=slurm-24317042_*.out
SLURMOUT_LUCERNE=slurm-??????_*.out
SLURMOUT_SOYBEAN=slurm-24317065_*.out
grep -n -i "err" slurm-24317065*_*.out | grep -v "mean absolute"
tail slurm-24317065*_*.out
ls -lh *-performance_assessment-maf_*missing_rate_*.csv
ls -lhtr
time Rscript perf_plot.R ${DIR}
# scancel -u jp3h
# rm slurm-* soybean-*.csv lucerne-*.csv zucchini-*.csv apple-*.csv grape-*.csv LINKIMPUTE* AOPT*-maf0.* AFIXED*-maf0.* MVI-maf0.* ploidy_vcf-* SIMULATED_MISSING-0.*

### After all jobs have finished, move the output and plot:
mkdir output
mv *-performance_assessment-maf_*-missing_rate_*.csv output/
time Rscript summary_plot.R \
    ${DIR}/output
```

## Miscellaneous: sensitivity analysis

Using a single small chromosome from the Lucerne dataset we will explore the entire parameter spaces across sparsity and marker density. We will test all possible combinations of the 4 parameters (minimum loci correlation, maximum pool distance, minimum number of linked loci, and minimum number of k-nearest neighbours).

We'll start with subsetting `lucerne.vcf` so that we only arbitrarily include *chromome 1* loci (chromsome ID: chr1.4) and *population 2*:

```shell
DIR=/group/pasture/Jeff/imputef/misc
VCF=${DIR}/lucerne.vcf
cd $DIR
echo "Looking at the number of loci covered per chromosome:"
for i in $(seq 1 8)
do
N=$(grep "chr${i}.4" ${VCF} | wc -l)
echo CHR${i}: ${N}
done
echo "Extracting chromosome 1 loci:"
head -n6 $VCF > ${DIR}/lucerne_chromosome1_population2.vcf.tmp ### header lines
grep -m1 "#CHR" $VCF >> ${DIR}/lucerne_chromosome1_population2.vcf.tmp ### column labels including sample names
grep "chr1.4" ${VCF} | grep -v "^##contig" >> ${DIR}/lucerne_chromosome1_population2.vcf.tmp
echo "Extracting populations 2 samples:"
for i in $(seq 1 7)
do
N=$(grep -m1 "^#CHR" ${DIR}/lucerne_chromosome1_population2.vcf.tmp | sed -z "s/\t/\n/g" | grep "DB-MS-31-22-00${i}" | wc -l)
echo POP_${i}: ${N}
done
IDX=$(echo 1-9,$(grep -m1 "^#CHR" ${DIR}/lucerne_chromosome1_population2.vcf.tmp | sed -z "s/\t/\n/g" | grep -n "DB-MS-31-22-002" | cut -d':' -f1 | sed -z 's/\n/,/g' | sed 's/,$//g'))
cut -f${IDX} ${DIR}/lucerne_chromosome1_population2.vcf.tmp > ${DIR}/lucerne_chromosome1_population2.vcf
rm ${DIR}/lucerne_chromosome1_population2.vcf.tmp
wc -l ${DIR}/lucerne_chromosome1_population2.vcf
bat -l tsv --wrap never ${DIR}/lucerne_chromosome1_population2.vcf
```

Now, we will assess imputation accuracy across various combinations of the 4 parameters, as well as across 10 sparsity levels, and 10 marker density levels:

```shell
DIR=/group/pasture/Jeff/imputef/res
cd $DIR
sbatch perf_misc.slurm

### Monitoring:
conda activate rustenv
DIR=/group/pasture/Jeff/imputef/res
cd $DIR
squeue -u jp3h | sort
SLURMOUT=slurm-24312170.out
grep -n -i "err" ${SLURMOUT} | grep -v "mean absolute"
wc -l ${DIR}/sensitivity_analysis_output.csv
bat --wrap never ${DIR}/sensitivity_analysis_output.csv
```

### Questions we wish to anwer:

1. Is the 4D parameter space (*min_loci_corr*, *max_pool_dist*, *min_l_loci*, and *min_k_neighbours*) actually smooth and unimodal?
2. How does *marker_density* affect imputation accuracy?
3. How does *marker_density* affect the optimum parameters?
4. Is there a single combination of parameters which yields reasonably high imputation accuracy across marker density and sparsity levels?

```R
dir = "/group/pasture/Jeff/imputef/res"
setwd(dir)
df = read.csv("sensitivity_analysis_output.csv")
str(df)

### Evaluate the data
summary(df)
sum(is.na(df))
sum(is.na(df[, grepl("mae_q", colnames(df))]))
sum(is.na(df[, !grepl("mae_q", colnames(df))]))

### Linear modelling
mod = lm(mae_frequencies ~ marker_density + sparsity + min_loci_corr + max_pool_dist + min_l_loci + min_k_neighbours, data=df)
summary(mod)
mod_complete = lm(mae_frequencies ~ marker_density * sparsity * min_loci_corr * max_pool_dist * min_l_loci * min_k_neighbours, data=df)
summary(mod_complete)
mod_null = lm(mae_frequencies ~ 1, data=df)
summary(mod_null)

anova(mod_null, mod, mod_complete)

txtplot::txtplot(x=df$marker_density, y=df$mae_frequencies, xlab="marker_density", ylab="mae_frequencies")
txtplot::txtplot(x=df$sparsity, y=df$mae_frequencies, xlab="sparsity", ylab="mae_frequencies")
txtplot::txtplot(x=df$min_loci_corr, y=df$mae_frequencies, xlab="min_loci_corr", ylab="mae_frequencies")

### Random forest
data = data.frame(mae_frequencies=df$mae_frequencies, 
                  marker_density=df$marker_density, 
                  sparsity=df$sparsity, 
                  min_loci_corr=df$min_loci_corr, 
                  max_pool_dist=df$max_pool_dist, 
                  min_l_loci=df$min_l_loci, 
                  min_k_neighbours=df$min_k_neighbours)
# data$marker_density_X_sparsity = as.factor(paste0(data$marker_density, "_X_", data$sparsity))
data$marker_density_X_sparsity = data$marker_density * data$sparsity

data$min_loci_corr_X_max_pool_dist = data$min_loci_corr * data$max_pool_dist
data$min_l_loci_X_min_k_neighbours = data$min_l_loci * data$min_k_neighbours

data$min_loci_corr_X_min_k_neighbours = data$min_loci_corr * data$min_k_neighbours
data$min_l_loci_X_max_pool_dist = data$min_l_loci * data$max_pool_dist

data$min_loci_corr_X_min_l_loci = data$min_loci_corr * data$min_l_loci
data$max_pool_dist_X_min_k_neighbours = data$max_pool_dist * data$min_k_neighbours

data$min_loci_corr_X_max_pool_dist_X_min_l_loci_X_min_k_neighbours = data$min_loci_corr * data$max_pool_dist * data$min_l_loci * data$min_k_neighbours

# rf = partykit::cforest(mae_frequencies ~ ., data=data, trace=TRUE, cores=32)
# summary(rf)
# VARIMP = partykit::varimp(rf, conditional=TRUE, cores=2) ### uses mclapply - hence will be cloning stuff and eating a lot of RAM - therefore reduces the number of cores as required

# summary(VARIMP)
# str(VARIMP)
# print(VARIMP)
# txtplot::txtplot(VARIMP)

# ### TESTING partykit::cforest
# # tmp_testing_cforest.R
# # tmp_testing_cforest.slurm
# # DIR=/group/pasture/Jeff/imputef/res
# # cd $DIR
# # cat slurm-24309707.out

# # rf = randomForest::randomForest(mae_frequencies ~ ., data=train, proximity=TRUE)
# # print(rf)
# # p1 = predict(rf, train)
# # cor(train$mae_frequencies, p1)
# # txtplot::txtplot(x=train$mae_frequencies, y=p1)
# # p2 = predict(rf, test)
# # cor(test$mae_frequencies, p2)
# # txtplot::txtplot(x=test$mae_frequencies, y=p2)

random_forest_output = randomForest::randomForest(mae_frequencies ~ ., data=data, proximity=TRUE)
cor(data$mae_frequencies, predict(random_forest_output, data))
str(random_forest_output)

IMPORTANCE = random_forest_output$importance[order(random_forest_output$importance, decreasing=TRUE), , drop=FALSE]
print(IMPORTANCE)

print("The choice of the nearest neighbours appears to be the most important variables, i.e. max_pool_dist and min_k_neighbours.")
print("This is followed by the dataset marker density and sparsity combination.")
print("Then by the correlation-by-distance threshold combinations and minimum l-by-k combinations.")

print("Question 2: How does marker density affect imputation accuracy?")
print("Answer 2: Marker density per se has the least effect on imputation accuracy. 
       However, the combination of marker density and sparsity has a magnitude greater effect on imputation accuracy than marker density alone.")

### For each combination of marker density and sparsity levels ...
###     assess imputation accuracy modality across the 4 parameter spaces, and ...
###     identify the optimum parameter combination.
vec_dataset_combinations = sort(unique(paste0(df$marker_density, "_X_", df$sparsity)))
for (dataset in vec_dataset_combinations) {
    # dataset = vec_dataset_combinations[1]
    print("#######################################################")
    print(dataset)
    marker_density = as.numeric(unlist(strsplit(dataset, "_X_"))[1])
    sparsity = as.numeric(unlist(strsplit(dataset, "_X_"))[2])
    subdf = df[(df$marker_density==marker_density) & (df$sparsity==sparsity), ]
    # txtplot::txtplot(x=df$min_loci_corr, y=df$mae_frequencies, xlab="min_loci_corr", ylab="mae")
    # txtplot::txtplot(x=df$max_pool_dist, y=df$mae_frequencies, xlab="max_pool_dist", ylab="mae")
    # txtplot::txtplot(x=df$min_l_loci, y=df$mae_frequencies, xlab="min_l_loci", ylab="mae")
    # txtplot::txtplot(x=df$min_k_neighbours, y=df$mae_frequencies, xlab="min_k_neighbours", ylab="mae")

    ### Describe the parameter spaces
    ### Are the trends in MAE of a parameter the same across levels of all other parameters? Not likely!
    vec_min_loci_corr = sort(unique(subdf$min_loci_corr))
    vec_max_pool_dist = sort(unique(subdf$max_pool_dist))
    vec_min_l_loci = sort(unique(subdf$min_l_loci))
    vec_min_k_neighbours = sort(unique(subdf$min_k_neighbours))
    for (min_loci_corr in vec_min_loci_corr) {
        for (max_pool_dist in vec_max_pool_dist) {
            for (min_l_loci in vec_min_l_loci) {
                # min_loci_corr=0.7; max_pool_dist=0.8; min_l_loci=3
                idx = which((subdf$min_loci_corr==min_loci_corr) & (subdf$max_pool_dist==max_pool_dist) & (subdf$min_l_loci==min_l_loci))
                if (length(idx)==0) {
                    next
                }
                dat = subdf[idx, ]
                # txtplot::txtplot(x=dat$min_k_neighbours, y=dat$mae_frequencies, xlab="min_k_neighbours", ylab="MAE")
                agg = aggregate(mae_frequencies ~ min_k_neighbours, data=dat, FUN=mean)
            }
        }
    }

    ### Despite these, maybe the current output is too granular, and maybe as we get more info it all smooths out and the parameter spaces become more or less smooth and unimodal?
    

    ### Find the very best parameter combination at the current dataset
    agg = aggregate(mae_frequencies ~ min_loci_corr + max_pool_dist + min_l_loci + min_k_neighbours, data=subdf, FUN=mean)
    agg = agg[order(agg$mae_frequencies, decreasing=FALSE), ]
    optim = data.frame(marker_density=marker_density, sparsity=sparsity, agg[1, ])
    if (dataset == vec_dataset_combinations[1]) {
        df_best_params = optim
    } else {
        df_best_params = rbind(df_best_params, optim)
    }
    agg_min_loci_corr = aggregate(mae_frequencies ~ min_loci_corr, data=subdf, FUN=mean)
    txtplot::txtplot(x=agg_min_loci_corr$min_loci_corr, y=agg_min_loci_corr$mae_frequencies, xlab="min_loci_corr", ylab="mae")
    agg_max_pool_dist = aggregate(mae_frequencies ~ max_pool_dist, data=subdf, FUN=mean)
    txtplot::txtplot(x=agg_max_pool_dist$max_pool_dist, y=agg_max_pool_dist$mae_frequencies, xlab="max_pool_dist", ylab="mae")
    agg_min_l_loci = aggregate(mae_frequencies ~ min_l_loci, data=subdf, FUN=mean)
    txtplot::txtplot(x=agg_min_l_loci$min_l_loci, y=agg_min_l_loci$mae_frequencies, xlab="min_l_loci", ylab="mae")
    agg_min_k_neighbours = aggregate(mae_frequencies ~ min_k_neighbours, data=subdf, FUN=mean)
    txtplot::txtplot(x=agg_min_k_neighbours$min_k_neighbours, y=agg_min_k_neighbours$mae_frequencies, xlab="min_k_neighbours", ylab="mae")
}

print("Question 1: Is MAE smooth and/or unimodal per parameter space?")
print("Answer 1: No. The parameter spaces are neither smooth nor unimodal. The parameter spaces are complex with some local minima far from the global minimum.")

print(df_best_params)

txtplot::txtplot(x=df_best_params$marker_density, y=df_best_params$mae_frequencies, xlab="marker_density", ylab="mae")
agg_marker_density = aggregate(mae_frequencies ~ marker_density, data=df_best_params, FUN=mean)
print(agg_marker_density[order(agg_marker_density$mae_frequencies, decreasing=FALSE), ])

print("Question 3: How does marker density affect the optimum parameter combinations?")
print("Answer 3: On average using the best imputation accuracy per marker-density-by-sparsity dataset, higher marker density improves imputation accuracy.")

txtplot::txtdensity(df_best_params$min_loci_corr, xlab="min_loci_corr")
txtplot::txtdensity(df_best_params$max_pool_dist, xlab="max_pool_dist")
txtplot::txtdensity(df_best_params$min_l_loci, xlab="min_l_loci")
txtplot::txtdensity(df_best_params$min_k_neighbours, xlab="min_k_neighbours")

mean_optim_params = as.data.frame(t(apply(df_best_params, MAR=2, FUN=mean)))
print(mean_optim_params)
pred_optim_min_loci_corr = round(mean_optim_params$min_loci_corr, 1)
pred_optim_max_pool_dist = round(mean_optim_params$max_pool_dist, 1)

vec_min_l_loci = sort(unique(df$min_l_loci))
delta_min_l_loci = abs(vec_min_l_loci - round(mean_optim_params$min_l_loci))
pred_optim_min_l_loci = vec_min_l_loci[which(delta_min_l_loci == min(delta_min_l_loci))]

vec_min_k_neighbours = sort(unique(df$min_k_neighbours))
delta_min_k_neighbours = abs(vec_min_k_neighbours - round(mean_optim_params$min_k_neighbours))
pred_optim_min_k_neighbours = vec_min_k_neighbours[which(delta_min_k_neighbours == min(delta_min_k_neighbours))]

idx = which((df$min_loci_corr == pred_optim_min_loci_corr) & (df$max_pool_dist == pred_optim_max_pool_dist) & (df$min_l_loci == pred_optim_min_l_loci) & (df$min_k_neighbours == pred_optim_min_k_neighbours))

df[idx, ]

print("Question 4: Is there a single combination of parameters which performs reasonably well across datasets?")
print("Answer 4: We are not certain yet as we are only capturing a single dataset for which the predicted optimum parameter combination is present. 
       However, given the complexity of the parameter spaces, this is unlikely. Again, blame the no-free-lunch theorem.
       But wait! Some unimodality in min_loci_corr is starting to show centered on 0.5, and bimodalities in the rest, i.e.,
       max_pool_dist near both extremes, min_l_loci at <10 and ~17, and min_k_neighbours at ~5 and ~17. These are very curiosome indeed!")

```

### Preliminary results (20240201) --> being addressed with the most recent changes/merge which should improve optimisation where we move forwards then backwards across the parameter spaces (2024/02/02)

*Answers to [question](#questions-we-wish-to-anwer), possible consequences and follow-up questions:*

1. No. The parameter spaces are neither smooth nor unimodal (*Figure 1*). The parameter spaces are complex with some local minima far from the global minimum. This means we have to improve our optimisation algorithm in [`optim.rs`](../src/rust/src/optim.rs). Also, note that the choice of the nearest neighbours (`max_pool_dist_X_min_k_neighbours`) appears to be the most important set of variables in terms imputation accuracy (*Table 1*). Will this mean that we can set `min_loci_corr` to some low value and likewise `min_l_loci` to an arbitrary value that is not too high as low minimum correlation should capture sufficiently numerous loci for genetic distance estimation; and then just optimise for `min_k_neighbours` and `max_pool_dist`?
2. Marker density per se has the least effect on imputation accuracy (*Table 1*). However, the combination of marker density and sparsity has a magnitude greater effect on imputation accuracy than marker density alone. But what does this mean? Are they antagonistic, i.e. low sparsity but high density means better accuracy and vice-versa as expected? What is the accuracy in the middle, i.e. medium sparsity and medium density? How about when both are high or both are low?
3. On average using the best imputation accuracy per marker-density-by-sparsity dataset, higher marker density improves imputation accuracy (*Table 2*).
4. We are not certain yet as we are only capturing a single dataset for which the predicted optimum parameter combination is present. However, given the complexity of the parameter spaces, this is unlikely. Again, blame the no-free-lunch theorem. For the current findings, please see *Table 3* for the mean optimum parameter values across datasets (marker-density-by-sparsity combinations).  But wait, look at *Figure 2* and notice some unimodality in `min_loci_corr` centered on 0.5, and bimodalities in the rest, i.e., `max_pool_dist` near both extremes, `min_l_loci` at <10 and ~17, and `min_k_neighbours` at ~5 and ~17. These are very curiosome indeed!


*Table 1*. Variable importance based on random forest regression

| Variable                                                        | Importance  |
| :-------------------------------------------------------------- | ----------: |
| max_pool_dist_X_min_k_neighbours                                | 0.058613063 |
| min_k_neighbours                                                | 0.038643309 |
| max_pool_dist                                                   | 0.035782306 |
| min_loci_corr_X_min_k_neighbours                                | 0.023709804 |
| marker_density_X_sparsity                                       | 0.016431124 |
| min_loci_corr_X_max_pool_dist                                   | 0.014462863 |
| min_l_loci_X_min_k_neighbours                                   | 0.007870009 |
| min_loci_corr                                                   | 0.007014149 |
| min_loci_corr_X_min_l_loci                                      | 0.006842891 |
| min_l_loci_X_max_pool_dist                                      | 0.006793837 |
| sparsity                                                        | 0.004732896 |
| min_loci_corr_X_max_pool_dist_X_min_l_loci_X_min_k_neighbours   | 0.004610137 |
| min_l_loci                                                      | 0.004595450 |
| marker_density                                                  | 0.001893720 |


*Table 2*. Effect of marker density on imputation accuracy (measured by mean absolute error, MAE)

| marker_density | mae_frequencies |
| :------------: | --------------: |
|            1.0 |      0.06718490 |
|            0.8 |      0.06766566 |
|            0.2 |      0.06786764 |
|            0.6 |      0.06796051 |
|            0.4 |      0.06828012 |


*Table 3*. Mean optimum parameter values across datasets, i.e. marker-density-by-sparsity combinations

| Variable          | Optimum value |
| :---------------- | ------------: |
| marker_density    |           0.4 |
| sparsity          |           0.3 |
| min_loci_corr     |           0.4 |
| max_pool_dist     |           0.4 |
| min_l_loci        |            17 |
| min_k_neighbours  |            17 |


*Figure 1*. Representative parameter spaces across datasets (marker-density-by-sparsity combinations)

```
         +-+---------+---------+---------+---------+--------+--+        |               +-+---------+---------+---------+---------+---------+--+
   0.071 + *                                                   +        |         0.074 +                                                   *  +
         |                                                     |        |               | *                                                    |
         |                                                     |        |               |                                                      |
  0.0705 +                                                     +        |         0.073 +                                                      +
         |           *                             *           |        |               |                                              *       |
m        |                                                  *  |        |       m       |                                                      |
a   0.07 +                               *                     +        |       a       |                                                      |
e        |                *         *                          |        |       e 0.072 +                                         *            +
         |      *                                       *      |        |               |           *    *                                     |
  0.0695 +                                    *                +        |               |      *                                               |
         |                                                     |        |         0.071 +                                    *                 +
         |                     *                               |        |               |                     *    *    *                      |
   0.069 +-+---------+---------+---------+---------+--------+--+        |               +-+---------+---------+---------+---------+---------+--+
           0        0.2       0.4       0.6       0.8       1           |                 0        0.2       0.4       0.6       0.8        1   
                              min_loci_corr                             |                                     min_loci_corr                     
        +-+---------+---------+---------+---------+---------+--+        |               +-+---------+---------+---------+---------+---------+--+
        | *                                                 *  |        |               |                                                   *  |
  0.073 +                                                      +        |         0.076 +                                                      +
        |                                                      |        |               |                                                      |
  0.072 +                                                      +        |         0.075 + *                                                    +
        |                                                      |        |               |                                                      |
m       |                                                      |        |       m 0.074 +                                                      +
a 0.071 +                                                      +        |       a       |                                                      |
e       |                                                      |        |       e 0.073 +                                                      +
   0.07 +                                                      +        |               |                                                      |
        |                *    *    *    *    *    *    *       |        |         0.072 +                                                      +
  0.069 +           *                                          +        |               |                *    *    *    *         *            |
        |      *                                               |        |         0.071 +      *    *                        *         *       +
        +-+---------+---------+---------+---------+---------+--+        |               +-+---------+---------+---------+---------+---------+--+
          0        0.2       0.4       0.6       0.8        1           |                 0        0.2       0.4       0.6       0.8        1   
                              max_pool_dist                             |                                     max_pool_dist                     
         ++---------+---------+---------+---------+---------+--+        |                ++---------+---------+---------+---------+---------+--+
   0.071 +    *                                                +        |                |  *                                                  |
         |                                                     |        |                |                            *                        |
  0.0705 +                                                     +        |         0.0725 +                                                     +
         |   *                                                 |        |                |                                                     |
         |     *                      *                     *  |        |                |                                                     |
m   0.07 +                                                     +        |       m        |      *                                              |
a        |      *          *                     *             |        |       a  0.072 + *               *                                *  +
e 0.0695 +                                                     +        |       e        |    *                                                |
         |  *                                                  |        |                |   *                                                 |
         |                                                     |        |         0.0715 +                                                     +
   0.069 +                                                     +        |                |                                                     |
         | *                                                   |        |                |     *                                 *             |
  0.0685 ++---------+---------+---------+---------+---------+--+        |                ++---------+---------+---------+---------+---------+--+
          0        10        20        30        40        50           |                 0        10        20        30        40        50   
                               min_l_loci                               |                                      min_l_loci                       
  0.0725 ++---------+---------+---------+---------+---------+--+        |               ++----------+---------+---------+---------+---------+--+
         | *                                                   |        |         0.075 + *                                                    +
   0.072 +                                                     +        |               |                                                      |
         |                                                     |        |               |                                                      |
  0.0715 +                                                     +        |         0.074 +                                                      +
         |                                                     |        |               |                                                      |
m  0.071 +                                                     +        |       m 0.073 +                                                      +
a        |  *                                                  |        |       a       |                                                      |
e 0.0705 +                                                     +        |       e       |   *                                                  |
         |                                                     |        |         0.072 +    **                                                +
    0.07 +     *                                            *  +        |               |      *                                               |
  0.0695 +   ** *                                              +        |         0.071 +       *                     *                        +
         |                 *          *          *             |        |               |                  *                     *          *  |
         ++---------+---------+---------+---------+---------+--+        |               ++----------+---------+---------+---------+---------+--+
          0        10        20        30        40        50           |                0         10        20        30        40        50   
                            min_k_neighbours                            |                                   min_k_neighbours
```

*Figure 2*. Distribution of the best parameters across datasets (marker-density-by-sparsity combinations)

```
    +--+---------+----------+----------+---------+----------+--+
    |                          *******                         |
    |                       ***       **                       |
1.5 +                     ***          ***                     +
    |                    **              **                    |
    |                  **                 ***                  |
    |                 **                    **                 |
  1 +               ***                      ***               +
    |              **                          ***             |
    |             **                             ***           |
0.5 +            **                                *****       +
    |          ***                                     ******  |
    |  *********                                               |
    +--+---------+----------+----------+---------+----------+--+
       0        0.2        0.4        0.6       0.8         1   
                            min_loci_corr                       
    +--+---------+----------+----------+---------+----------+--+
    |    *******                                               |
    |   **     **                                              |
1.5 +  **       **                                             +
    |            **                                            |
    |              **                                          |
    |               **                                         |
  1 +                **                                        +
    |                 ***                                      |
    |                   ***                                    |
    |                     ***                          *****   |
0.5 +                       *****                *******   **  +
    |                           ******************             |
    +--+---------+----------+----------+---------+----------+--+
       0        0.2        0.4        0.6       0.8         1   
                            max_pool_dist                       
     +-+---------+----------+----------+----------+---------+--+
0.15 +  **                                                     +
     |  ***                                                    |
     |    *                                                    |
     |     *                                                   |
 0.1 +     *                                                   +
     |     **                                                  |
     |      *                                                  |
     |      **                                                 |
0.05 +       *                                                 +
     |       **          *                                     |
     |        *        ****        ****                        |
   0 +         *********   *********  ***********************  +
     +-+---------+----------+----------+----------+---------+--+
       0        10         20         30         40        50   
                             min_l_loci                         
     +-+---------+----------+----------+----------+---------+--+
0.04 +   *****                                                 +
     |  **   **        ****                                    |
     |  *     **      **   *                                   |
0.03 +         **   ***     *                                  +
     |          *****       **                                 |
     |                       **                                |
0.02 +                        **                               +
     |                         ***                             |
     |                           *****                         |
0.01 +                               ****                      +
     |                                  ****                   |
     |                                     ******************  |
     +-+---------+----------+----------+----------+---------+--+
       0        10         20         30         40        50   
                          min_k_neighbours                      
```
































































