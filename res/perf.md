# Assessing imputation accuracies

## Variables

### 3 imputation algorithms

1. **MVI**: mean value imputation
2. **AFIXED**: adaptive LD-kNN imputation using fixed/default minimum loci correlation (0.9), and maximum pool distance (0.1)
3. **AOPTIM**: adaptive LD-kNN imputation with optimisation for the minimum loci correlation and maximum pool distance per locus requiring imputation
4. **LINKIMPUTE**: LD-kNN imputation designed for diploids (results presented only for the individual diploid dataset, i.e. grape data)

### 2 minor allele frequency thresholds:

1. 0.01
2. 0.05

### 10 sparsity levels (missing rate):
1. 0.01
2. 0.1
3. 0.2
4. 0.3
5. 0.4
6. 0.5
7. 0.6
8. 0.7
9. 0.8
10. 0.9

## Datasets

1. autotetraploid *Medicago sativa* (2n=4x=32; 2.74 Gb genome; 155 samples x 124,151 biallelic loci; in-house source)
2. pools of diploid *Glycine max* (2n=2x=20; 1.15 Gb genome; 478 pools (each pool comprised of 42 individuals) x 39,636 biallelic loci; source: [http://gong_lab.hzau.edu.cn/Plant_imputeDB/#!/download_soybean](http://gong_lab.hzau.edu.cn/Plant_imputeDB/#!/download_soybean))
3. diploid *Vitis vinifera* (2n=2x=38; 0.5 Gb genome; 77 samples x 8,506 biallelic loci; source: [021667_FileS1 - zip file](https://academic.oup.com/g3journal/article/5/11/2383/6025349#supplementary-data))


### prepping_datasets.sh

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
    soybean.vcf
mv LinkImpute.data.grape.num.raw.txt.vcf grape.vcf
```

#### poolify.R

Generate pools from soybean individual genotype data, where we assume each individual per pool are perfectly equally represented (best-case scenario in the real-world).

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

#### ssv2vcf.R

Convert space-delimited genotype data from LinkImpute paper.

```R
args = commandArgs(trailingOnly = TRUE)
# args = c("/group/pasture/Jeff/imputef/misc/LinkImpute.data.apple.num.raw.txt", "/group/pasture/Jeff/imputef/misc/soybean.vcf")
# args = c("/group/pasture/Jeff/imputef/misc/LinkImpute.data.grape.num.raw.txt", "/group/pasture/Jeff/imputef/misc/soybean.vcf")
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


### Assess the genetic relationships between samples per dataset

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

1. Concordance: $c = {{1 \over n} \Sigma_{i=1}^{n} p}$, where: $p=
\begin{cases}
0 \text{ if } \hat g \ne g_{true}\\
1 \text{ if } \hat g = g_{true}
\end{cases}
$.
This is used for genotype classes, i.e., binned allele frequencies: $g = {{1 \over {ploidy}} round(q*ploidy)}$, here $q = P(allele)$. Note that there is alternative way of defining these genotype classes with strict boundaries, i.e., homozygotes have fixed allele frequencies.
2. Mean absolute error: $mae = {{1 \over n} \Sigma_{i=1}^{n}|\hat q - q_{true}|}$.
3. Coefficient of determination: $R^2 = { 1 - {{\Sigma_{}^{}(\hat q - q_{true})^2} \over {\Sigma_{}^{}(\hat q_{true} - \bar q_{true})^2}} }$


## Execution

```shell
### Create slurm scripts with specific memory and time limits per dataset
DIR=/group/pasture/Jeff/imputef/res
cd $DIR
for DATASET in "grape" "lucerne" "soybean"
do
    if [ $DATASET == "grape" ]
    then
        sed 's/--job-name="imputef"/--job-name="grapeImp"/g' perf.slurm | \
            sed 's/--mem=250G/--mem=100G/g' | \
            sed 's/--time=14-0:0:00/--time=0-0:30:00/g' > perf_${DATASET}.slurm
    elif [ $DATASET == "lucerne" ]
    then
        sed 's/--job-name="imputef"/--job-name="lucerImp"/g' perf.slurm | \
            sed 's/--time=14-0:0:00/--time=15-0:0:00/g' > perf_${DATASET}.slurm
    else
        sed 's/--job-name="imputef"/--job-name="soyImp"/g' perf.slurm | \
            sed 's/--mem=250G/--mem=200G/g' > perf_${DATASET}.slurm
    fi
done

### Submit array jobs for each dataset
for DATASET in "grape" "lucerne" "soybean"
do
    if [ $DATASET == "grape" ]
    then
        INI=1
        FIN=20
    elif [ $DATASET == "lucerne" ]
    then
        INI=21
        FIN=40
    else
        INI=41
        FIN=60
    fi
    echo ${DATASET}: ${INI}-${FIN}
    sbatch --array=${INI}-${FIN} perf_${DATASET}.slurm
done

### Monitor the jobs
conda activate rustenv
DIR=/group/pasture/Jeff/imputef/res
cd $DIR
squeue -u jp3h | sort
# tail slurm-2668461*_*.out
tail slurm-27726996_*.out
# grep -n -i "err" slurm-2668461*_*.out | grep -v "mean absolute" | grep -v "slurm_get_node_energy" | grep -v "_get_joules_task"
grep -n -i "err" slurm-27726996_*.out | grep -v "mean absolute" | grep -v "slurm_get_node_energy" | grep -v "_get_joules_task"
wc -l *-performance_assessment-maf_*missing_rate_*.csv
ls -lhtr
time Rscript perf_plot.R ${DIR}

### After all jobs have finished, move the output and plot:
mkdir output
mv *-performance_assessment-maf_*-missing_rate_*.csv output/
time Rscript perf_plot.R \
    ${DIR}/output
```

## Take-home message

The adaptive LD-kNN imputation algorithm works reasonably well even across the entire range of sparsity levels (0.1% to 90% missing data).

Note that the discrepancy between our imputation algorithm and LinkImpute's algorithm in the context of imputing binary genotypes is attributed to 3 differences:
- the main one is our optimisation per locus, i.e. per locus with at least one pool missing allele frequencies, we simulate pools to be missing data and perform a pseudo-grid search across combinations of the `min_loci_corr` and `max_pool_dist` threshold values which breaks out of the inner loop (`max_pool_dist`) if estimated MAE is larger than the current lowest MAE, then continues the outer loop to restart the inner looop and so on.
- the use of mean weighted allele frequencies in our case and weighted mode genotype class in LinkImpute, and
- the use of mean absolute error to measure accuracy in our case and concordance in LinkImpute.


## Miscellaneous

We are finding that the estiamtes of the imputation accuracy during imputation is underestimated at low sparsity and overestimated at high sparsity. Below we will try to fit a model to correct for these discrepancies.

First, let's extract the actual and estimated imputation error from the all the datasets:

```shell
conda activate rustenv
DIR=/group/pasture/Jeff/imputef/res
cd $DIR
echo "$(head -n1 $(ls *-performance_assessment-*.csv | head -n1)),estimated_mae" > misc_getting_more_accurate_imputation_error_estimates.csv
for f in $(ls *-performance_assessment-*.csv)
do
    # f=$(ls *-performance_assessment-*.csv | head -n1)
    maf_sparsity=$(echo $f | cut -d'-' -f3-4 | sed 's/maf_/maf/g' | sed 's/rate_/rate/g' | sed 's/.csv//g')
    slurm_out=$(grep "$maf_sparsity" slurm-26684613_*.out | head -n1 | cut -d':' -f1)
    grep -i "err" $slurm_out | cut -d':' -f2 | sed 's/ //g' > tmp_0.tmp
    tail -n+2  $f | grep -v -i "linkimpute" > tmp_1.tmp
    paste -d',' tmp_1.tmp tmp_0.tmp > tmp_2.tmp
    cat tmp_2.tmp >> misc_getting_more_accurate_imputation_error_estimates.csv
done
bat --wrap never misc_getting_more_accurate_imputation_error_estimates.csv
```

Now, let us proceed to modelling the imputation accuracy discrepancies:

```R
df = read.csv("misc_getting_more_accurate_imputation_error_estimates.csv")
str(df)

for (maf in unique(df$maf)) {
    # maf = unique(df$maf)[1]
    for (algo in unique(df$algorithm)) {
        # algo = unique(df$algorithm)[1]
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        print(paste0(algo, "-maf=", maf))
        idx = which((df$algorithm == algo) & (df$maf == maf))
        dat = data.frame(true=df$mae_frequencies[idx], pred=df$estimated_mae[idx])
        txtplot::txtplot(dat$true, dat$pred, xlab="true", ylab="estimated")
        vec_i = c(1:10)
        vec_r2 = c()
        for (i in vec_i) {
            mod = lm(true ~poly(pred, i, raw=TRUE), data=dat)
            # summary(mod)
            y_hat = predict(mod, newx=data.frame(pred=dat$pred))
            # txtplot::txtplot(dat$true, y_hat)
            vec_r2 = c(vec_r2, summary(mod)$adj.r.sq)
        }
        # cbind(vec_i, vec_r2)
        txtplot::txtplot(vec_i, vec_r2, xlab="degree", ylab="corr")
        degree = which(round(floor(vec_r2*100)/100, 2) == 0.99)[1]
        if (is.na(degree)) {
            degree = which(vec_r2 == max(vec_r2))[1]
        }
        mod = lm(true ~poly(pred, degree, raw=TRUE), data=dat)
        # summary(mod)
        y_hat = predict(mod, newx=data.frame(pred=dat$pred))
        txtplot::txtplot(dat$true, y_hat, xlab="true", ylab="y_hat")
        cor(dat$true, y_hat)
        cbind(true=dat$true, y_hat=y_hat, pred=dat$pred)

        b_hat = coef(mod)
        
        #######################################################################
        ### These below work well enough across maf and sparsity
        b_hat = c(0.05993495, -0.02031803, 0.49929972, -1.28096468, 1.01944518)
        degree = 4
        #######################################################################

        X = rep(1, nrow(dat))
        for (i in 1:degree) {
            X = cbind(X, dat$pred^i)
        }
        y_hat - (X %*% b_hat)
        print(cor(dat$true, y_hat))
        print(degree)
        print(b_hat)
    }
}
```
