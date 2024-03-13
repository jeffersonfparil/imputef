# imputef

Impute allele frequencies to reduce sparsity of genotype data from polyploids, pooled individuals, and populations.

|**Build Status**|**License**|
|:--------------:|:---------:|
| <a href="https://github.com/jeffersonfparil/imputef/actions"><img src="https://github.com/jeffersonfparil/imputef/actions/workflows/rust.yml/badge.svg"></a> | [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) |

## Overview

- [Installation](#installation)
- [Usage](#usage)
- [Details](#details)
- [Optimisation](#optimisation-approaches)
- [Performance evaluation](#performance-evaluation)
- [Troubleshooting](#troubleshooting)
- [References](#references)
- [Acknowledgements](#acknowledgements)

## Installation

### Pre-compiled binaries

1. Download the appropriate executable binary compatible with your system

- [GNU/Linux (x86 64-bit)](https://github.com/jeffersonfparil/imputef/releases/download/v1.0.0/imputef-x86_64-linux)
- [Windows 10 (x86 64-bit)](https://github.com/jeffersonfparil/imputef/releases/download/v1.0.0/imputef-x86_64-windows.exe)
- ~~[macOS Catalina (x86 64-bit)](https://github.com/jeffersonfparil/imputef/releases/download/v1.0.0/DOESNOTEXIST)~~ (macOS binary pending)

2. Configure and execute

- In GNU/Linux (macOS binary pending):

```shell
chmod +x imputef
./imputef
```

- In Windows, open command prompt via: *Win + R*, type "cmd" and press enter. Navigate to your download folder and execute, e.g. `imputef-x86_64-windows.exe -h`.

### Compile from source

1. Clone the repository

```shell
git clone https://jeffersonfparil:<API_KEY>@github.com/jeffersonfparil/imputef.git main
```

2. Load the Rust development environment via Conda (please see [Conda installation instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) if you do not have Conda pre-installed)

```shell
cd imputef/
conda env create --file res/rustenv.yml
conda activate rustenv
```

3. Compile and optionally create an alias or a symbolic link or add to $PATH

```shell
cargo build --release
target/release/imputef -h
# Option 1: alias
echo alias imputef="$(pwd)/target/release/imputef" >> ~/.bashrc
source ~/.bashrc
# Option 2: symlink
sudo ln -s $(pwd)/target/release/imputef /usr/bin/imputef
# Option 3: add to $PATH (Note that you need to place this in your ~/.bashrc or the appropriate shell initialisation file)
export PATH=${PATH}:$(pwd)/target/release
### Check
type -a imputef
cd ~; imputef -h; cd -
```

## Usage

### Quickstart

```shell
imputef -h
imputef -f tests/test.csv   # allele frequency table as input
imputef -f tests/test.sync  # synchronised pileup file as input
imputef -f tests/test.vcf   # variant call format as input without missing data
imputef -f tests/test_2.vcf   # variant call format as input
imputef -f tests/test_2.vcf --method mean # use mean value imputation
imputef -f tests/test_2.vcf --min-loci-corr=0.75 --max-pool-dist=0.25 # define some minimum loci correlation and maximum genetic distance thresholds
```

### Arguments

- **fname**: Filename of the genotype file to be imputed in uncompressed [vcf](#variant-call-format-vcf), [sync](#synchronised-pileup-sync), or [allele frequency table](#allele-frequency-table-csv). Details on these genotype formats are available below.
- **method**: Imputation method. Use "mean" for mean value imputation or "aldknni" for adaptive LD-kNN imputation. [Default = "aldknni"]
- **min-coverage**: Minimum coverage per locus, i.e. if at a locus, a pool falls below this value (does not skip missing data, i.e. missing locus has a depth of zero), then the whole locus is omitted. Set this to zero if the vcf has been filtered and contains missing values, i.e. `./.` or `.|.`. [Default = 0]
- **min-allele-frequency**: Minimum allele frequency per locus, i.e. if at a locus, a pool has all its alleles below this value and/or above the additive complement of this value (skipping missing data), then the entire locus is omitted. [Default = 0.0001]
- **max-missingness-rate-per-locus**: Maximum fraction of pools missing per locus, i.e. if at a locus, there were more pools missing than the coverage dictated by this threshold, then the locus is omitted. [Default = 1.00]
- **pool-sizes**: Vector of pool sizes, i.e. the number of individuals included in each pool. Enter the pool sizes separated by commas `,`. This can also be set to a single arbitrarily large value like 100 for individual polyploids or if allele frequency estimates are expected to be accurate. [Default = "100.0"]
- **min-depth-below-which-are-missing**: Minimum depth at which loci with depth below this threshold are set to missing. Set to one if the input vcf has already been filtered and the loci beyond the depth thresholds have been set to missing, otherwise set to an integer above zero. [Default = 1.00]
- **max-depth-above-which-are-missing**: Maximum depth at which loci with depth above this threshold are set to missing. Set to some large arbitrarily large value (e.g. 1000000) if the input vcf has already been filtered and the loci beyond the depth thresholds have been set to missing, otherwise set to an integer above zero. [Default = 1000000.0]
- **frac-top-missing-pools**: Fraction of pools with the highest number of missing loci to be omitted. Set to zero if the input vcf has already been filtered and the loci beyond the depth thresholds have been set to missing, otherwise set to a decimal number between zero and one. [Default = 0.0]
- **frac-top-missing-loci**: Fraction of loci with the highest number of pools with missing data to be omitted. Set to zero if the input vcf has already been filtered and the loci beyond the depth thresholds have been set to missing, otherwise set to an decimal number between zero and one. [Default = 0.0]
- **min-loci-corr**: Minimum correlation (Pearson's correlation) between the locus requiring imputation and other loci deemed to be in linkage with it. Ranges from 0.0 to 1.0, but use -1.0 or any negative number to perform per locus optimisations to find the best value minimising imputation. [Default = 0.9]
- **max-pool-dist**: Maximum genetic distance (mean absolute difference in allele frequencies) between the pool or sample requiring imputation and pools or samples deemed to be the closest neighbours. Ranges from 0.0 to 1.0, but use -1.0 or any negative number to perform per locus optimisations to find the best value minimising imputation. [Default = 0.1]
- **min-l-loci**: Minimum number of linked loci to be used in estimating genetic distances between the pool or sample requiring imputation and other pools or samples (minimum value of 1). This argument overrides `min-loci-corr`, i.e. the minimum number of loci will be met regardless of the minimum loci correlation threshold. [Default = 20]
- **min-k-neighbours**: Minimum number of k-nearest neighbours of the pool or sample requiring imputation (minimum value of 1). This argument overrides `max-pool-dist`, i.e. the minimum number of k-nearest neighbours will be met regardless of the maximum genetic distance threshold. [Default = 5]
- **restrict-linked-loci-per-chromosome**: Restrict the choice of linked loci to within the chromosome the locus requiring imputation belongs to? [default: false] [Default = false; i.e. no flag]
- **n-reps**: Number of replications for the estimation of imputation accuracy in terms of mean absolute error (MAE). It is used to define the number of random non-missing samples to use as replicates for the estimation of MAE and optimisation (minimum value of 1). [Default = 10]
- **n-threads**: Number of computing threads or processor cores to use in the computations. [Default = 2]
- **fname-out-prefix**: Prefix of the output files including the [imputed allele frequency table](#allele-frequency-table-csv) (`<fname-out-prefix>-<time>-<random-id>-IMPUTED.csv`). [Default = ""; which corresponds to the name of the input genotype file]


## Genotype data formats

### Note

**Header line/s and comments should be prepended by '#'.**

### Variant call format (vcf)

Canonical variant calling or genotype data format for individual samples. This should include the `AD` field (allele depth), and may or may not have genotypes called (e.g. generated via bctools mpileup -a AD,DP ...). If the `GT` field is present but the `AD` field is absent, then each sample is assumed to be an individual diploid, i.e., neither a polyploid nor a pool. See [VCFv4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf) and [VCFv4.3](https://samtools.github.io/hts-specs/VCFv4.3.pdf) for details in the format specifications.

### Synchronised pileup (sync)
- an extension of [popoolation2's](https://academic.oup.com/bioinformatics/article/27/24/3435/306737) sync or synchronised pileup file format, which includes a header line prepended by '#' showing the names of each column including the names of each pool. Additional header line/s and comments prepended with '#' may be added anywhere within the file.
- tab-delimited
- *Header line/s*:  optional header line/s including the names of the pools, e.g. `# chr pos ref pool1 pool2 pool3 pool4 pool5`
- *Column 1*:       chromosome or scaffold name
- *Column 2*:       locus position 
- *Column 3*:       reference allele, e.g. A, T, C, G 
- *Column/s 4 to n*:  colon-delimited allele counts: A:T:C:G:DEL:N, where "DEL" refers to insertion/deletion, and "N" is unclassified. A pool or population or polyploid individual is represented by a single column of this colon-delimited allele counts.

### Allele frequency table (csv)
- comma-delimited
- *Header line*: ` #chr,pos,allele,<pool_name_1>,...,<pool_name_n>`
- each locus is represented by 2 or more rows, i.e. 2 for biallelic loci, and >2 for multi-allelic loci


## Details

Imputation of genotype data from sequencing of more than 2 sets of genomes, i.e. polyploid individuals, population samples, or pools of individuals. This tool can also perform simple genotype data filtering prior to imputation. Two imputation methods are available: (1) mean value imputation which uses the arithmetic mean of the locus across non-missing pools; (2) adaptive linkage-informed k-nearest neighbour imputation. This is an attempt to extend the [LD-kNNi method of Money et al, 2015, i.e. LinkImpute](https://doi.org/10.1534/g3.115.021667), which was an extension of the [kNN imputation of Troyanskaya et al, 2001](https://doi.org/10.1093/bioinformatics/17.6.520). Similar to LD-kNNi, LD is estimated using Pearson's product moment correlation across loci per pair of samples. Mean absolute difference in allele frequencies is used to define genetic distance between samples, instead of taxicab or Manhattan distance in LD-kNNi. Four parameters can be set by the user, (1) minimum loci correlation threshold: dictates the minimum LD between the locus requiring imputation and other loci which will be used to estimate genetic distance between samples; (2) maximum genetic distance threshold: sets the maximum genetic distance between the sample requiring imputation and the samples (i.e. nearest neighbours) to be used in weighted mean imputation of missing allele frequencies; (3) minimum number of loci linked to the locus requiring imputation: overrides minimum loci correlation threshold if this minimum is not met; and (4) minimum k-nearest neighbours: overrides maximum genetic distance threshold if this minimum is not met. The first two parameters (minimum loci correlation and maximum genetic distance thresholds) can be optimised per locus requiring imputation using non-missing samples as replicates simulating missing data to minimum the mean absolute error in imputation.

### `--method="mean"`: mean value imputation

This imputation uses the arithmetic mean of the observed allele frequencies across all samples where the locus was genotyped:

$$
\hat q_{r,j} = { {1 \over (n-m)} { \sum_{i \ne r}^{n} q_{i,j} } }
$$

where:

- $\hat q_{r,j}$ is the imputed allele frequency of sample $r$ at the $j^{\text {th}}$ locus,
- $n$ is the total number of samples,
- $m$ is the number of samples which are missing data at the $j^{\text {th}}$ locus, and
- $q_{i,j}$ is the known allele frequency of the $i^{\text {th}}$ sample at the $j^{\text {th}}$ locus.

### `--method="aldknni"`: adaptive linkage disequilibrium (LD)-based k-nearest neighbour imputation of genotype data

This is an attempt to extend the [LD-kNNi method of Money et al, 2015, i.e. LinkImpute](https://doi.org/10.1534/g3.115.021667), which was an extension of the [kNN imputation of Troyanskaya et al, 2001](https://doi.org/10.1093/bioinformatics/17.6.520). Similar to LD-kNNi, linkage disequilibrium (LD) is estimated using Pearson's product moment correlation per pair of loci, which is computed per chromosome by default, but can be computed across the entire genome. We use the mean absolute difference/error (MAE) between allele frequencies among linked loci as an estimate of genetic distance between samples. Fixed values for the minimum correlation to identify loci used in distance estimation, and maximum genetic distance to select the k-nearest neighbours can be defined. Additionally, minimum number of loci to include in distance estimation, and minimum number of nearest neighbours can be set. Moreover, all four parameters can be optimised, i.e. the minimum correlation and/or maximum distance and/or minimum number of loci and/or minimum number of nearest neighbours which minimises the MAE between predicted and expected allele frequencies after simulating 10% missing data are identified.

The allele depth information (`AD`), i.e. the unfiltered allele depth which includes the reads which did not pass the variant caller filters are used to calculate allele frequencies. If the `GT` field is present but the `AD` field is absent, then each sample is assumed to be an individual diploid, i.e., neither a polyploid nor a pool. Optional filtering steps based on minimum depth, minimum allele frequency, and maximum sparsity are available. LD estimation and imputation are multi-threaded, and imputation output is written into disk as an [allele frequency table](#allele-frequency-table-csv). The structs, traits, methods, and functions defined in this library are subsets of [poolgen](https://github.com/jeffersonfparil/poolgen), and will eventually be merged. 

The imputed allele frequency is computed as:

$$
\hat q_{r,j} = { \sum_{i \ne r}^{k} q_{i,j} (1 - \delta_{i,r}) }
$$

with:

$$
\delta_{i,r} = { {1 \over \sum d_{i,r}} d_{i,r} }
$$

and

$$
d_{i,r} = { {1 \over c} { \sum_{j=1}^{c} |q_{i,j} - q_{r,j}| } }
$$

where:

- $\hat q_{r,j}$ is the imputed allele frequency of sample $r$ at the $j^{\text {th}}$ locus,
- $n$ is the total number of samples,
- $m$ is the number of samples which are missing data at the $j^{\text {th}}$ locus,
- $q_{i,j}$ is the known allele frequency of the $i^{\text {th}}$ sample at the $j^{\text {th}}$ locus,
- $k$ is the number of nearest neighbours or the samples most closely related to the sample requiring imputation, i.e. sample $r$ at locus $j$, and
- $\delta_{i,r}$ is scaled $d_{i,r}$ which is the genetic distance between the $i^{\text {th}}$ sample and sample $r$. This distance is the mean absolute difference in allele frequencies between the two samples across $c$ linked loci.

The variables $k$ and $c$ are proportional to the user inputs `max_pool_dist` (default=0.1) and `min_loci_corr` (default=0.9), respectively. The former defines the maximum distance of samples to be considered as one of the k-nearest neighbours, while the latter refers to the minimum correlation with the locus requiring imputation to be included in the estimation of the genetic distance.

## Optimisation approaches

1. Using the built-in grid-search optimisation for minimum loci correlation and maximum genetic distance thresholds per locus:

```shell
imputef -f tests/test_2.vcf --min-loci-corr=-1 --max-pool-dist=-1
```

2. Grid-search optimisation for all four parameters which assumes a common set of optimal parameters across all loci

```shell
### Find the optimal l-loci to use for distance estimation, k-nearest neighbours, minimum loci correlation, and maximum distance
### i.e. the combination of these 4 parameters which minimises the mean absolute error between expected and predicted allele frequencies.
echo 'l,k,corr,dist,mae' > grid_search.csv
for l in 5 10 15
do
    for k in 1 2 3 4 5
    do
        for corr in 0.75 0.95 1.00
        do
            for dist in 0.0 0.10 0.25
            do
                echo "@@@@@@@@@@@@@@@@@@@@@@@@"
                echo ${l},${k},${corr},${dist}
                imputef -f tests/test_2.vcf \
                    --min-l-loci=${l} \
                    --min-k-neighbours=${k} \
                    --min-loci-corr=${corr} \
                    --max-pool-dist=${dist} \
                    --n-reps=3 > log.tmp
                fname_out=$(tail -n1 log.tmp | cut -d':' -f2 | tr -d ' ')
                mae=$(grep "Expected imputation accuracy in terms of mean absolute error:" log.tmp | cut -d':' -f2 | tr -d ' ')
                echo ${l},${k},${corr},${dist},${mae} >> grid_search.csv
                rm $fname_out log.tmp
            done
        done
    done
done
awk -F',' 'NR == 1 || $5 < min {min = $5; min_line = $0} END {print min_line}' grid_search.csv
```

## Performance evaluation

Datasets: 

- autotetraploid *Medicago sativa* (2n=4x=32; 2.74 Gb genome; 155 samples x 124,151 biallelic loci; in-house source)
- diploid *Vitis vinifera* (2n=2x=38; 0.5 Gb genome; 77 samples x 8,506 loci biallelic; [Money et al., 2015](https://doi.org/10.1534/g3.115.021667)) with the 2.90% missing data prior to sparsity simulations
- pools of diploid *Glycine max* (2n=2x=20; 1.15 Gb genome; 478 pools (each pool comprised of 42 individuals) x 39,636 biallelic loci; source: [http://gong_lab.hzau.edu.cn/Plant_imputeDB/#!/download_soybean](http://gong_lab.hzau.edu.cn/Plant_imputeDB/#!/download_soybean))

Performance metrics:

- Concordance: $c = {{1 \over n} \Sigma_{i=1}^{n} p}$, where: $p=
\begin{cases}
0 \text{ if } \hat g \ne g_{true}\\
1 \text{ if } \hat g = g_{true}
\end{cases}
$.
This is used for genotype classes, i.e., binned allele frequencies: $g = {{1 \over {ploidy}} round(q*ploidy)}$, here $q = P(allele)$. Note that there is alternative way of defining these genotype classes with strict boundaries, i.e., homozygotes have fixed allele frequencies.
- Mean absolute error: $mae = {{1 \over n} \Sigma_{i=1}^{n}|\hat q - q_{true}|}$.
- Coefficient of determination: $R^2 = { 1 - {{\Sigma_{}^{}(\hat q - q_{true})^2} \over {\Sigma_{}^{}(\hat q_{true} - \bar q_{true})^2}} }$

------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------

### Autotetraploid (Lucerne) mean absolute error

![mae_barplots](./res/lucerne-Mean_absolute_error.svg)

### Pool (Soybean pools) mean absolute error

![mae_barplots](./res/soybean-Mean_absolute_error.svg)

### Diploid (Grape) mean absolute error

![mae_barplots](./res/grape-Mean_absolute_error.svg)

------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------

### Autotetraploid (Lucerne) concordance of observed and imputed genotype classes

![concordance_genotype_classes_barplots](./res/lucerne-Concordance.svg)

### Pool (Soybean pools) concordance of observed and imputed genotype classes

![concordance_genotype_classes_barplots](./res/soybean-Concordance.svg)

### Diploid (Grape) concordance of observed and imputed genotype classes

![concordance_genotype_classes_barplots](./res/grape-Concordance.svg)

------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------

## Troubleshooting

- Out-of-memory (OOM) error is likely due to the pairwise LD estimation across the genome. We have taken steps to reduce the likelihood of this happening by using `u8` instead of `f64` in calculating Pearson's correlation between pairs of loci. However, if you encountered this OOM error, please consider using the flag `--restrict-linked-loci-per-chromosome` to estimate pairwise LD per chromosome only. This assumes that you have a dense coverage of the genome, i.e., there are enough markers with non-missing data to accurately determine relationships between loci and samples to yield good imputation accuracies.

## References

- Money D, Gardner K, Migicovsky Z, Schwaninger H, Zhong GY, Myles S. LinkImpute: fast and accurate genotype imputation for nonmodel organisms. G3: Genes|Genomes|Genetics. 2015;5(11):2383–90. doi:10.1534/g3.115.021667.
- Troyanskaya O, Cantor M, Sherlock G, Brown P, Hastie T et al. , 2001 Missing value estimation methods for DNA microarrays. Bioinformatics 17: 520–525.
- Schwender H, 2012 Imputing missing genotypes with weighted k nearest neighbors. J. Toxicol. Environ. Health A 75: 438–446.

## Acknowledgements

This work was conceived and developed during my employment in Agriculture Victoria. The imputation algorithm in this repository was inspired by the algorithms presented in the 3 papers above and [Luke Pembletton](https://github.com/lpembleton)'s tetraploid imputation algorithm written in R. The core data structures, traits, and methods are largely shared with my open-source (GPLv3) project [poolgen](https://github.com/jeffersonfparil/poolgen).
