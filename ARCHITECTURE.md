# Architecture of imputef

```shell
src/
├── main.rs
│   ├── structs_and_traits.rs
│   ├── helpers.rs
│   ├── sync.rs
│   ├── vcf.rs
│   ├── geno.rs
│   ├── mvi.rs
│   ├── aldknni.rs
│   └── filter_missing.rs
└── phen.rs
```

## Module description

1. [`main.rs`](src/main.rs) - the main script handling the user inputs via [`clap`](https://github.com/clap-rs/clap).
2. [`structs_and_traits.rs`](src/structs_and_traits.rs) - one-stop-shop for all the structs and traits used in the entire project
3. [`helpers.rs`](src/helpers.rs) - helper functions used throughout the project
4. [`sync.rs`](src/sync.rs) - main genotype data parsing, filtering, and writing methods
5. [`vcf.rs`](src/vcf.rs) - vcf file parsing
6. [`geno.rs`](src/geno.rs) - allele frequency table file parsing
7. [`mvi.rs`](src/mvi.rs) - mean value imputation including missing data simulation to estimate expected imputation accuracy for this imputation method
8. [`aldknni.rs`](src/aldknni.rs) - allele frequency LD-kNN imputation methods including linkage and genetic distance estimation
9. [`filter_missing.rs`](src/filter_missing.rs) - genotype data filtering by depth, and sparsity per locus or sample
10. [`phen.rs`](src/phen.rs) - non-critical in this project as we do not need phenotype data for imputation. This is reserved for quantitative and population genetics analyses in [poolgen](https://github.com/jeffersonfparil/poolgen).

