# Core architecture of the genomic_selection code repository

## Dependency map

```shell
src/
├── main.rs
|   ├── structs_and_traits.rs
|   ├── sync.rs
|   │   ├── structs_and_traits.rs
|   │   └── helpers.rs
|   ├── vcf.rs
|   │   ├── structs_and_traits.rs
|   │   └── helpers.rs
|   ├── geno.rs
|   │   ├── structs_and_traits.rs
|   │   └── helpers.rs
|   ├── mvi.rs
|   │   ├── structs_and_traits.rs
|   │   └── helpers.rs
|   ├── aldknni.rs
|   │   ├── structs_and_traits.rs
|   │   └── helpers.rs
|   └── filter_missing.rs
|       ├── structs_and_traits.rs
|       └── helpers.rs
└── phen.rs
    └── structs_and_traits.rs
```

## Module description

1. [`main.R`](src/main.R) - the user-interface script with arguments parsing and help docs (i.e. `Rscript src/main.R -h`)
2. [`load.R`](src/load.R) - loads phenotype and covariate data in character-delimited text formats, as well as genotype data in variant call format (vcf), and binary R data format (Rds)


0. [`phen.rs`](src/phen.rs) - non-critical in this project as we do not need phenotype data for imputation. This is reserved for quantitative and population genetics analyses in [poolgen](https://github.com/jeffersonfparil/poolgen).

