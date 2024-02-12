library(testthat)
# library(imputef)
rextendr::document()
devtools::load_all()

fn_extract_missing = function(fname_imputed) {
    df = read.csv(fname_imputed)
    df = df[rowSums(df[4:ncol(df)])>0, ]
    idx = which(((df$X.chr=="chrA") & (df$pos==2309794) & (df$allele=="A")) | ((df$X.chr=="chrA") & (df$pos==2309996) & (df$allele=="A")))
    return(c(unlist(df[idx, 4:5]), dim(df)))
}

fname_vcf = "tests/test.vcf"
fname_sync = "tests/test.sync"
fname_csv = "tests/test.csv"

test_that(
    "mvi", {
        print("mvi:")
        vcf = fn_extract_missing(mvi(fname=fname_vcf)) - c(0, 0, 0, 0, 21, 0)
        sync = fn_extract_missing(mvi(fname=fname_sync)) - c(0, 0, 0, 0, 0, 0)
        csv = fn_extract_missing(mvi(fname=fname_csv)) - c(0, 0, 0, 0, 1, 0)
        expect_equal(vcf, sync, tolerance=0.1)
        expect_equal(vcf, csv, tolerance=0.1)
    }
)

test_that(
    "aldknni_fixed", {
        print("aldknni_fixed:")
        vcf = fn_extract_missing(aldknni(fname=fname_vcf, min_loci_corr=0.9, max_pool_dist=0.1))
        sync = fn_extract_missing(aldknni(fname=fname_sync, min_loci_corr=0.9, max_pool_dist=0.1))
        csv = fn_extract_missing(aldknni(fname=fname_csv, min_loci_corr=0.9, max_pool_dist=0.1)) - c(0, 0, 0, 0, 1, 0)
        expect_equal(vcf, sync, tolerance=0.1)
        expect_equal(vcf, csv, tolerance=0.1)
    }
)

test_that(
    "aldknni_optim", {
        print("aldknni_optim:")
        vcf = fn_extract_missing(aldknni(fname=fname_vcf))
        sync = fn_extract_missing(aldknni(fname=fname_sync))
        csv = fn_extract_missing(aldknni(fname=fname_csv)) - c(0, 0, 0, 0, 1, 0)
        expect_equal(vcf, sync, tolerance=0.1)
        expect_equal(vcf, csv, tolerance=0.1)
    }
)
