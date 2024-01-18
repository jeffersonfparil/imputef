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
        expect_equal(vcf, sync, tolerance=0.05)
        expect_equal(vcf, csv, tolerance=0.05)
    }
)

test_that(
    "aldknni_fixed", {
        print("aldknni_fixed:")
        vcf = fn_extract_missing(aldknni(fname=fname_vcf)) - c(0, 0, 0, 0, 21, 0)
        sync = fn_extract_missing(aldknni(fname=fname_sync)) - c(0, 0, 0, 0, 0, 0)
        csv = fn_extract_missing(aldknni(fname=fname_csv)) - c(0, 0, 0, 0, 1, 0)
        expect_equal(vcf, sync, tolerance=0.05)
        expect_equal(vcf, csv, tolerance=0.05)
    }
)

test_that(
    "aldknni_optim_cd", {
        print("aldknni_optim_cd:")
        vcf = fn_extract_missing(aldknni(fname=fname_vcf, optimise_n_steps_min_loci_corr=10, optimise_n_steps_max_pool_dist=10, min_l_loci=1, min_k_neighbours=1, optimise_max_l_loci=1, optimise_max_k_neighbours=1)) - c(0, 0, 0, 0, 21, 0)
        sync = fn_extract_missing(aldknni(fname=fname_sync, optimise_n_steps_min_loci_corr=10, optimise_n_steps_max_pool_dist=10, min_l_loci=1, min_k_neighbours=1, optimise_max_l_loci=1, optimise_max_k_neighbours=1)) - c(0, 0, 0, 0, 0, 0)
        csv = fn_extract_missing(aldknni(fname=fname_csv, optimise_n_steps_min_loci_corr=10, optimise_n_steps_max_pool_dist=10, min_l_loci=1, min_k_neighbours=1, optimise_max_l_loci=1, optimise_max_k_neighbours=1)) - c(0, 0, 0, 0, 1, 0)
        expect_equal(vcf, sync, tolerance=0.05)
        expect_equal(vcf, csv, tolerance=0.05)
    }
)

test_that(
    "aldknni_optim_lk", {
        print("aldknni_optim_lk:")
        vcf = fn_extract_missing(aldknni(fname=fname_vcf, min_loci_corr=0.0, max_pool_dist=1.0, optimise_n_steps_min_loci_corr=1, optimise_n_steps_max_pool_dist=1, optimise_max_l_loci=100, optimise_max_k_neighbours=100)) - c(0, 0, 0, 0, 21, 0)
        sync = fn_extract_missing(aldknni(fname=fname_sync, min_loci_corr=0.0, max_pool_dist=1.0, optimise_n_steps_min_loci_corr=1, optimise_n_steps_max_pool_dist=1, optimise_max_l_loci=100, optimise_max_k_neighbours=100)) - c(0, 0, 0, 0, 0, 0)
        csv = fn_extract_missing(aldknni(fname=fname_csv, min_loci_corr=0.0, max_pool_dist=1.0, optimise_n_steps_min_loci_corr=1, optimise_n_steps_max_pool_dist=1, optimise_max_l_loci=100, optimise_max_k_neighbours=100)) - c(0, 0, 0, 0, 1, 0)
        expect_equal(vcf, sync, tolerance=0.05)
        expect_equal(vcf, csv, tolerance=0.05)
    }
)

test_that(
    "aldknni_optim_all", {
        print("aldknni_optim_all:")
        vcf = fn_extract_missing(aldknni(fname=fname_vcf, optimise_n_steps_min_loci_corr=10, optimise_n_steps_max_pool_dist=10, optimise_max_l_loci=100, optimise_max_k_neighbours=100)) - c(0, 0, 0, 0, 21, 0)
        sync = fn_extract_missing(aldknni(fname=fname_sync, optimise_n_steps_min_loci_corr=10, optimise_n_steps_max_pool_dist=10, optimise_max_l_loci=100, optimise_max_k_neighbours=100)) - c(0, 0, 0, 0, 0, 0)
        csv = fn_extract_missing(aldknni(fname=fname_csv, optimise_n_steps_min_loci_corr=10, optimise_n_steps_max_pool_dist=10, optimise_max_l_loci=100, optimise_max_k_neighbours=100)) - c(0, 0, 0, 0, 1, 0)
        expect_equal(vcf, sync, tolerance=0.05)
        expect_equal(vcf, csv, tolerance=0.05)
    }
)
