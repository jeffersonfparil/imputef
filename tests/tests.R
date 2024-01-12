library(testthat)
# library(imputef)
rextendr::document()
devtools::load_all()

MVI_impute_2_loci = function(fname_vcf_sync_csv) {
    out = mvi(fname=fname_vcf_sync_csv)
    df = read.csv(out)
    idx_1 = which((df$X.chr=="chrA") & (df$pos==2309794))
    idx_2 = which((df$X.chr=="chrA") & (df$pos==2309996))
    vcf_idx_1 = df$Entry.0[idx_1]
    vcf_idx_2 = df$Entry.1[idx_2]
    return(c(vcf_idx_1, vcf_idx_2, nrow(df), ncol(df)))
}

ALDKNNI_FIXED_LK_impute_2_loci = function(fname_vcf_sync_csv) {
    out = aldknni(fname=fname_vcf_sync_csv, min_loci_corr=10, max_pool_dist=3)
    df = read.csv(out)
    idx_1 = which((df$X.chr=="chrA") & (df$pos==2309794))
    idx_2 = which((df$X.chr=="chrA") & (df$pos==2309996))
    vcf_idx_1 = df$Entry.0[idx_1]
    vcf_idx_2 = df$Entry.1[idx_2]
    return(c(vcf_idx_1, vcf_idx_2, nrow(df), ncol(df)))
}

ALDKNNI_FIXED_CORRDIST_impute_2_loci = function(fname_vcf_sync_csv) {
    out = aldknni(fname=fname_vcf_sync_csv, min_loci_corr=0.75, max_pool_dist=0.25)
    df = read.csv(out)
    idx_1 = which((df$X.chr=="chrA") & (df$pos==2309794))
    idx_2 = which((df$X.chr=="chrA") & (df$pos==2309996))
    vcf_idx_1 = df$Entry.0[idx_1]
    vcf_idx_2 = df$Entry.1[idx_2]
    return(c(vcf_idx_1, vcf_idx_2, nrow(df), ncol(df)))
}

ALDKNNI_OPTIM_LK_impute_2_loci = function(fname_vcf_sync_csv) {
    out = aldknni(fname=fname_vcf_sync_csv, min_loci_corr=0, max_pool_dist=0, optimise_for_thresholds=FALSE, optimise_n_reps=1)
    df = read.csv(out)
    idx_1 = which((df$X.chr=="chrA") & (df$pos==2309794))
    idx_2 = which((df$X.chr=="chrA") & (df$pos==2309996))
    vcf_idx_1 = df$Entry.0[idx_1]
    vcf_idx_2 = df$Entry.1[idx_2]
    return(c(vcf_idx_1, vcf_idx_2, nrow(df), ncol(df)))
}

ALDKNNI_OPTIM_CORRDIST_impute_2_loci = function(fname_vcf_sync_csv) {
    out = aldknni(fname=fname_vcf_sync_csv, min_loci_corr=0, max_pool_dist=0, optimise_for_thresholds=TRUE, optimise_n_reps=1)
    df = read.csv(out)
    idx_1 = which((df$X.chr=="chrA") & (df$pos==2309794))
    idx_2 = which((df$X.chr=="chrA") & (df$pos==2309996))
    vcf_idx_1 = df$Entry.0[idx_1]
    vcf_idx_2 = df$Entry.1[idx_2]
    return(c(vcf_idx_1, vcf_idx_2, nrow(df), ncol(df)))
}

tests = function() {
    # test_that(
    #     "mvi", {
    #         print("mvi:")
    #         vcf = MVI_impute_2_loci("tests/test.vcf")
    #         sync = MVI_impute_2_loci("tests/test.sync")
    #         csv = MVI_impute_2_loci("tests/test.csv")
    #         expect_equal(vcf, c(0.018519, 0.981481, 0.084808, 0.915192, 753, 13))
    #         expect_equal(vcf, sync)
    #         expect_equal(vcf, csv - c(0, 0, 0, 0, 1, 0)) ### reduce the number of alleles across loci by one as the csv or geno format file had one locus appended with an additional missing alternative allele
    #     }
    # )
    test_that(
        "aldknni_fixed_lk", {
            print("aldknni_fixed_lk:")
            vcf = ALDKNNI_FIXED_LK_impute_2_loci("tests/test.vcf")
            sync = ALDKNNI_FIXED_LK_impute_2_loci("tests/test.sync")
            csv = ALDKNNI_FIXED_LK_impute_2_loci("tests/test.csv")
            expect_equal(vcf, c(0.019371, 0.980629, 0.086328, 0.913672, 753, 13))
            expect_equal(vcf, sync)
            expect_equal(vcf, csv - c(0, 0, 0, 0, 1, 0)) ### reduce the number of alleles across loci by one as the csv or geno format file had one locus appended with an additional missing alternative allele
        }
    )
    test_that(
        "aldknni_fixed_corrdist", {
            print("aldknni_fixed_corrdist:")
            vcf = ALDKNNI_FIXED_CORRDIST_impute_2_loci("tests/test.vcf")
            sync = ALDKNNI_FIXED_CORRDIST_impute_2_loci("tests/test.sync")
            csv = ALDKNNI_FIXED_CORRDIST_impute_2_loci("tests/test.csv")
            expect_equal(vcf, c(0.017927, 0.982073, 0.084668, 0.915332, 753, 13))
            expect_equal(vcf, sync)
            expect_equal(vcf, csv - c(0, 0, 0, 0, 1, 0), tolerance=1e-5) ### reduce the number of alleles across loci by one as the csv or geno format file had one locus appended with an additional missing alternative allele
        }
    )
    test_that(
        "aldknni_optim_lk", {
            print("aldknni_optim_lk:")
            vcf = ALDKNNI_OPTIM_LK_impute_2_loci("tests/test.vcf")
            sync = ALDKNNI_OPTIM_LK_impute_2_loci("tests/test.sync")
            csv = ALDKNNI_OPTIM_LK_impute_2_loci("tests/test.csv")
            expect_equal(vcf, c(0.0, 1.0, 0.0, 1.0, 753, 13), tolerance=0.1)
            expect_equal(vcf, sync, tolerance=0.2)
            expect_equal(vcf, csv - c(0, 0, 0, 0, 1, 0), tolerance=0.2) ### reduce the number of alleles across loci by one as the csv or geno format file had one locus appended with an additional missing alternative allele
        }
    )
    test_that(
        "aldknni_optim_corrdist", {
            print("aldknni_optim_corrdist:")
            vcf = ALDKNNI_OPTIM_CORRDIST_impute_2_loci("tests/test.vcf")
            sync = ALDKNNI_OPTIM_CORRDIST_impute_2_loci("tests/test.sync")
            csv = ALDKNNI_OPTIM_CORRDIST_impute_2_loci("tests/test.csv")
            expect_equal(vcf, c(0.0, 1.0, 0.0, 1.0, 753, 13), tolerance=0.1)
            expect_equal(vcf, sync, tolerance=0.2)
            expect_equal(vcf, csv - c(0, 0, 0, 0, 1, 0), tolerance=0.2) ### reduce the number of alleles across loci by one as the csv or geno format file had one locus appended with an additional missing alternative allele
        }
    )
}

tests()
