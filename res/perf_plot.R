args = commandArgs(trailingOnly=TRUE)
# args = c("/group/pasture/Jeff/imputef/res")
dir_ouput = args[1]
setwd(dir_ouput)

### Plotting function per dataset
plot_metrics = function(df, dataset, vec_2_metrics=c("mae_frequencies", "concordance_classes"), vec_2_metrics_labels=c("Mean absolute error", "Concordance")) {
  # vec_fnames = list.files(path=".", pattern="-performance_assessment-maf_")
  # vec_fnames = vec_fnames[grepl("-missing_rate_", vec_fnames)]
  # vec_fnames = vec_fnames[grepl(".csv$", vec_fnames)]
  # dataset = "grape"
  # vec_files = vec_fnames[grepl(paste0("^", dataset), vec_fnames)]
  # for (i in 1:length(vec_files)) {
  #   df_tmp = read.csv(vec_files[i])
  #   if (i==1) {
  #     df = df_tmp
  #   } else {
  #     df = rbind(df, df_tmp)
  #   }
  # }
  # vec_2_metrics = c("mae_frequencies", "concordance_classes"); vec_2_metrics_labels = c("Mean absolute error", "Concordance")
  ### Load allele frequencies
  source("perf_functions.R")
  vcf = vcfR::read.vcfR(paste0("../misc/", dataset, ".vcf"), verbose=TRUE)
  mat_genotypes = fn_extract_allele_frequencies(vcf=vcf, min_depth=10, max_depth=200)$mat_genotypes
  vec_mu_af = colMeans(mat_genotypes, na.rm=TRUE)
  idx_af = which((vec_mu_af>=0.01) & (vec_mu_af<=(1-0.01)))
  vec_af = as.vector(mat_genotypes[, idx_af])
  vec_af = vec_af[!is.na(vec_af)]
  vec_af = c(vec_af, 1-vec_af)
  txtplot::txtdensity(vec_af)
  ### Sort algorithms according to increasing complexity and with LinkImpute at the bottom as it will not be assessed for polyploid and pool datasets
  df$algorithm = as.character(df$algorithm)
  df$algorithm[df$algorithm=="mvi"] = "a"
  df$algorithm[df$algorithm=="aldknni_fixed"] = "b"
  df$algorithm[df$algorithm=="aldknni_optim"] = "c"
  df$algorithm[df$algorithm=="linkimpute"] = "d"
  df$algorithm = as.factor(df$algorithm)
  levels(df$algorithm)[levels(df$algorithm)=="a"] = "MVI"
  levels(df$algorithm)[levels(df$algorithm)=="b"] = "AFIXED"
  levels(df$algorithm)[levels(df$algorithm)=="c"] = "AOPTIM"
  levels(df$algorithm)[levels(df$algorithm)=="d"] = "LINKIMPUTE"
  ### Mean concordance
  agg_concordance = aggregate(concordance_classes ~ algorithm, data=df, FUN=mean)
  agg_concordance = agg_concordance[order(agg_concordance$concordance_classes, decreasing=TRUE), ]
  print(agg_concordance)
  ### Mean computation time
  agg_duration = aggregate(duration_mins ~ algorithm, data=df, FUN=mean)
  agg_duration = agg_duration[order(agg_duration$duration_mins, decreasing=FALSE), ]
  print(agg_duration)
  ### Plot colours
  vec_algorithm = levels(df$algorithm)
  print(paste(vec_algorithm, collapse=" | "))
  vec_maf = sort(unique(df$maf))
  vec_missing_rate = sort(unique(df$missing_rate))
  vec_colours = c("#b2df8a", "#33a02c", "#1f78b4", "#a6cee3")
  vec_colours = rep(vec_colours, times=ceiling(length(vec_algorithm)/length(vec_colours)))[1:length(vec_algorithm)]
  #### Plot per dataset
  vec_fnames_svg = c()
  for (i in 1:length(vec_2_metrics)) {
    # i = 1
    metric = vec_2_metrics[i]
    metric_label = vec_2_metrics_labels[i]
    fname_svg = paste0(dataset, "-", gsub(" ", "_", metric_label), ".svg")
    vec_fnames_svg = c(vec_fnames_svg, fname_svg)
    svg(fname_svg, width=15, height=9)
    layout(matrix(c(
                    1, 1, 1, 2,
                    1, 1, 1, 2,
                    1, 1, 1, 2,
                    1, 1, 1, 2,
                    3, 3, 3, 4,
                    3, 3, 3, 4,
                    3, 3, 3, 4,
                    3, 3, 3, 4,
                    5, 6, 7, 8,
                    5, 6, 7, 8,
                    5, 6, 7, 8,
                    9,10,11,12,
                    9,10,11,12,
                    9,10,11,12
                   ), byrow=TRUE, ncol=4))
    par(mar=c(5,5,5,1))
    for (maf in vec_maf) {
      # maf = vec_maf[1]
      subdf = df[df$maf == maf, ]
      df_algo_miss = expand.grid(algorithm=vec_algorithm, missing_rate=vec_missing_rate)
      eval(parse(text=paste0("agg_mu = aggregate(", metric, " ~ algorithm + missing_rate, data=subdf, FUN=mean, na.rm=FALSE)")))
      eval(parse(text=paste0("agg_sd = aggregate(", metric, " ~ algorithm + missing_rate, data=subdf, FUN=function(x){3*sd(x,na.rm=TRUE)}); agg_sd$", metric, "[is.na(agg_sd$", metric, ")] = 0")))
      agg_mu = merge(df_algo_miss, agg_mu, by=c("algorithm", "missing_rate"), all=TRUE)
      agg_sd = merge(df_algo_miss, agg_sd, by=c("algorithm", "missing_rate"), all=TRUE)
      idx_sort = order(agg_mu$missing_rate)
      agg_mu = agg_mu[idx_sort, ]
      agg_sd = agg_sd[idx_sort, ]
      eval(parse(text=paste0("mat_mu = matrix(agg_mu$", metric, ", nrow=length(unique(agg_mu$algorithm)), byrow=FALSE); rownames(mat_mu) = agg_mu$algorithm[1:length(unique(agg_mu$algorithm))]; colnames(mat_mu) = round(sort(unique(agg_mu$missing_rate)), 4)")))
      eval(parse(text=paste0("mat_sd = matrix(agg_sd$", metric, ", nrow=length(unique(agg_sd$algorithm)), byrow=FALSE); rownames(mat_sd) = agg_sd$algorithm[1:length(unique(agg_sd$algorithm))]; colnames(mat_sd) = round(sort(unique(agg_sd$missing_rate)), 4)")))
      idx_sort = c()
      for (i in 1:length(vec_algorithm)) {
        idx_sort = c(idx_sort, which(rownames(mat_mu) == vec_algorithm[i]))
      }
      mat_mu = mat_mu[idx_sort, ]
      mat_sd = mat_sd[idx_sort, ]
      ### Barplot
      par(xpd=TRUE) ### xpd=TRUE allows us to place the legend outside the plot area
      if (grepl("concordance", metric) == TRUE) {
        y_lim = c(0, 1.00+max(c(0.01, max(agg_sd[,3], na.rm=TRUE))))
      } else {
        y_lim = c(0, max(mat_mu, na.rm=TRUE)+max(c(0.01, max(agg_sd[,3], na.rm=TRUE))))
      }
      bplot = barplot(mat_mu, beside=TRUE, col=vec_colours, ylim=y_lim, border=NA, main=paste0("maf = ", maf), xlab="Sparsity (missing/total)", ylab=metric_label, las=1)
      if (maf==min(vec_maf)) {
        legend("topleft", inset=c(-0.001, -0.5), legend=vec_algorithm, fill=vec_colours, bty="n", horiz=TRUE)
        legend("topright", inset=c(-0.001, -0.5), legend="Error bars refer to ± 3 standard deviations", bty="n", horiz=TRUE)
        par(xpd=FALSE)
      }
      arrows(bplot, mat_mu-mat_sd,
            bplot, mat_mu+mat_sd, length=0.05, angle=90, code=3)
      # text_pos = mat_mu+mat_sd
      # text_pos = text_pos + (0.1*text_pos)
      # if ((grepl("mae", metric) | grepl("rmse", metric)) == FALSE) {
      #   text_pos[text_pos < 0.1] = 0.1 ### So that the labels are still visible
      # }
      par(xpd=TRUE)
      text((bplot-0.25), 0.0, labels=round(mat_mu,4), cex=0.9, srt=90, pos=4)
      par(xpd=FALSE)
      ### Allele frequency distribution
      idx_af = which((vec_mu_af>=maf) & (vec_mu_af<=(1-maf)))
      vec_af = as.vector(mat_genotypes[, idx_af])
      vec_af = vec_af[!is.na(vec_af)]
      vec_af = c(vec_af, 1-vec_af)
      d = density(vec_af)
      plot(d, main=paste0("maf = ", maf), xlab="Allele frequency", ylab="Density", las=1, bty="n")
    }

    ### FOR MAF=5% ONLY:

    ### Plot mae across expected freqs at 5% MAF
    subdf = df[df$maf == 0.05, ]
    q = seq(0, 1, by=0.1)
    for (algo in vec_algorithm) {
      # algo = vec_algorithm[1]
      par(mar=c(5,5,1,6))
      idx_col_mae = which(grepl("^mae_", colnames(subdf)))
      plot(x=c(0,1), y=range(subdf[, idx_col_mae], na.rm=TRUE), type="n", xlab="Observed allele frequencies (5% MAF)", ylab="Imputation error\n(mean absolute error)", main=algo, las=1)
      colour = vec_colours[vec_algorithm==algo]
      for (i in 1:6) {
        # i = 1
        missing_rate = vec_missing_rate[i]
        # eval(parse(text=paste0("colour = rgb(", paste(as.vector(col2rgb(colour)), collapse = "/256,"), "/256, alpha=", missing_rate, ")")))
        idx = which((subdf$algorithm == algo) & (subdf$missing_rate == missing_rate))
        if (length(idx) > 1) {
          mae_across_freqs = NULL
          for (ix in idx) {
            # ix = idx[1]
            if (is.null(mae_across_freqs)) {
              eval(parse(text=paste0("mae_across_freqs = c(", paste(paste0("subdf$mae_", q, "[ix]"), collapse=", "), ")")))
            } else {
              eval(parse(text=paste0("mae_across_freqs = rbind(mae_across_freqs, c(", paste(paste0("subdf$`mae_", q, "`[ix]"), collapse=", "), "))")))
            }
          }
          mae_across_freqs = colMeans(mae_across_freqs, na.rm=TRUE)
        } else {
          eval(parse(text=paste0("mae_across_freqs = c(", paste(paste0("subdf$mae_", q, "[idx]"), collapse=", "), ")")))  
        }
        idx = which(!is.na(mae_across_freqs))

        if (dataset == "grape") {
          # ### The grape dataset only includes the minor alleles, hence to get the full picture we need to unfold the alleles to represent the reference alleles
          # lines(x=q[idx], y=(mae_across_freqs[idx] + rev(mae_across_freqs[idx]))/2, lwd=7*missing_rate, col=colour)
          ### The grape dataset uses the minor allele frequencies, hence to be uniform across all the datasets we used we take the additive inverse
          lines(x=1.00-q[idx], y=mae_across_freqs[idx], lty=i, col=colour)
        } else {
          ### The other datasets have all the alleles present and should expectedly generate a bell-ish curve or MAE across allele frequencies
          lines(x=q[idx], y=mae_across_freqs[idx], lty=i, col=colour)
        }
      }
      grid()
      if (algo == vec_algorithm[1]) {
        legend("topleft", legend=c("Sparsity", vec_missing_rate[1:6]), col=c(0, rep("black", 6)), lty=c(1,1:6), bty="n", cex=0.75)
      }
      par(mar=c(5, 4, 4, 2) + 0.1)
    }
    if (length(vec_algorithm) < 4) {
      plot(0, type="n", xlab="", ylab="", bty="n", xaxt="n", yaxt="n")
    }

    ### Plot allele frequency variances per locus x mae (MAF=5% only)
    subdf = df[df$maf == 0.05, ]
    for (algo in vec_algorithm) {
      # algo = vec_algorithm[1]
      colour = vec_colours[vec_algorithm==algo]
      idx_col = which(grepl("var_x_mae_", colnames(subdf)))
      n_bins = length(idx_col) / 2
      par(mar=c(5,5,1,6))
      idx_temp = which(subdf$algo == algo) 
      y_temp = data.frame(y=unlist(subdf[idx_temp, idx_col[seq(2, length(idx_col), by=2)]]), b=rep(1:n_bins, each=length(idx_temp)), m=rep(subdf$missing_rate[idx_temp], time=n_bins))
      y_temp = aggregate(y ~ b + m, data=y_temp, FUN=mean, na.rm=TRUE)[,3]
      plot(x=range(subdf[, idx_col[c(1, (length(idx_col)-1))]], na.rm=TRUE), y=range(y_temp, na.rm=TRUE), type="n", xlab="Observed variance in allele frequency per locus (5% MAF)", ylab="Imputation error\n(mean absolute error)", main=algo, las=1)
      for (i in 1:6) {
        # i = 1
        missing_rate = vec_missing_rate[i]
        # eval(parse(text=paste0("colour = rgb(", paste(as.vector(col2rgb(colour)), collapse = "/256,"), "/256, alpha=", missing_rate, ")")))
        idx_row = which((subdf$algo == algo) & (subdf$missing_rate == missing_rate))
        vec_bin = c()
        vec_mae = c()
        for (j in 1:n_bins) {
          # j = 1
          vec_bin = c(vec_bin, mean(subdf[idx_row, idx_col[2*(j-1)+1]]))
          vec_mae = c(vec_mae, mean(subdf[idx_row, idx_col[2*j]]))
        }
        lines(x=vec_bin, y=vec_mae, lty=i, col=colour)
      }
      grid()
      if (algo == vec_algorithm[1]) {
        legend("topleft", legend=c("Sparsity", vec_missing_rate[1:6]), col=c(0, rep("black", 6)), lty=c(1, 1:6), bty="n", cex=0.75)
      }
      par(mar=c(5, 4, 4, 2) +0.1)
    }
    dev.off()
  }
  return(vec_fnames_svg)
}

### Plot
vec_fnames = list.files(path=".", pattern="-performance_assessment-maf_")
vec_fnames = vec_fnames[grepl("-missing_rate_", vec_fnames)]
vec_fnames = vec_fnames[grepl(".csv$", vec_fnames)]
vec_datasets = unique(unlist(lapply(strsplit(vec_fnames, "-"), FUN=function(x){x[[1]]})))
for (dataset in vec_datasets) {
  # dataset = vec_datasets[1]
  vec_files = vec_fnames[grepl(paste0("^", dataset), vec_fnames)]
  for (i in 1:length(vec_files)) {
    df_tmp = read.csv(vec_files[i])
    if (i==1) {
      df = df_tmp
    } else {
      df = rbind(df, df_tmp)
    }
  }
  if ((dataset == "cocksfoot") | (dataset == "soybean")) {
    df = droplevels(df[df$algorithm != "linkimpute", ])
  }
  print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
  print(dataset)
  if (dataset == "grape") {
    vec_fnames_svg = plot_metrics(df=df, dataset=dataset, vec_2_metrics=c("mae_classes", "concordance_classes"))
  } else if (dataset == "cocksfoot") {
    vec_fnames_svg = plot_metrics(df=df, dataset=dataset)
    vec_fnames_svg_2 = plot_metrics(df=df, dataset=dataset, vec_2_metrics=c("highConf_mae_frequencies", "highConf_concordance_classes"), vec_2_metrics_labels=c("Mean absolute error high confidence data", "Concordance high confidence data"))
  } else {
    vec_fnames_svg = plot_metrics(df=df, dataset=dataset)
  }
  print(vec_fnames_svg)
  # if (dataset == "cocksfoot") {
  #   vec_fnames_svg = plot_metrics(df=df, dataset=dataset, vec_2_metrics=c("highConf_mae_frequencies", "highConf_concordance_classes"), vec_2_metrics_labels=c("Mean absolute error high confidence data", "Concordance high confidence data"))
  #   print(vec_fnames_svg)
  # }
}
