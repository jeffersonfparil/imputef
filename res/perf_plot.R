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
  # Sort algorithms according to increasing complexity and with LinkImpute at the bottom as it will not be assessed for polyploid and pool datasets
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

  # agg_concordance = aggregate(concordance_classes ~ algorithm + maf + missing_rate, data=df, FUN=mean)
  agg_concordance = aggregate(concordance_classes ~ algorithm, data=df, FUN=mean)
  agg_concordance = agg_concordance[order(agg_concordance$concordance_classes, decreasing=TRUE), ]
  print(agg_concordance)

  # agg_duration = aggregate(duration_mins ~ algorithm + maf + missing_rate, data=df, FUN=mean)
  agg_duration = aggregate(duration_mins ~ algorithm, data=df, FUN=mean)
  agg_duration = agg_duration[order(agg_duration$duration_mins, decreasing=FALSE), ]
  print(agg_duration)

  vec_algorithm = levels(df$algorithm)
  print(paste(vec_algorithm, collapse=" | "))
  vec_maf = sort(unique(df$maf))
  vec_missing_rate = sort(unique(df$missing_rate))
  vec_colours = c("#b2df8a", "#33a02c", "#1f78b4", "#a6cee3")
  vec_colours = rep(vec_colours, times=ceiling(length(vec_algorithm)/length(vec_colours)))[1:length(vec_algorithm)]
  
  vec_fnames_svg = c()
  for (i in 1:length(vec_2_metrics)) {
    # i = 1
    metric = vec_2_metrics[i]
    metric_label = vec_2_metrics_labels[i]
    if (grepl("mae", metric) | grepl("rmse", metric)) {
      eval(parse(text=paste0("ylim = c(0, max(df$", metric, ", na.rm=TRUE)+(2*sd(df$", metric, ", na.rm=TRUE)))")))
    } else {
      ylim = c(0, 1.17)
    }
    fname_svg = paste0(dataset, "-", gsub(" ", "_", metric_label), ".svg")
    vec_fnames_svg = c(vec_fnames_svg, fname_svg)
    svg(fname_svg, width=15, height=9)
    n_cols = 48
    layout(matrix(c(rep(1, n_cols*(5/8)), rep(2, n_cols*(3/8)),
                    rep(3, n_cols*(5/8)), rep(4, n_cols*(3/8)),
                    rep(5:8, each=n_cols/4),
                    rep(9:12, each=n_cols/4)
                   ), byrow=TRUE, ncol=n_cols))
    for (maf in vec_maf) {
      # maf = vec_maf[1]
      subdf = df[df$maf == maf, ]
      df_algo_miss = expand.grid(algorithm=vec_algorithm, missing_rate=vec_missing_rate)
      eval(parse(text=paste0("agg_mu = aggregate(", metric, " ~ algorithm + missing_rate, data=subdf, FUN=mean, na.rm=FALSE)")))
      eval(parse(text=paste0("agg_sd = aggregate(", metric, " ~ algorithm + missing_rate, data=subdf, FUN=sd, na.rm=TRUE); agg_sd$", metric, "[is.na(agg_sd$", metric, ")] = 0")))
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
      bplot = barplot(mat_mu, beside=TRUE, col=vec_colours, border=NA, ylim=ylim, main=paste0("maf = ", maf), xlab="Sparsity (missing/total)", ylab=metric_label, las=1)
      if (maf==min(vec_maf)) {
        legend("topright", inset=c(-0.001, -0.5), legend=vec_algorithm, fill=vec_colours, bty="n", horiz=TRUE)
        par(xpd=FALSE)
      }
      arrows(bplot, mat_mu-mat_sd,
            bplot, mat_mu+mat_sd, length=0.05, angle=90, code=3)
      grid()
      text_pos = mat_mu+mat_sd
      text_pos = text_pos + (0.1*text_pos)
      if ((grepl("mae", metric) | grepl("rmse", metric)) == FALSE) {
        text_pos[text_pos < 0.1] = 0.1 ### So that the labels are still visible
      }
      par(xpd=TRUE)
      text((bplot-0.25), text_pos, labels=round(mat_mu,4), cex=0.9, srt=90, pos=4)
      par(xpd=FALSE)
      ### Line plot
      plot(x=range(agg_mu[,2], na.rm=TRUE), y=range(c(agg_mu[,3]-(2*agg_sd[,3]), agg_mu[,3]+(2*agg_sd[,3])), na.rm=TRUE), 
        main=paste0("maf = ", maf), xlab="Sparsity (missing/total)", ylab=metric_label, las=1, type="n")
      grid()
      for (algo in vec_algorithm) {
        # algo = vec_algorithm[1]
        idx = agg_mu$algorithm == algo
        agg_mu_sub = agg_mu[idx, ]
        agg_sd_sub = agg_sd[idx, ]
        lines(x=agg_mu_sub[,2], y=agg_mu_sub[,3], col=vec_colours[vec_algorithm==algo], lwd=2)
        points(x=agg_mu_sub[,2], y=agg_mu_sub[,3], col="black", bg=vec_colours[vec_algorithm==algo], pch=23)
        arrows(x0=agg_mu_sub[,2], y0=agg_mu_sub[,3]-agg_sd_sub[,3],
              x1=agg_mu_sub[,2], y1=agg_mu_sub[,3]+agg_sd_sub[,3], length=0.05, angle=90, code=3)
      }
    }
    ### Plot mae across expected freqs at 5% MAF
    subdf = df[df$maf == 0.05, ]
    q = seq(0, 1, by=0.1)
    for (algo in vec_algorithm) {
      # algo = vec_algorithm[1]
      par(mar=c(5,5,1,1))
      idx_col_mae = which(grepl("mae_", colnames(subdf)))
      plot(x=c(0,1), y=range(subdf[, idx_col_mae], na.rm=TRUE), type="n", xlab="Observed allele frequencies (5% MAF)", ylab="Imputation error\n(mean absolute error)", main=algo, las=1)
      colour = vec_colours[vec_algorithm==algo]
      for (missing_rate in vec_missing_rate) {
        # missing_rate = vec_missing_rate[1]
        eval(parse(text=paste0("colour = rgb(", paste(as.vector(col2rgb(colour)), collapse = "/256,"), "/256, alpha=", missing_rate, ")")))
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
          lines(x=1.00-q[idx], y=mae_across_freqs[idx], lwd=7*missing_rate, col=colour)
        } else {
          ### The other datasets have all the alleles present and should expectedly generate a bell-ish curve or MAE across allele frequencies
          lines(x=q[idx], y=mae_across_freqs[idx], lwd=7*missing_rate, col=colour)
        }
      }
      grid()
      if (algo == vec_algorithm[1]) {
        n0 = 1
        n2 = length(unique(subdf$missing_rate))
        n1 = floor(n2/2)
        legend("topleft", legend=c("Sparsity", vec_missing_rate[n0:n1]), col=c(0, rgb(0,0,0,alpha=vec_missing_rate[n0:n1])), lwd=c(0, 5*vec_missing_rate[n0:n1]), bty="n")
        legend("top", legend=c("", vec_missing_rate[(n1+1):n2]), col=c(0, rgb(0,0,0,alpha=vec_missing_rate[(n1+1):n2])), lwd=c(0, 5*vec_missing_rate[(n1+1):n2]), bty="n")
      }
      par(mar=c(5, 4, 4, 2) +0.1)
    }

    ### Plot allele frequency variances per locus x mae (MAF=5% only)
    subdf = df[df$maf == 0.05, ]
    for (algo in vec_algorithm) {
      # algo = vec_algorithm[1]
      colour = vec_colours[vec_algorithm==algo]
      idx_col = which(grepl("var_x_mae_", colnames(subdf)))
      n_bins = length(idx_col) / 2
      par(mar=c(5,5,1,1))
      idx_temp = which(subdf$algo == algo) 
      y_temp = data.frame(y=unlist(subdf[idx_temp, idx_col[seq(2, length(idx_col), by=2)]]), b=rep(1:n_bins, each=length(idx_temp)), m=rep(subdf$missing_rate[idx_temp], time=n_bins))
      y_temp = aggregate(y ~ b + m, data=y_temp, FUN=mean, na.rm=TRUE)[,3]
      plot(x=range(subdf[, idx_col[c(1, (length(idx_col)-1))]], na.rm=TRUE), y=range(y_temp, na.rm=TRUE), type="n", xlab="Observed variance in allele frequency per locus (5% MAF)", ylab="Imputation error\n(mean absolute error)", main=algo, las=1)
      for (missing_rate in vec_missing_rate) {
        # missing_rate = vec_missing_rate[1]
        eval(parse(text=paste0("colour = rgb(", paste(as.vector(col2rgb(colour)), collapse = "/256,"), "/256, alpha=", missing_rate, ")")))
        idx_row = which((subdf$algo == algo) & (subdf$missing_rate == missing_rate))
        vec_bin = c()
        vec_mae = c()
        for (j in 1:n_bins) {
          # j = 1
          vec_bin = c(vec_bin, mean(subdf[idx_row, idx_col[2*(j-1)+1]]))
          vec_mae = c(vec_mae, mean(subdf[idx_row, idx_col[2*j]]))
        }
        lines(x=vec_bin, y=vec_mae, lwd=7*missing_rate, col=colour)
      }
      grid()
      if (algo == vec_algorithm[1]) {
        n0 = 1
        n2 = length(unique(subdf$missing_rate))
        n1 = floor(n2/2)
        legend("topleft", legend=c("Sparsity", vec_missing_rate[n0:n1]), col=c(0, rgb(0,0,0,alpha=vec_missing_rate[n0:n1])), lwd=c(0, 5*vec_missing_rate[n0:n1]), bty="n")
        legend("top", legend=c("", vec_missing_rate[(n1+1):n2]), col=c(0, rgb(0,0,0,alpha=vec_missing_rate[(n1+1):n2])), lwd=c(0, 5*vec_missing_rate[(n1+1):n2]), bty="n")
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
  if ((dataset == "lucerne") | (dataset == "soybean")) {
    df = droplevels(df[df$algorithm != "linkimpute", ])
  }
  print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
  print(dataset)
  if (dataset == "grape") {
    vec_fnames_svg = plot_metrics(df=df, dataset=dataset, vec_2_metrics=c("mae_classes", "concordance_classes"))
  } else if (dataset == "lucerne") {
    vec_fnames_svg = plot_metrics(df=df, dataset=dataset)
    vec_fnames_svg_2 = plot_metrics(df=df, dataset=dataset, vec_2_metrics=c("highConf_mae_frequencies", "highConf_concordance_classes"), vec_2_metrics_labels=c("Mean absolute error high confidence data", "Concordance high confidence data"))
  } else {
    vec_fnames_svg = plot_metrics(df=df, dataset=dataset)
  }
  print(vec_fnames_svg)
  # if (dataset == "lucerne") {
  #   vec_fnames_svg = plot_metrics(df=df, dataset=dataset, vec_2_metrics=c("highConf_mae_frequencies", "highConf_concordance_classes"), vec_2_metrics_labels=c("Mean absolute error high confidence data", "Concordance high confidence data"))
  #   print(vec_fnames_svg)
  # }
}
