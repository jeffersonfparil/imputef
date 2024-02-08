args = commandArgs(trailingOnly=TRUE)
# args = c("/group/pasture/Jeff/imputef/res")
dir_ouput = args[1]
setwd(dir_ouput)

### Plotting function per dataset
plot_metrics = function(df, dataset, vec_2_metrics=c("mae_frequencies", "concordance_classes"), vec_2_metrics_labels=c("Mean absolute error", "Concordance")) {
  # vec_2_metrics = c("mae_frequencies", "concordance_classes")
  # vec_2_metrics_labels = c("Mean absolute error", "Concordance")
  # Sort algorithms according to increasing complexity and with LinkImpute at the bottom as it will not be assessed for polyploid and pool datasets
  df$algorithm = as.character(df$algorithm)
  df$algorithm[df$algorithm=="mvi"] = "a"
  df$algorithm[df$algorithm=="aldknni_fixed"] = "b"
  df$algorithm[df$algorithm=="aldknni_optim_cd"] = "c"
  df$algorithm[df$algorithm=="aldknni_optim_lk"] = "d"
  df$algorithm[df$algorithm=="aldknni_optim_all"] = "e"
  df$algorithm[df$algorithm=="linkimpute"] = "f"
  df$algorithm = as.factor(df$algorithm)
  levels(df$algorithm)[levels(df$algorithm)=="a"] = "MVI"
  levels(df$algorithm)[levels(df$algorithm)=="b"] = "AFIXED"
  levels(df$algorithm)[levels(df$algorithm)=="c"] = "AOPTIMCD"
  levels(df$algorithm)[levels(df$algorithm)=="d"] = "AOPTIMLK"
  levels(df$algorithm)[levels(df$algorithm)=="e"] = "AOPTIMAL"
  levels(df$algorithm)[levels(df$algorithm)=="f"] = "LINKIMPUTE"

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
  vec_colours = c("#88CCEE", "#44AA99", "#117733", "#CC6677", "#DDCC77", "#AA4499")
  vec_colours = rep(vec_colours, times=ceiling(length(vec_algorithm)/length(vec_colours)))[1:length(vec_algorithm)]
  
  vec_fnames_svg = c()
  # n_plots = 2*length(vec_maf)
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
    svg(fname_svg, width=15, height=7)
    # mat_layout_base = matrix(1:n_plots, nrow=(n_plots/2), byrow=TRUE)
    # layout(cbind(mat_layout_base[,1], mat_layout_base[,1], mat_layout_base))
    layout(matrix(c(1, 1, 1, 1, 2,  2,
                    3, 3, 3, 3, 4,  4,
                    5, 6, 7, 8, 9, 10), byrow=TRUE, ncol=6))
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
      plot(x=c(0,1), y=c(0,1), type="n", xlab="Expected allele frequencies (5% MAF)", ylab="Imputation error\n(mean absolute error)", main=algo, las=1)
      colour = vec_colours[vec_algorithm==algo]
      for (missing_rate in vec_missing_rate) {
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
          mae_across_freqs = colMeans(mae_across_freqs)
        } else {
          eval(parse(text=paste0("mae_across_freqs = c(", paste(paste0("subdf$mae_", q, "[idx]"), collapse=", "), ")")))  
        }
        idx = which(!is.na(mae_across_freqs))
        lines(x=q[idx], y=mae_across_freqs[idx], lwd=5*missing_rate, col=colour)
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
  # dataset = vec_datasets[3]
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
  vec_fnames_svg = plot_metrics(df=df, dataset=dataset)
  print(vec_fnames_svg)
  if (dataset == "lucerne") {
    vec_fnames_svg = plot_metrics(df=df, dataset=dataset, vec_2_metrics=c("highConf_mae_frequencies", "highConf_concordance_classes"), vec_2_metrics_labels=c("Mean absolute error high confidence data", "Concordance high confidence data"))
    print(vec_fnames_svg)
  }
}
