args = commandArgs(trailingOnly=TRUE)
# args = c("/group/pasture/Jeff/imputef/res")
dir_ouput = args[1]
setwd(dir_ouput)

### Plotting function per dataset
plot_metrics = function(df, dataset) {
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
  # vec_colours = c("#b2df8a", "#33a02c", "#a6cee3", "#1f78b4")
  # vec_colours = c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6")
  vec_colours = c("#88CCEE", "#44AA99", "#117733", "#CC6677", "#DDCC77", "#AA4499")
  # vec_colours = c("#88CCEE", "#CC6677", "#44AA99", "#DDCC77", "#AA4499", "#117733")
  # vec_colours = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")
  # vec_colours = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
  # vec_colours = c("#00bf7d", "#00b4c5", "#0073e6", "#2546f0", "#5928ed", "#FFC20A")
  # vec_colours = c("#b3c7f7", "#8babf1", "#0073e6", "#0461cf", "#054fb9", "#0C7BDC")
  # vec_colours = c("#c44601", "#f57600", "#8babf1", "#0073e6", "#054fb9", "#994F00")
  # vec_colours = c("#38A3A5", "#BDB76B", "#F4D03F", "#89CEDB", "#E0293F", "#333333")
  # vec_colours = c("#0089BF", "#00B4C5", "#D2B48C", "#F06292", "#A29BDB", "#2546F0")
  # vec_colours = c("#FDB400", "#E5603D", "#A4C139", "#00B2DF", "#9B97D3", "#62368F")
  # vec_colours = c("#91D2BD", "#FBE9D1", "#F5E6CE", "#D8C092", "#A28F9D", "#E1F5FE")
  # vec_colours = c("#FFC0A8", "#D3AEB4", "#E5C171", "#9BD6B0", "#789DA7", "#4D4D4F")
  # vec_colours = c("#00B894", "#F27032", "#F06292", "#A4C139", "#A29BDB", "#001F3F")
  # vec_colours = c("#001B43", "#13CFE9", "#A9D04F", "#F8B195", "#D3D3D3", "#2E343A")
  vec_colours = rep(vec_colours, times=ceiling(length(vec_algorithm)/length(vec_colours)))[1:length(vec_algorithm)]
  vec_metrics = c("mae_frequencies", "concordance_classes")
  vec_metrics_labels = c("Mean absolute error", "Concordance")
  vec_fnames_svg = c()
  # n_plots = 2*length(vec_maf)
  for (i in 1:length(vec_metrics)) {
    # i = 2
    metric = vec_metrics[i]
    metric_label = vec_metrics_labels[i]
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
      text((bplot-0.25), text_pos, labels=round(mat_mu,2), cex=0.9, srt=90, pos=4)
      par(xpd=FALSE)
      ### Line plot
      plot(x=range(agg_mu[,2], na.rm=TRUE), y=range(agg_mu[,3]+(2*agg_sd[,3]), na.rm=TRUE), 
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


    ### Tests: plotting mae across freqs
    q = c(0.0, 0.01, 0.05, seq(0.1, 0.9, by=0.1), 0.95, 0.99, 1.00)
    for (algo in vec_algorithm) {
      par(mar=c(5,5,1,1))
      plot(x=c(0,1), y=c(0,1), type="n", xlab="Expected allele frequencies", ylab="Imputation error\n(mean absolute error)", main=algo, las=1)
      colour = vec_colours[vec_algorithm==algo]
      for (missing_rate in vec_missing_rate) {
        eval(parse(text=paste0("colour = rgb(", paste(as.vector(col2rgb(colour)), collapse = "/256,"), "/256, alpha=", missing_rate, ")")))
        idx = which((subdf$algorithm == algo) & (subdf$missing_rate == missing_rate))
        eval(parse(text=paste0("mae_across_freqs = c(", paste(paste0("subdf$mae_", q, "[idx]"), collapse=", "), ")")))
        idx = which(!is.na(mae_across_freqs))
        lines(x=q[idx], y=mae_across_freqs[idx], lwd=5*missing_rate, col=colour)
      }
      grid()
      if (algo == vec_algorithm[1]) {
        legend("topleft", legend=c("Sparsity", vec_missing_rate[1:5]), col=c(0, rgb(0,0,0,alpha=vec_missing_rate[1:5])), lwd=c(0, 5*vec_missing_rate[1:5]), bty="n")
        legend("top", legend=c("", vec_missing_rate[6:10]), col=c(0, rgb(0,0,0,alpha=vec_missing_rate[6:10])), lwd=c(0, 5*vec_missing_rate[6:10]), bty="n")
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
  vec_fnames_svg = plot_metrics(df=df, dataset=dataset)
  print(vec_fnames_svg)
}
