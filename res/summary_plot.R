args = commandArgs(trailingOnly=TRUE)
# args = c("/group/pasture/Jeff/imputef/res")
dir_ouput = args[1]
setwd(dir_ouput)

### Plotting function per dataset
plot_metrics = function(df, dataset) {
  df$algorithm[df$algorithm=="mvi"] = "MVI"
  # df$algorithm[df$algorithm=="lukes"] = "SAMP"
  df$algorithm[df$algorithm=="aldknni"] = "ALDKNNI"
  df$algorithm[df$algorithm=="aldknni_optim_thresholds"] = "ALDKNNI_OPTIM_THRESHOLDS"
  df$algorithm[df$algorithm=="aldknni_optim_counts"] = "ALDKNNI_OPTIM_COUNTS"
  df$algorithm[df$algorithm=="linkimpute"] = "LINKIMPUTE"

  # agg_concordance = aggregate(concordance_classes ~ algorithm + maf + missing_rate, data=df, FUN=mean)
  agg_concordance = aggregate(concordance_classes ~ algorithm, data=df, FUN=mean)
  agg_concordance = agg_concordance[order(agg_concordance$concordance_classes, decreasing=TRUE), ]
  print(agg_concordance)

  # agg_duration = aggregate(duration_mins ~ algorithm + maf + missing_rate, data=df, FUN=mean)
  agg_duration = aggregate(duration_mins ~ algorithm, data=df, FUN=mean)
  agg_duration = agg_duration[order(agg_duration$duration_mins, decreasing=FALSE), ]
  print(agg_duration)

  vec_algorithm = sort(unique(df$algorithm))
  print(paste(vec_algorithm, collapse=" | "))
  vec_maf = sort(unique(df$maf))
  # vec_colours = c("#b2df8a", "#33a02c", "#a6cee3", "#1f78b4")
  # vec_colours = c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6")
  # vec_colours = c("#88CCEE", "#44AA99", "#117733", "#CC6677", "#DDCC77", "#AA4499")
  # vec_colours = c("#88CCEE", "#CC6677", "#44AA99", "#DDCC77", "#AA4499", "#117733")
  vec_colours = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")
  if (length(vec_algorithm)==4) {
    vec_colours = vec_colours[c(1, 2, 3, 5)]
  }
  vec_colours = rep(vec_colours, times=ceiling(length(vec_algorithm)/length(vec_colours)))[1:length(vec_algorithm)]
  vec_metrics = c("mae_freqs", "r2_freqs", "concordance_classes")
  vec_metrics_labels = c("Mean absolute error", "Coefficient of determination", "Concordance")
  vec_fnames_svg = c()
  n_plots = 2*length(vec_maf)
  for (i in 1:length(vec_metrics)) {
    # i = 1
    metric = vec_metrics[i]
    metric_label = vec_metrics_labels[i]
    if (grepl("mae", metric) | grepl("rmse", metric)) {
      eval(parse(text=paste0("ylim = c(0, max(df$", metric, ", na.rm=TRUE)+(2*sd(df$", metric, ", na.rm=TRUE)))")))
    } else {
      ylim = c(0, 1.1)
    }
    fname_svg = paste0(dataset, "-", gsub(" ", "_", metric_label), ".svg")
    vec_fnames_svg = c(vec_fnames_svg, fname_svg)
    svg(fname_svg, width=12, height=9)
    mat_layout_base = matrix(1:n_plots, nrow=(n_plots/2), byrow=TRUE)
    layout(cbind(mat_layout_base[,1], mat_layout_base[,1], mat_layout_base))
    for (maf in vec_maf) {
      # maf = vec_maf[1]
      subdf = df[df$maf == maf, ]
      eval(parse(text=paste0("agg_mu = aggregate(", metric, " ~ algorithm + missing_rate, data=subdf, FUN=mean)")))
      eval(parse(text=paste0("agg_sd = aggregate(", metric, " ~ algorithm + missing_rate, data=subdf, FUN=sd, na.rm=TRUE); agg_sd$", metric, "[is.na(agg_sd$", metric, ")] = 0")))
      eval(parse(text=paste0("mat_mu = matrix(agg_mu$", metric, ", nrow=length(unique(agg_mu$algorithm)), byrow=FALSE); rownames(mat_mu) = agg_mu$algorithm[1:length(unique(agg_mu$algorithm))]; colnames(mat_mu) = round(sort(unique(agg_mu$missing_rate)), 4)")))
      eval(parse(text=paste0("mat_sd = matrix(agg_sd$", metric, ", nrow=length(unique(agg_sd$algorithm)), byrow=FALSE); rownames(mat_sd) = agg_sd$algorithm[1:length(unique(agg_sd$algorithm))]; colnames(mat_sd) = round(sort(unique(agg_sd$missing_rate)), 4)")))
      ### Barplot
      par(xpd=TRUE) ### xpd=TRUE allows us to place the legend outside the plot area
      bplot = barplot(mat_mu, beside=TRUE, col=vec_colours, border=NA, ylim=ylim, main=paste0("maf = ", maf), xlab="Sparsity (missing/total)", ylab=metric_label, las=1)
      if (maf==min(vec_maf)) {
        legend("topright", inset=c(-0.001, -0.5), legend=vec_algorithm, fill=vec_colours, bty="n")
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
