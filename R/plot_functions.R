#' @title Writes out a violin plot of -log10(p) values for quantiles of q
#'
#' @description Can be used to investigate the relationship between p and q
#'
#' @param p p values for principal trait (vector of length n)
#' @param q auxiliary data values (vector of length n)
#' @param ylim y-axis limits (-log10)
#'
#' @return ggplot object
#' @import ggplot2
#' @import cowplot
#' @export
corr_plot <- function(p, q, ylim = c(0, 1.5)){
  
  quantiles <- quantile(q, probs = seq(0, 1, 0.2))
  
  quans <- cut(q, c(unname(quantiles)))
  
  quans = factor(quans,levels(quans))
  
  w1 <- which(quans==levels(quans)[[1]])
  w2 <- which(quans==levels(quans)[[2]])
  w3 <- which(quans==levels(quans)[[3]])
  w4 <- which(quans==levels(quans)[[4]])
  w5 <- which(quans==levels(quans)[[5]])
  
  n1 <- length(w1)
  n2 <- length(w2)
  n3 <- length(w3)
  n4 <- length(w4)
  n5 <- length(w5)
  
  df <- data.frame(p = -log10(c(p[w1], p[w2], p[w3], p[w4], p[w5])), q = c(q[w1], q[w2], q[w3], q[w4], q[w5]), quantiles = c(rep(levels(quans)[[1]], n1), rep(levels(quans)[[2]], n2), rep(levels(quans)[[3]], n3), rep(levels(quans)[[4]], n4), rep(levels(quans)[[5]], n5)))
  
  df$quantiles = factor(df$quantiles, levels(quans))
  
  ggplot(df, aes(x = quantiles, y = p, group = quantiles)) + geom_violin(aes(fill = quantiles))+
    theme_cowplot(12) + background_grid(major = "xy", minor = "xy") + xlab("q") + ylab("p (-log10)") + theme(legend.text=element_text(size=8)) + geom_boxplot(aes(fill=quantiles), width = 0.1) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_cartesian(ylim = ylim)
}

#' @title Writes out a plot of p against v coloured by q
#'
#' @description Can be used to visualise the results from Flexible cFDR
#' 
#' @param p p values for principal trait (vector of length n)
#' @param q auxiliary data values (vector of length n)
#' @param v v values from cFDR
#' @param axis_lim Optional axis limits
#'
#' @return ggplot object
#' @import ggplot2
#' @import cowplot
#' @export
pv_plot <- function(p, q, v, axis_lim = c(0, 1)){

  mid <- median(q)
  
  df <- data.frame(p, q, v)

  ggplot(df, aes(x = p, y = v, col = q)) + geom_point() + theme_cowplot(12) + background_grid(major = "xy", minor = "none") + geom_abline(intercept = 0, slope = 1,  linetype="dashed") + xlab("raw p-values") + ylab("v-values") + ggtitle(paste0("Flexible cFDR results")) + scale_color_gradient2(midpoint = mid, low = "blue", mid = "white", high = "red", space = "Lab") + coord_cartesian(ylim = axis_lim, xlim = axis_lim)
}

#' @title Writes out a plot of -log10(p) against -log10(v) coloured by q
#'
#' @description Can be used to visualise the results from Flexible cFDR
#'
#' @param p p values for principal trait (vector of length n)
#' @param q auxiliary data values (vector of length n)
#' @param v v values from cFDR
#' @param axis_lim Optional axis limits
#'
#' @return ggplot object
#' @import ggplot2
#' @import cowplot
#' @export
log10pv_plot <- function(p, q, v, axis_lim = c(0, 20)){

  mid <- median(q)
  
  df <- data.frame(p, q, v)

  ggplot(df, aes(x = -log10(p), y = -log10(v), col = q)) + geom_point() + theme_cowplot(12) + background_grid(major = "xy", minor = "none") + geom_abline(intercept = 0, slope = 1,  linetype="dashed") + xlab("raw p-values (-log10)") + ylab("v-values (-log10)") + ggtitle(paste0("Flexible cFDR results")) + scale_color_gradient2(midpoint = mid, low = "blue", mid = "white", high = "red", space = "Lab") + coord_cartesian(ylim = axis_lim, xlim = axis_lim)
}

#' @title Writes out a stratified Q-Q plot.
#'
#' @description Can be used to investigate the relationship between p and q
#'
#' @details Note that this function does not do the heavy lifting of styling the plot's aesthetics.
#' 
#' @param data_frame \code{data.frame} containing p-values and auxiliary data values
#' @param prin_value_label label of principal p-value column in \code{data_frame}
#' @param cond_value_label label of conditional trait column in \code{data_frame}
#' @param thresholds threshold values to define strata
#'
#' @import ggplot2 
#' @import cowplot
#'
#' @return ggplot object 
#' @export
#' 
stratified_qqplot <- function(data_frame, prin_value_label, cond_value_label = NULL, thresholds = c(1, 1e-1, 1e-2, 1e-3, 1e-4)) {
  
  data_frame$negLogP <- -log10(data_frame[, prin_value_label])
  
  if(is.null(cond_value_label)) {
    daf <- data_frame[, c(prin_value_label, 'negLogP')]
    daf <- daf[order(daf[,prin_value_label]), ]
    daf$pp <- -log10(ppoints(nrow(daf)))
    daf$threshold <- factor(c(1))
  } else {
    data_frame <- data_frame[, c(prin_value_label, cond_value_label, 'negLogP')]
    
    dafs <- list()
    
    for(i in seq_along(thresholds)) {
      daf <- subset(data_frame, get(cond_value_label) < thresholds[i])
      daf <- daf[order(daf[ , prin_value_label]) , ]
      daf$pp <- -log10(ppoints(nrow(daf)))
      daf$threshold <- factor(thresholds)[i]
      dafs[[i]] <- daf
    }
    
    daf <- do.call(rbind, dafs)
  }
  
  ggplot(data=daf) + geom_line(aes(x = pp, y = negLogP, group = threshold, colour = threshold)) + geom_abline(intercept=0,slope=1, linetype="dashed") + theme_cowplot(12) + background_grid(major = "xy", minor = "none")

}
