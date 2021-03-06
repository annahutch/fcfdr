% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{stratified_qqplot}
\alias{stratified_qqplot}
\title{Stratified Q-Q plot.}
\usage{
stratified_qqplot(
  data_frame,
  prin_value_label,
  cond_value_label = NULL,
  thresholds = c(1, 0.1, 0.01, 0.001, 1e-04)
)
}
\arguments{
\item{data_frame}{\code{data.frame} containing p-values and auxiliary data values}

\item{prin_value_label}{label of principal p-value column in \code{data_frame}}

\item{cond_value_label}{label of conditional trait column in \code{data_frame}}

\item{thresholds}{threshold values to define strata}
}
\value{
ggplot object
}
\description{
Stratified Q-Q plot.
}
\details{
Can be used to investigate the relationship between p and q

Note that this function does not do the heavy lifting of styling the plot's aesthetics.
}
