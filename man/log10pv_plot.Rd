% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{log10pv_plot}
\alias{log10pv_plot}
\title{Plot -log10(p) against -log10(v) and colour by q}
\usage{
log10pv_plot(p, q, v, axis_lim = c(0, 20))
}
\arguments{
\item{p}{p values for principal trait (vector of length n)}

\item{q}{auxiliary data values (vector of length n)}

\item{v}{v values from cFDR}

\item{axis_lim}{Optional axis limits}
}
\value{
ggplot object
}
\description{
Plot -log10(p) against -log10(v) and colour by q
}
\details{
Can be used to visualise the results from Flexible cFDR
}
