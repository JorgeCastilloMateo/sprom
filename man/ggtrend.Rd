% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggtrend.R
\name{ggtrend}
\alias{ggtrend}
\title{Plot trend}
\usage{
ggtrend(vol3D, yint, ychar, xchar, xnumbreaks, xlabbreaks)
}
\arguments{
\item{vol3D}{Input 3D volume.}

\item{yint}{y-axis limits.}

\item{ychar}{y-axis label.}

\item{xchar}{x-axis label.}

\item{xnumbreaks}{x-axis breaks.}

\item{xlabbreaks}{x-axis label breaks.}
}
\value{
A ggplot2 object depicting the t-times weighted trend.
}
\description{
This function generates a ggplot graph with a weighted average
of indicators with an stationary decay of 1/t. Being t (time) the x-axis.
}
\details{
The weighted trend will be summarized over the first dimension of the input 3D volume.
}
\examples{

# load("data/Y365.rda")
# ggtrend(vol3D = Y365,
#       yint = c(0, 3),
#       ychar = expression(t \%*\% hat(p)[t]),
#       xchar = "t (year)",
#       xnumbreaks = c(1, 21, 41, 61),
#       xlabbreaks = c("1 (1960)", "21 (1980)", "41 (2000)", "61 (2020)"))

}
\author{
Zeus Gracia-Tabuenca
}
