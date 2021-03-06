% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ANCOM.R
\name{ANCOM}
\alias{ANCOM}
\title{ANCOM, the traditional method for analysis of composition of microbiomes}
\usage{
ANCOM(
  Y,
  x,
  tfun = "t2",
  pseudo = 0.5,
  sig = 0.05,
  zero_cut = 0.9,
  lib_cut = 0,
  struc_zero = F
)
}
\arguments{
\item{Y}{the abundance table,  count/relative or any smoothed data matrix are supported}

\item{x}{the variable of interest}

\item{tfun}{the test function used, including t.test, wilcox.test, kruskal.test, and a speed version for wilcox.test based on the Mann-Whitney U statistic which denote as 't2'}

\item{pseudo}{the pseudo number for zero smoothing, default is 0.5}

\item{sig}{the level of significance to count the reject number, default is 0.05}

\item{zero_cut}{a numerical fraction between 0 and 1. Taxa with proportion of zeroes greater than zero_cut will be excluded in the analysis. Default is 0.90}

\item{lib_cut}{a numerical threshold for filtering samples based on library sizes. Samples with library sizes less than lib_cut will be excluded in the analysis. Default is 0, i.e. do not filter any sample.}

\item{struc_zero}{whether to detect structural zeros. Default is FALSE}
}
\value{
a dataframe with 3 columns,
\itemize{
\item Reject.number, the reject number of each microbe
\item REJECT, the final result whether each microbe is associated with the variable of interest or not
\item detect0.7 hard threshold 0.7 to identify the differential microbes
}
}
\description{
The traditional one uses non-parameteric methods, Wilcoxon or Kruskal-Wallis rank-sum test for two or more groups, to test each joint model directly.
And the emperical distribution of the reject numbers are used to define a cutoff for concluding each microbe to be differential or not. This is same function which implements in QIIME2 with python.
If there is other covariates, as implement in https://github.com/FrederickHuangLin/ANCOM who used ANOVA. But it could be used fastANCOM to achieve the same goal with the log-linear model.
}
\examples{
data <- matrix(rpois(100*60, 10), 60)
b <- rnorm(10, 4)
data[31:60,1:10] <- t(t(data[31:60,1:10])*2^b)
group <- rep(c(0, 1), each=30)
fit <- ANCOM(Y=data, x=group)
head(fit)
}
