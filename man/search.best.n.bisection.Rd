% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/powercalc.R
\name{search.best.n.bisection}
\alias{search.best.n.bisection}
\title{Method of bisection for estimating optimal sample size}
\usage{
search.best.n.bisection(
  R,
  betaSE,
  lambda,
  pos,
  nstart,
  nmax,
  weights = 1,
  typeII = 0.2,
  alpha = 0.01
)
}
\arguments{
\item{R}{(p x p) matrix containing theoretical correlation between SNP pairs}

\item{betaSE}{effect size relative to residual standard deviation}

\item{lambda}{shrinkage parameter}

\item{pos}{vector (LEN nqtl) of SNP indices for assumed QTL positions}

\item{nstart}{minimum value for grid search}

\item{nmax}{maximum value for grid search}

\item{weights}{vector (LEN p) of SNP-specific weights or scalar if weights
are equal for all SNPs; default value 1}

\item{typeII}{type-II error level; default value 0.2}

\item{alpha}{type-I error level; default value 0.01}
}
\value{
integer of optimal sample size
}
\description{
A grid \code{[nstart, nmax]} for possible sample size is
  considered. Instead of executing a time-consuming grid search, the method
  of bisection is applied to this interval. For each step, the function
  \code{pwr.normtest} is called for the given set of parameters.
}
\examples{
  ### correlation matrix (should depend on sire haplotypes)
  R <- AR1(100, rho = 0.1)
  ### positions of putative QTL signals
  pos <- c(14, 75)
  ### optimal sample size
  search.best.n.bisection(R, 0.35, 1200, pos, 10, 5000)
}
