#' @title Correlation matrix of an autoregressive model of order 1
#' @name AR1
#' @description An order-1 autoregressive correlation matrix is set up which is
#'   used for the examples on power and sample size calculation.
#' @param p dimension
#' @param rho correlation
#' @return (p x p) matrix
#' @examples
#'   AR1(10, 0.2)
#' @export
AR1 <- function(p, rho){
  A <- diag(1, p)
  if(p > 1){
    for(i in 1:p){
      if(i < p) A[i, (i + 1):p] <- sapply(1:(p - i), function(k){rho ^ k})
      if(i > 1) A[i, 1:(i - 1)] <- sapply((i - 1):1, function(k){rho ^k})
    }
  }
  return(A)
}


#' @title Variance of estimator
#' @name calcvar
#' @description Calculation of variance of estimator and residual degrees of
#'   freedom
#' @details The variance of estimator beta (regression coefficient of SNP-BLUP
#'   approach) and the residual degrees of freedom are calculated based on the
#'   eigenvalue decomposition of correlation matrix \code{R}
#' @param lambda shrinkage parameter
#' @param eigendec eigenvalue decomposition of (p x p) correlation matrix
#'   \code{R}
#' @param n sample size
#' @param weights vector (LEN p) of SNP-specific weights or scalar if weights
#'   are equal for all SNPs; default value 1
#' @return
#' \describe{
#'   \item{\code{df}}{residual degrees of freedom}
#'   \item{\code{var.beta}}{vector (LEN p) of variance of estimator beta up to a
#'     constant (i.e. residual variance / n)}
#' }
#' @examples
#'   ### correlation matrix (should depend on sire haplotypes)
#'   R <- AR1(100, rho = 0.1)
#'   eigendec <- eigen(R)
#'   out <- calcvar(1200, eigendec, 100)
#' @export
calcvar <- function(lambda, eigendec, n, weights = 1){
  p <- length(eigendec$values)
  if (length(weights) < p) weights <- rep(weights[1], p)

  V <- eigendec$vectors
  df.res <- sum(eigendec$values * (eigendec$values + 2 * lambda / n) / (eigendec$values + lambda / n) ^2)
  vb <- NA
  if(df.res >= n) {
    message("WARNING: df.res > n ", round(df.res, 4)," ", n)
  } else{
    vb <- sapply(1:p, function(k) {1 / weights[k] * sum(V[k, ]^2 * eigendec$values / (eigendec$values + lambda / n) ^ 2)}) # up to a constant sigmaE2 / n
  }
  out <- list(df = df.res, var.beta = vb)
  return(out)
}


#' @title Ratio of expected value to variance of estimator
#' @name coeff.beta.k
#' @description The ratio of expected value to standard deviation is calculated
#'   for the estimator of a selected regression coefficient.
#' @param k index of selected regression coefficient
#' @param beta.true (LEN p) vector of regression coefficients
#' @param lambda shrinkage parameter
#' @param eigendec eigenvalue decomposition of (p x p) correlation matrix
#'   \code{R}
#' @param n sample size
#' @param weights vector (LEN p) of SNP-specific weights or scalar if weights
#'   are equal for all SNPs; default value 1
#' @return ratio
coeff.beta.k <- function(k, beta.true, lambda, eigendec, n, weights = 1){
  p <- length(eigendec$values)
  if (length(weights) < p) weights <- rep(weights[1], p)

  V <- eigendec$vectors
  d <- eigendec$values / (eigendec$values + lambda / n)
  z <- c()
  for(i in 1:length(d)) {
    z[i] <- sum(d * V[i, ] * V[k, ]) # coefficient for beta_i, then E(beta_k)=sum(coeff_(i,k)*beta_i)
  }
  return(1 / sqrt(weights[k]) * sum( z * sqrt(weights) * beta.true))
}

#' @title Start value for estimating optimal sample size
#' @name startvalue
#' @description Calculation of start value for estimating optimal sample size
#' @details Minimum sample size that exceeds residual degrees of freedom; this
#'   value can be used as start value in grid search for optimal sample size
#' @param lambda shrinkage parameter
#' @param R (p x p) matrix containing theoretical correlation between SNP pairs
#' @param nfam number of half-sib families
#' @param weights vector (LEN p) of SNP-specific weights or scalar if weights
#'   are equal for all SNPs; default value 1
#' @examples
#'   ### correlation matrix (should depend on sire haplotypes)
#'   R <- AR1(100, rho = 0.1)
#'   startvalue(1200, R, 10)
#' @return start value
#' @export
startvalue <- function(lambda, R, nfam, weights = 1){
  p <- nrow(R)
  if (length(weights) < p) weights <- rep(weights[1], p)
  if (sum(weights == 1) != p ){
    R <- t(apply(R, 1, function(x) {x / sqrt(weights)}))
    R <- apply(R, 2, function(x) {x / sqrt(weights)})
  } # leads to "near-symmetric" matrix
  e <- eigen(R)

  repeat{
    res <- calcvar(lambda, e, nfam, weights)
    if(res$df < nfam) break
    nfam <- nfam + 10
    if(nfam == 1000) {
      message('No start value obtained')
      break
    }
  }
  return(nfam)
}

#' @title Probability under alternative hypothesis (power)
#' @name pwr.normtest
#' @description Calculation of power is based on normal distribution.
#'  At each selected QTL position, the probability of the corresponding
#'  regression coefficient being different from zero is calculated using a
#'  t-like test statistic which has normal distribution with mean
#'  \code{E(beta_k)/sqrt{Var(beta_k)}}
#'  and variance 1. Under the null hypothesis \code{beta_k = 0},
#'  \code{E(beta_k) = 0}.
#'  Then, the mean value is returned as power.
#' @param R (p x p) matrix containing theoretical correlation between SNP pairs
#' @param n sample size
#' @param betaSE effect size relative to residual standard deviation
#' @param lambda shrinkage parameter
#' @param pos vector (LEN nqtl) of SNP indices for assumed QTL positions
#' @param weights weights vector (LEN p) of SNP-specific weights or scalar if weights
#'   are equal for all SNPs; default value 1
#' @param alpha type-I error level; default value 0.01
#' @return
#' \describe{
#'   \item{\code{result}}{mean power at selected QTL positions}
#'   \item{\code{h2.le}}{QTL heritability under linkage-equilibrium assumption}
#'   \item{\code{h2.ld}}{QTL heritability under linkage-disequilibrium
#'     assumption}
#' }
#' @references Wittenburg, Bonk, Doschoris, Reyer (2020) Design of Experiments
#'   for Fine-Mapping Quantitative Trait Loci in Livestock Populations.
#'   BMC Genetics 21:66. \doi{10.1186/s12863-020-00871-1}
#' @examples
#'   ### correlation matrix (should depend on sire haplotypes)
#'   R <- AR1(100, rho = 0.1)
#'   ### positions of putative QTL signals
#'   pos <- c(14, 75)
#'   ### power at given sample size and other parameters
#'   pwr.normtest(R, 100, 0.35, 1200, pos)
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @export
pwr.normtest <- function(R, n, betaSE, lambda, pos, weights = 1, alpha = 0.01){

  p <- nrow(R)
  if (length(weights) < p) weights <- rep(weights[1], p)
  if (sum(weights == 1) != p ){
    R <- t(apply(R, 1, function(x) {x / sqrt(weights)}))
    R <- apply(R, 2, function(x) {x / sqrt(weights)})
  } # leads to "near-symmetric" matrix
  e <- eigen(R)

  beta.true <- rep(0, p); beta.true[pos] <- betaSE
  x <- t(beta.true) %*% R %*% beta.true
  h2.ld <- x / (1 + x)
  h2.le <- sum(beta.true^2) / (1 + sum(beta.true^2))

  out <- calcvar(lambda, e, n, weights)
  threshold <- qnorm(1 - alpha / 2)
  exp.beta <- sapply(1:p, function(k) {coeff.beta.k(k, beta.true, lambda, e, n, weights)})
  term <- exp.beta / sqrt(out$var.beta) * sqrt(n)
  pwr <- pnorm(threshold - term, lower.tail = F) + pnorm(-threshold - term)
  result <- mean(pwr[pos])

  return(c(result, h2.le, h2.ld))
}



#' @title Method of bisection for estimating optimal sample size
#' @name search.best.n.bisection
#' @description A grid \code{[nstart, nmax]} for possible sample size is
#'   considered. Instead of executing a time-consuming grid search, the method
#'   of bisection is applied to this interval. For each step, the function
#'   \code{pwr.normtest} is called for the given set of parameters.
#' @param R (p x p) matrix containing theoretical correlation between SNP pairs
#' @param betaSE effect size relative to residual standard deviation
#' @param lambda shrinkage parameter
#' @param pos vector (LEN nqtl) of SNP indices for assumed QTL positions
#' @param nstart minimum value for grid search
#' @param nmax maximum value for grid search
#' @param weights vector (LEN p) of SNP-specific weights or scalar if weights
#'   are equal for all SNPs; default value 1
#' @param typeII type-II error level; default value 0.2
#' @param alpha type-I error level; default value 0.01
#' @return integer of optimal sample size
#' @examples
#'   ### correlation matrix (should depend on sire haplotypes)
#'   R <- AR1(100, rho = 0.1)
#'   ### positions of putative QTL signals
#'   pos <- c(14, 75)
#'   ### optimal sample size
#'   search.best.n.bisection(R, 0.35, 1200, pos, 10, 5000)
#' @export
search.best.n.bisection <- function(R, betaSE, lambda, pos, nstart, nmax, weights = 1, typeII = 0.2, alpha = 0.01){
  y <- c(); k <- 0
  n1 <- nstart
  n3 <- nmax
  n2 <- ceiling(mean(c(n1, n3)))

  y[n1] <- pwr.normtest(R, n1, betaSE, lambda, pos, weights, alpha)[1]
  y[n2] <- pwr.normtest(R, n2, betaSE, lambda, pos, weights, alpha)[1]
  y[n3] <- pwr.normtest(R, n3, betaSE, lambda, pos, weights, alpha)[1]

  repeat{
    if(y[n3] < 1 - typeII) {
      message('WARNING: power ', 1- typeII, ' can not be achieved. Increase nmax.')
      return(n3)
      break
    }

    if(all(y[c(n1, n2, n3)] >= 1 - typeII)){
      return(n1)
      break

    } else if(all(y[c(n2, n3)] >= 1 - typeII)){
      n3 <- n2
      n2 <- ceiling(mean(c(n1, n3)))
      y[n2] <- pwr.normtest(R, n2, betaSE, lambda, pos, weights, alpha)[1]

    } else if(y[n3] >= 1 - typeII){
      n1 <- n2
      n2 <- ceiling(mean(c(n1, n3)))
      y[n2] <- pwr.normtest(R, n2, betaSE, lambda, pos, weights, alpha)[1]
    }

    k <- k + 1
    if(k > log2(nmax - nstart)) {
      return(n2)
      break
    }
  }
}

#' @title Wrapper function for sample size calculation
#' @name pwr.snpblup
#' @description Given parameters specified by the experimenter, optimal sample
#'   size is estimated by repeatedly applying \code{search.best.n.bisection}.
#' @details  Sample size depends on parameters specified by the experimenter
#'   (number of half-sib families, number of QTL, heritability, correlation
#'   matrix). These values are converted into parameters required for the
#'   probability density function under the alternative hypothesis (beta_k !=0,
#'   for k selected QTL positions). As power depends on the selected QTL
#'   positions, these are sampled at random and power calculations are repeated.
#'   Afterwards the mean value is a plausible estimate of optimal sample size.
#'
#'   Linear model for SNP-BLUP approach:
#'   \code{y = X beta + e}
#'   with \code{t(beta) = (beta_1, ldots, beta_p)}
#'   Ridge approach:
#'   \code{hat{beta} = (Xt X + I lambda)^{-1} Xt y}
#'
#'   The identity matrix \code{I} can be replaced by a diagonal matrix
#'   containing SNP-specific weights yielding a generalised ridge approach.
#'
#' @param nfathers number of half-sib families
#' @param nqtl number of QTL assumed
#' @param h2 heritability captured by QTL
#' @param R (p x p) matrix containing theoretical correlation between SNP pairs
#' @param weights vector (LEN p) of SNP-specific weights or scalar if weights
#'   are equal for all SNPs; default value 1
#' @param rep number of repetitions; default value 10
#' @param nmax maximum value for grid search; default value 5000
#' @param alpha type-I error level; default value 0.01
#' @param typeII type-II error level; default value 0.2
#' @return vector of optimal sample size over all repetitions
#' @references Wittenburg, Bonk, Doschoris, Reyer (2020) Design of Experiments
#'   for Fine-Mapping Quantitative Trait Loci in Livestock Populations. BMC
#'   Genetics 21:66. \doi{10.1186/s12863-020-00871-1}
#' @examples
#'   ### input parameters specified by experimenter
#'   # number of half-sib families
#'   nfathers <- 10
#'   # number of assumed QTL
#'   nqtl <- 2
#'   # QTL heritability
#'   h2 <- 0.2
#'   ### correlation matrix (should depend on sire haplotypes)
#'   R <- AR1(100, rho = 0.1)
#'   ### optimal sample size in a multi-marker approach
#'   set.seed(11)
#'   pwr.snpblup(nfathers, nqtl, h2, R, rep = 1)
#' @export
pwr.snpblup <- function(nfathers, nqtl, h2, R, rep = 10, nmax = 5000, weights = 1, typeII = 0.2, alpha = 0.01){
  q <- nrow(R)

  ### derived input parameters for statistics
  lambda <- round((1 - h2) / h2 * q)
  betaSE <- sqrt(h2 / (nqtl - nqtl * h2))
  ### start value for sample size (to ensure that degrees of freedom are positive)
  nstart <- startvalue(lambda, R, nfathers, weights)

  nopt <- c()
  for(i in 1:rep){
    ### random QTL positions
    pos <- sample(1:nrow(R), nqtl)
    nopt[i] <- search.best.n.bisection(R, betaSE, lambda, pos, nstart, nmax, weights, typeII, alpha)
  }
  return(nopt)
}


#' @title Calculation of effective number of independent tests
#' @name simpleM
#' @description Adapted simpleM method which considers theoretical correlation
#'   between SNP pairs instead of composite LD values. Principal component
#'   decomposition yields the effective number of independent tests. This value
#'   is needed for the Bonferroni correction of type-I error when testing SNP
#'   effects based on a single-marker model.
#' @references
#'   Gao, Starmer & Martin (2008) A multiple testing correction method for
#'   genetic association studies using correlated single nucleotide
#'   polymorphisms. Genetic Epidemiology 32:361-369.
#' @param mat correlation matrix
#' @param quant percentage cutoff, variation of SNP data explained by
#'   eigenvalues; default value 0.995
#' @return effective number of independent tests
#' @importFrom pwr pwr.t.test
#' @examples
#'   ### correlation matrix (should depend on sire haplotypes)
#'   R <- AR1(100, rho = 0.1)
#'   ### effective number of tests
#'   Meff <- simpleM(R)
#'   ### relative effect size given heritability and number of QTL signals
#'   h2 <- 0.2
#'   nqtl <- 2
#'   betaSE <- sqrt(h2 / (nqtl - nqtl * h2))
#'   ### optimal sample size in a single-marker approach
#'   pwr::pwr.t.test(d = betaSE, sig.level = 0.01 / Meff, power = 0.8,
#'    alternative = "two.sided", type = "one.sample")
#' @export
simpleM <- function(mat, quant = 0.995){
  e <- eigen(mat, only.values = T)$values
  s <- sum(e)
  p <- length(e)
  ss <- sapply(1:p, function(x) {sum(e[1:x])/ s})
  m <- min(which(ss >= quant))
  return(m)
}
