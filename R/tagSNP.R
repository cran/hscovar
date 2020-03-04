#' @title tagSNP
#' @description Grouping of markers depending on correlation structure
#' @details Grouping of markers is based on the  correlation matrix. Apart from
#'   this, the strategy for grouping is similar to Carlson et al. (2004). A
#'   representative marker is suggested for each group.
#' @param mat (p x p) correlation matrix
#' @param threshold lower value of correlation considered for grouping
#' @return list (LEN number of groups) of lists (LEN 2); marker names correspond
#'   to column names of mat
#' \describe{
#'   \item{\code{snps}}{vector of marker IDs in group}
#'   \item{\code{tagsnp}}{representative marker suggested for this group}
#' }
#' @examples
#'   ### 1: INPUT DATA
#'   data(testdata)
#'   ### 2: COVARIANCE/CORRELATION MATRIX
#'   corrmat <- CovMat(matLD, H.sire, 100, pos.chr, corr = TRUE)
#'   ### 3: TAGSNPS FROM CORRELATION MATRIX
#'   bin <- tagSNP(corrmat$R)
#'   bin <- tagSNP(corrmat$R, 0.5)
#'   as.numeric(unlist(rlist::list.select(bin, tagsnp)))
#' @references Carlson, C. S., Eberle, M. A., Rieder, M. J., Yi, Q., Kruglyak,
#'   L. & Nickerson, D. A. Selecting a maximally informative set of single-
#'   nucleotide polymorphisms for association analyses using linkage
#'   disequilibrium. Am. J. Hum. Genet., 2004, 74:106-120.
#' @import rlist
#' @export
tagSNP <- function(mat, threshold = 0.8){
  if(max(abs(mat), na.rm = T) > 1 + 1e-6) stop("Correlation or R-squared matrix is expected.")
  # substitute conspicious values by zero
  mat <- apply(mat, 1, function(x) {x[!is.finite(x)] <- 0; return(x)}); diag(mat) <- 1
  if ((threshold >= 1) | (threshold <= 0)) stop("Threshold out of range")
  # initialise
  p <- nrow(mat)
  bin <- list()
  counter <- 0
  snpset <- 1:p
  coln <- colnames(mat)
  repeat{
    counter <- counter + 1
    # preliminary bins
    ls <- lapply(snpset, function(i){which(abs(mat[i, ]) > threshold)})
    # select largest bin only
    m <- which.max(lapply(ls, length))
    # tagSNP within largest bin
    if (length(ls[[m]]) == 1) {
      ts <- ls[[m]]
    } else {
      candidate <- apply(mat[ls[[m]], ls[[m]]], 1, function(x){sum(abs(x) > threshold) == length(ls[[m]])})
      ts <- ls[[m]][candidate][ceiling(sum(candidate) / 2)]
    }
    # output list uses original names of markers
    snps <- NULL
    bin[[counter]] <- list(snps = coln[ls[[m]]], tagsnp = coln[ts])
    # remaining SNPs for next iteration
    snpset <- setdiff(snpset, ls[[m]])
    if ((!length(snpset) > 0) || (counter == p)) break
  }
  # feedback
  z <- length(unique(unlist(list.select(bin, tp = snps))))
  message(paste(z, "SNPs have been grouped into", counter, "bins"))
  return(bin)
}
