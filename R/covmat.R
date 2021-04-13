#' @title Calculation of covariance or correlation matrix
#' @name CovMat
#' @description The theoretical covariance between pairs of markers is
#'   calculated from either paternal haplotypes and maternal linkage
#'   disequilibrium (LD) or vise versa. A genetic map is required. The
#'   implementation relies on paternal half-sib families and biallelic markers
#'   such as single nucleotide polymorphisms (SNP). If parental haplotypes are
#'   incomplete (i.e., SNP alleles are missing), those parents will be discarded
#'   at the corresponding pairs of SNPs. If maternal half-sib families are used,
#'   the roles of sire/dam are swapped. Multiple families can be considered.
#' @note Family size is used for weighting covariance terms in case of multiple
#'   half-sib families. It only matters if number of progeny differs.
#'
#'   If maternal haplotypes (H.mothers) are used instead of maternal LD (matLD),
#'   LD can be estimated from counting haplotype frequencies as:
#'
#'   \code{matLD <- LDdam(inMat = H.mother, pos.chr)}
#'
#'   If multiple chromosomes are considered, SNP positions are provided as, e.g.
#'
#'   \code{pos.chr <- list(pos.snp.chr1, pos.snp.chr2, pos.snp.chr3)}
#' @param linkMat (p x p) matrix of maternal LD between pairs of p markers;
#'   matrix is block diagonal in case of multiple chromosomes and must not
#'   contain missing values; use zeros if LD is uncertain
#' @param haploMat (2N x p) matrix of sires haplotypes for all chromosomes
#'   (2 lines per sire); coding with 0's and 1's reflecting reference and
#'   alternate alleles, respectively; missing values can be coded as NA or any
#'   integer but not 0 and 1
#' @param nfam vector (LEN N) containing number of progeny per family or
#'   scalar value in case of equal family size
#' @param pos_chr list (LEN number of chromosomes) of vectors (LEN number of
#'   markers) of genetic positions in Morgan per chromosome
#' @param map_fun character string of mapping function used; so far "haldane"
#'   (default) and "kosambi" are enabled
#' @param corr logical; \code{FALSE} (default) if output is covariance matrix or
#'   \code{TRUE} if output is correlation matrix
#' @return list (LEN 2) of matrix (DIM \eqn{p1} x \eqn{p1}) and vector
#'   (LEN \eqn{p1}) with \eqn{p1 \le p}
#' \describe{
#'   \item{\code{K}}{covariance matrix OR}
#'   \item{\code{R}}{correlation matrix}
#'   \item{\code{valid.snps}}{vector of SNP indices considered for covariance/
#'   correlation matrix}
#' }
#' @examples
#'   ### 1: INPUT DATA
#'   data(testdata)
#'   ### 2: COVARIANCE/CORRELATION MATRIX
#'   corrmat <- CovMat(matLD, H.sire, 100, pos.chr, corr = TRUE)
#'   ### 3: TAGSNPS FROM CORRELATION MATRIX
#'   bin <- tagSNP(corrmat$R)
#'   bin <- tagSNP(corrmat$R, 0.5)
#' @references Wittenburg, Bonk, Doschoris, Reyer (2020) Design of Experiments
#'   for Fine-Mapping Quantitative Trait Loci in Livestock Populations. BMC
#'   Genetics 21:66. \doi{10.1186/s12863-020-00871-1}
#' @import foreach Matrix
#' @export
CovMat <- function(linkMat, haploMat, nfam, pos_chr, map_fun = 'haldane', corr = F){
  N <- nrow(haploMat) / 2
  p <- length(unlist(pos_chr))

  haploMat <- apply(haploMat, 2, function(z){z[!(z %in% c(0, 1))] <- NA; z})

  # genotypes of sires
  XSire <- matrix(Haplo2Geno(as.matrix(haploMat)), nrow = N, ncol = p)

  # expection of paternally inherited SNP allele
  expectationMat <- matrix(ExpectMat(XSire), nrow = N, ncol = p)

  # family index
  fam <- vector(mode = "list", length = length(1:N))
  for(l in 1:N) fam[[l]] <- (2 * l - 1 ):(2 * l)

  # LD of paternally inherited SNP alleles; output is list -> matrices for each family
  indexFam <- NULL
  linkSire <- foreach(indexFam = 1:N) %do% LDsire(inMat = haploMat, pos_chr, family = fam[[indexFam]], map_fun = map_fun)

  # Weighted average over paternal half-sib families plus maternal LD
  if(length(nfam) == N) Ns <- nfam else if (length(nfam) == 1) Ns <- rep(nfam, N) else stop("ERROR family size")
  K <- CovarMatrix(expectationMat, linkMat, linkSire, Ns)

  s <- Matrix::diag(K)
  id <- s > 1e-6

  # feedback
  if (sum(id) < p) message(paste(sum(id), 'SNPs have non-zero variance -> output matrix has reduced dimensions.'))

  colnames(K) <- rownames(K) <- 1:p
  if(corr){
    K <- K[id, id]; s <- s[id]
    R <- t(apply(K, 1, function(x) {x / sqrt(s)}))
    R <- apply(R, 2, function(x) {x / sqrt(s)})
    return(list(R = R, valid.snps = which(id)))
  } else return(list(K = K[id, id], valid.snps = which(id)))
}


