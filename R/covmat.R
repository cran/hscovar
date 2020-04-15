#' @title Calculation of covariance or correlation matrix
#' @name CovMat
#' @description The theoretical covariance between pairs of markers is calculated
#'   from either paternal haplotypes and maternal linkage disequilibrium (LD) or
#'   vise versa. A genetic map is required. The implementation relies on
#'   paternal half-sib families and biallelic markers such as single nucleotide
#'   polymorphisms (SNP).
#' @note If maternal half-sib families are used, the roles of sire/dam are swapped.
#'   Multiple families can be considered.
#'
#'   Family size is used for weighting covariance terms in case of multiple
#'   half-sib families. It only matters if number of progeny differs.
#' @param linkDam (p x p) matrix of maternal LD between pairs of p markers;
#'   matrix is block diagonal in case of multiple chromosomes
#' @param haploSire (2N x p) matrix of sires haplotypes for all chromosomes
#'   (2 lines per sire)
#' @param nfam vector (LEN N) containing number of progeny per sire or
#'   scalar value in case of equal family size
#' @param pos_chr list (LEN number of chromosomes) of vectors (LEN number of
#'   markers) of genetic positions in Morgan per chromosome
#' @param corr logical; \code{TRUE} (default) if output is correlation matrix or
#'   \code{FALSE} if output is covariance matrix
#' @return list (LEN 2) of matrix (DIM \eqn{p1} x \eqn{p1}) and vector
#'   (LEN \eqn{p1}) with \eqn{p1 \le p}
#' \describe{
#'   \item{\code{K}}{covariance matrix OR}
#'   \item{\code{R}}{correlation matrix}
#'   \item{\code{valid.snps}}{vector of SNP indices considered for covariance/
#'   correlation matrix}
#' }
#' @note{
#'   If you have maternal haplotypes (H.mothers; same format as H.sire)
#'   instead of maternal LD (matLD) then LD can be estimated from counting
#'   haplotype frequencies as:
#'
#'   \code{matLD <- LDdam(inMat = H.mother, pos.chr)}
#'
#'   If multiple chromosomes are considered, then, for instance:
#'
#'   \code{pos.chr <- list(pos.snp.chr1, pos.snp.chr2, pos.snp.chr3)}
#' }
#' @examples
#'   ### 1: INPUT DATA
#'   data(testdata)
#'   ### 2: COVARIANCE/CORRELATION MATRIX
#'   corrmat <- CovMat(matLD, H.sire, 100, pos.chr, corr = TRUE)
#'   ### 3: TAGSNPS FROM CORRELATION MATRIX
#'   bin <- tagSNP(corrmat$R)
#'   bin <- tagSNP(corrmat$R, 0.5)
#' @references Wittenburg, Bonk, Doschoris, Reyer (2019) "Design of Experiments
#'   for Fine-Mapping Quantitative Trait Loci in Livestock Populations"
#'   \url{https://doi.org/10.1101/2019.12.17.879106}
#' @import foreach
#' @export
CovMat <- function(linkDam, haploSire, nfam, pos_chr, corr = T){
  N <- nrow(haploSire) / 2
  p <- length(unlist(pos_chr))

  # genotypes of sires
  XSire <- matrix(Haplo2Geno(as.matrix(haploSire)), nrow = N, ncol = p)

  # expection of paternally inherited SNP allele
  expectationMat <- matrix(ExpectMat(XSire), nrow = N, ncol = p)

  # family index
  fam <- vector(mode = "list", length = length(1:N))
  for(l in 1:N) fam[[l]] <- (2 * l - 1 ):(2 * l)

  # LD of paternally inherited SNP alleles; output is list -> matrices for each family
  indexFam <- NULL
  linkSire <- foreach(indexFam = 1:N) %do% LDsire(inMat = haploSire, pos_chr, family = fam[[indexFam]], map_fun = "haldane")

  # Weighted average over paternal half-sib families plus maternal LD
  if(length(nfam) == N) Ns <- nfam else if (length(nfam) == 1) Ns <- rep(nfam, N) else stop("ERROR family size")
  K <- CovarMatrix(expectationMat, linkDam, linkSire, Ns)

  s <- diag(K)
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


