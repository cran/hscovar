#' @title Description of the testdata
#' @name testdata
#' @docType data
#' @description The data set contains paternal haplotypes, maternal LD and
#'   genetic map positions that are required to calculate the covariance between
#'   pairs of markers.
#' \describe{
#'   The raw data can be downloaded at the source given below. Then,
#'   executing the following R code leads to the data that have been provided as
#'   \code{testdata.RData}.
#' \item{H.sire}{(2N x p) haplotype matrix for sires for all chromosomes
#'   (2 lines per sire)}
#' \item{matLD}{(p x p) matrix of maternal LD between pairs of p markers;
#'   matrix is block diagonal in case of multiple chromosomes}
#' \item{pos.chr}{list of vectors of genetic map positions per chromosome}
#' }
#' @source The data are available from the RADAR repository
#'   \url{https://dx.doi.org/10.22000/280}
#' @examples
#' \donttest{
#' ## data.frame of estimates of paternal recombination rate and maternal LD
#' load('Result.RData')
#' ## list of haplotypes of sires for each chromosome
#' load('sire_haplotypes.RData')
#' ## physical map
#' map <- read.table('map50K_ARS_reordered.txt', header = T)
#' ## select target region
#' chr <- 1
#' window <- 301:600
#' ## map information of target region
#' map.target <- map[map$Chr == chr, ][window, ]
#' Result.target <- Result[(Result$Chr == chr) & (Result$SNP1 %in% window) &
#'   (Result$SNP2 %in% window), ]
#' ## SNP position in Morgan approximated from recombination rate
#' part <- Result.target[Result.target$SNP1 == window[1], ]
#' sp <- smooth.spline(x = map.target$locus_Mb[part$SNP2 - window[1] + 1], y = part$Theta, df = 4)
#' pos.snp <- predict(sp, x =  map.target$locus_Mb[window - window[1] + 1])$y
#' ## list of SNPs positions
#' pos.chr <- list(pos.snp)
#' ## haplotypes of sires (mating candidates) in target region
#' H.sire <- rlist::list.rbind(haps[[chr]])[, window]
#' ## matrix of maternal LD (block diagonal if multiple chromosome)
#' matLD <- matrix(0, ncol = length(window), nrow = length(window))
#' ## off-diagonal elements
#' for(l in 1:nrow(Result.target)){
#'   id1 <- Result.target$SNP1[l] - window[1] + 1
#'   id2 <- Result.target$SNP2[l] - window[1] + 1
#'   matLD[id1, id2] <- matLD[id2, id1] <- Result.target$D[l]
#' }
#' ## diagonal elements
#' for(k in unique(Result.target$SNP1)){
#'   id <- k - window[1] + 1
#'   p <- Result.target$fAA[Result.target$SNP1 == k] + Result.target$fAB[Result.target$SNP1 == k]
#'   matLD[id, id] <- max(p * (1 - p))
#' }
#' }
#' @importFrom rlist list.rbind
NULL
