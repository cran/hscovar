# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                               Haplo2Geno
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#' @title Conversion of haplotypes into genotypes
#' @name Haplo2Geno
#' @description Haplotypes are converted into into genotypes without checking
#'  for missing values.
#' @param inpMat [MATRIX] haplotype matrix (2 lines per individual)
#' @return
#' \describe{
#'   \item{\code{outMa}}{(N x p) genotype matrix}
#' }
#' @examples
#'  data(testdata)
#'  G <- Haplo2Geno(H.sire)
#' @export
Haplo2Geno = function( inpMat ) {
  n = dim( inpMat )[1]
  p = dim( inpMat )[2]
  outMat = matrix( NA, nrow = n/2, ncol = p)
  for( i in 1:(n/2) ) outMat[ i, ] = inpMat[ (2*i - 1), ] + inpMat[ 2*i, ]
  colnames(outMat) = colnames(inpMat)     # Retain original column names
  return( outMat )
}


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                               Expectation matrices
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#' @title Expected value of paternally inherited allele
#' @name ExpectMat
#' @description Expected value is +/-0.5 if sire is homozygous reference/
#'   alternate allele or 0 if sire is heterozygous at the investigated marker
#' @param inMat [MATRIX] The paternal genotype matrix
#' @return
#' \describe{
#'   \item{\code{ExP.Fa}}{(N x p) matrix of expected values}
#' }
#' @examples
#'  data(testdata)
#'  G <- Haplo2Geno(H.sire)
#'  E <- ExpectMat(G)
#' @export
ExpectMat = function( inMat) {
  rule = function( index ) ifelse( index == 2 , 0.5, ifelse( index == 0, -0.5, ifelse(index == 1, 0, NA)) )
  ExP.Fa = apply( inMat, 2, rule )
  return( ExP.Fa )
}


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                               LD matrices
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#' @title Calculation of maternal LD matrix
#' @name LDdam
#' @description Matrix containing linkage disequilibrium between marker pairs
#'   on maternal gametes is set up by counting haplotypes frequencies.
#' @details The function generates a block diagonal sparse matrix based on
#'   Matrix::bdiag. Use as.matrix() to obtain a regular one.
#' @param inMat [MATRIX] The maternal HAPLOTYPE matrix.
#' @param pos_chr [LIST] The marker positions in Morgan on chromosomes.
#' @return
#' \describe{
#'   \item{\code{Dd}}{(p x p) matrix of maternal LD}
#' }
#' @examples
#'  ## haplotype matrix of n individuals at p SNPs
#'  p <- 10; n <- 4
#'  mat <- matrix(ncol = p, nrow = 2 * n, sample(c(0, 1), size = 2 * n * p, replace = TRUE))
#'  LDdam(mat, list(1:p))
#' @import Matrix
#' @export
LDdam = function( inMat, pos_chr ){

  length_chr = vector( mode = "list", length = length( pos_chr ) )
  inMatByChr = vector( mode = "list", length = length( pos_chr ) )
  Dd = vector( mode = "list", length = length( pos_chr ) )

  for( nc in seq_along( pos_chr ) ) {

    length_chr[[ nc ]] = length( pos_chr[[ nc ]] )
    # Split inMat by chromosome
    if( nc == 1 ) inMatByChr[[ nc ]] = inMat[ , 1:sum( unlist( length_chr[ 1:nc ] ) ) ]
    else inMatByChr[[ nc ]] = inMat[ , ( sum( unlist( length_chr[ 1:( nc - 1 ) ] ) ) + 1 ):sum( unlist( length_chr[ 1:nc ] ) ) ]

    n = dim( inMatByChr[[ nc ]] )[1]
    p = dim( inMatByChr[[ nc ]] )[2]
    haplo.freq = array( numeric(), dim = c( 2, 2, p, p ) )
    Dd[[ nc ]] = matrix( NA, nrow = p, ncol = p )
    for( j in 1:p ) {
      for( k in j:p ) { # Calculates the upper triangular part plus diagonal
        zz = sum( !inMatByChr[[ nc ]][ , j ] & !inMatByChr[[ nc ]][ , k ] ) # zero-zero
        oz = sum( inMatByChr[[ nc ]][ , j ] & !inMatByChr[[ nc ]][ , k ] )  # one-zero
        zo = sum( !inMatByChr[[ nc ]][ , j ] &  inMatByChr[[ nc ]][ , k ] ) # zero-one
        oo = sum( inMatByChr[[ nc ]][ , j ] &  inMatByChr[[ nc ]][ , k ] )  # one-one
        haplo.freq[ , , j, k ] = 1/n * matrix( c( zz, oz, zo, oo ), ncol = 2 )
        Dd[[ nc ]][ j, k ] = det( haplo.freq[ , , j, k ] )
      }
    }
    Dd[[ nc ]] = Matrix::forceSymmetric( Dd[[ nc ]] ) # The whole matrix
  }
  return( Reduce( Matrix::bdiag, Dd ) ) # Creates a block diagonal sparse matrix
}


#' @title Calculation of paternal LD matrix
#' @name LDsire
#' @description Matrix containing linkage disequilibrium between marker pairs
#'   on paternal gametes is set up from sire haplotypes and genetic-map
#'   information for each half-sib family.
#' @details The function generates a block diagonal sparse matrix based on
#'   Matrix::bdiag. Use as.matrix() to obtain a regular one.
#' @param inMat [MATRIX] Haplotype matrix for sires for all chromosomes.
#' @param pos_chr [LIST] The marker positions in Morgan on chromosomes.
#' @param family [VECTOR] Which family (sire) should be processed?
#'   Vector with consecutive entries of the form 1:2, 3:4, 5:6 and so on,
#'   linking to haplotypes (rows in inMat) of the corresponding sire
#' @param map_fun ["haldane" or "kosambi"] The mapping function applied.
#' @return Ds
#' \describe{
#'   \item{\code{Ds}}{(p x p) matrix of paternal LD}
#' }
#' @examples
#'  data(testdata)
#'  LDfam2 <- LDsire(H.sire, pos.chr, family = 3:4)
#' @import Matrix
#' @export
LDsire = function( inMat, pos_chr, family, map_fun = "haldane" ) {

  # Mapping functions
  # 1. Haldane
  if( map_fun == "haldane" )
    theta = function( pos_chr, nc, j, k ) {
      0.5 * ( 1 - exp( -2 * ( abs( pos_chr[[ nc ]][k] - pos_chr[[ nc ]][j] ) ) ) )
    }
  # 2. Kosambi
  if( map_fun == "kosambi" )
    theta = function( pos_chr, nc, j, k ) {
      0.5 * tanh( 2 * ( abs( pos_chr[[ nc ]][k] - pos_chr[[ nc ]][j] ) ) )
    }

  if( length( family ) != 2 ) stop( "Number of rows not correct." )
  inMat = inMat[ family, ]

  length_chr = vector( mode = "list", length = length( pos_chr ) )
  inMatByChr = vector( mode = "list", length = length( pos_chr ) )
  Ds = vector( mode = "list", length = length( pos_chr ) )
  for( nc in seq_along( pos_chr ) ) {
    length_chr[[ nc ]] = length( pos_chr[[ nc ]] )
    # Split inMat by chromosome
    if( nc == 1 ) inMatByChr[[ nc ]] = inMat[ , 1:length( pos_chr[[ nc ]] ) ]
    else inMatByChr[[ nc ]] = inMat[ , ( sum( unlist( length_chr[ 1:( nc - 1 ) ] ) ) + 1 ):sum( unlist( length_chr[ 1:nc ] ) ) ]

    p = dim( inMatByChr[[ nc ]] )[2]
    Ds[[ nc ]] = matrix( NA, nrow = p, ncol = p )
    for( j in 1:p ) {
      for( k in j:p ) { # Calculates the upper triangular part plus diagonal
        if(!anyNA(inMatByChr[[ nc ]][ , c(j, k) ] )){
          if( sum( inMatByChr[[ nc ]][ , j ] ) == 1 & sum( inMatByChr[[ nc ]][ , k ] ) == 1 ) {
            if( sum( inMatByChr[[ nc ]][ 1, j ], inMatByChr[[ nc ]][ 1, k ] ) == 1 ) Ds[[ nc ]][ j, k ] = -0.25 * ( 1 - 2 * theta( pos_chr, nc, j, k ) )
            else Ds[[ nc ]][ j, k ] = 0.25 * ( 1 - 2 * theta( pos_chr, nc, j, k ) )
          } else Ds[[ nc ]][ j, k ] = 0
        } else Ds[[ nc ]][ j, k ] = NA
      }
    }
    Ds[[ nc ]] = Matrix::forceSymmetric( Ds[[ nc ]] ) # The whole matrix
  }
  return( Reduce( Matrix::bdiag, Ds ) ) # Creates a block diagonal sparse matrix
}


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                               Covariance matrices
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#' @title Calculation of covariance matrices from maternal and paternal LD
#' @name CovarMatrix
#' @description The covariance matrix is set as maternal plus paternal LD
#'   matrix where the paternal part is a weighted average of sire-specific LD
#'   matrices.
#' @details The internal suMM function works on lists!
#' @param exp_freq_mat [MATRIX] paternal EXPECTATION matrix
#' @param LDDam [MATRIX] maternal Linkage Disequilibrium matrix
#' @param LDSire [LIST] Linkage disequilibrium matrices for the sires; each
#'   element of the list corresponds to a family
#' @param Ns [VECTOR] family size for each sire s
#' @return
#' \describe{
#'   \item{\code{covK}}{(p x p) matrix of covariance between markers}
#' }
#' @examples
#'  data(testdata)
#'  G <- Haplo2Geno(H.sire)
#'  E <- ExpectMat(G)
#'  LDfam2 <- LDsire(H.sire, pos.chr, family = 3:4)
#'  LDfam3 <- LDsire(H.sire, pos.chr, family = 5:6)
#'  ## covariance matrix based on sires 2 and 3 only, each with 100 progeny
#'  K <- CovarMatrix(E[2:3, ], LDDam = matLD, LDSire = list(LDfam2, LDfam3), Ns = c(100, 100))
#' @import parallel
#' @export
CovarMatrix = function( exp_freq_mat, LDDam, LDSire, Ns ) {

  suMM = function( InsertList ) {
    InsertList = lapply(InsertList, function(z){z <- replace(z, is.na(z), 0)})
    Reduce( "+", Map( "*", InsertList, Ns/sum( Ns ) ) )
  }

  #        --- --- --- Process the EXPECTATION matrices --- --- ---
  e1 = list()
  for( i in 1:dim( exp_freq_mat )[1] ) { e1[[i]] = exp_freq_mat[ i, ] }
  e2 = mclapply( e1, tcrossprod )

  #                  --- --- --- Compute the COVARIANCE matrix --- --- ---
  sm = suMM( e1 )
  covK = LDDam + suMM( LDSire ) + suMM( e2 ) - tcrossprod( sm )
  return( covK )
}




