#### Package Documentation and Imports ####

#' The sumrep package
#' 
#' Description of the package.
#' 
#' @name     sumrep
#' 
#' @import ggplot2
#' @import methods
#' @import stats
#' @import utils
#' @importFrom  alakazam     checkColumns gravy aminoAcidProperties calcDiversity
#' @importFrom  ape          as.DNAbin GC.content
#' @importFrom  Biostrings   DNAStringSet vcountPattern
#' @importFrom  CollessLike  sackin.index colless.like.index cophen.index
#' @importFrom  crayon       green yellow
#' @importFrom  data.table   data.table fread fwrite as.data.table
#' @importFrom  dplyr        last
#' @importFrom  entropy      KL.empirical
#' @importFrom  grid         grid.newpage grid.layout pushViewport viewport
#' @importFrom  HDMD         FactorTransform
#' @importFrom  jsonlite     fromJSON
#' @importFrom  magrittr     %>% %$% extract
#' @importFrom  Peptides     kideraFactors aIndex
#' @importFrom  seqinr       read.fasta write.fasta translate
#' @importFrom  shazam       distToNearest findThreshold calcBaseline summarizeBaseline
#'                           createSubstitutionMatrix createMutabilityMatrix
#' @importFrom  stringdist   stringdistmatrix
#' @importFrom  stringr      str_sub str_split
#' @importFrom  textmineR    CalcJSDivergence
#' @importFrom  yaml         yaml.load yaml.load_file
NULL


#### Data ####

#' Test data
#'
#' Data set desription.
#'
#' @format Data format details.
#' \itemize{
#'   \item{value1}  {description 1}
#'   \item{value2}  {description 2}
#' }
#' 
#' @name test_dat
NULL

#' Test data
#'
#' Data set desription.
#'
#' @format Data format details.
#' \itemize{
#'   \item{value1}  {description 1}
#'   \item{value2}  {description 2}
#' }
#' 
#' @name test_dat_boot
NULL

#' Test simulations
#'
#' Data set desription.
#'
#' @format Data format details.
#' \itemize{
#'   \item{value1}  {description 1}
#'   \item{value2}  {description 2}
#' }
#'
#'@name test_simu
NULL