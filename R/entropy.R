#' Extract model entropy
#'
#' Extract the entropy of a Maxent model from a \code{lambdas} objects generated
#' by \code{\link{parse_lambdas}}.
#'
#' @param x A \code{lambdas} object returned by \code{\link{parse_lambdas}}.
#' @keywords maxent, lambdas
#' @seealso \code{\link{parse_lambdas}} \code{\link{project}}
#' @export
#' @examples
#' # Below we use the dismo::maxent example to fit a Maxent model:
#' if (require(dismo) && require(rJava) && 
#'     file.exists(system.file('java/maxent.jar', package='dismo'))) {
#'   fnames <- list.files(system.file('ex', package="dismo"), '\\.grd$', 
#'                        full.names=TRUE )
#'   predictors <- stack(fnames)
#'   occurrence <- system.file('ex/bradypus.csv', package='dismo')
#'   occ <- read.table(occurrence, header=TRUE, sep=',')[,-1]
#'   me <- maxent(predictors, occ, path=file.path(tempdir(), 'example'), 
#'                factors='biome')
#' 
#'   # ... and then parse the lambdas information:
#'   lam <- parse_lambdas(me)
#'   entropy(lam)
#' }
#' 
entropy <- function(x) x$entropy
