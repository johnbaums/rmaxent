#' Extract model entropy
#'
#' Extract the entropy of a Maxent model from a `lambdas` objects generated
#' by [parse_lambdas()].
#'
#' @param x A `lambdas` object returned by [parse_lambdas()].
#' @keywords maxent, lambdas
#' @seealso [parse_lambdas()] [project()]
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
