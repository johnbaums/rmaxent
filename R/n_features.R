#' Extract the number of features used in a Maxent model
#'
#' Calculate the number of non-zero features used in a Maxent model.
#'
#' @param lambdas A `vector` of one or more file paths to Maxent .lambdas 
#'   files, or a `list` of one or more `MaxEnt` fitted model objects.
#' @return An vector with one element for each element of `lambdas`.
#' @export
#' @examples
#' # Below we use the dismo::maxent example to fit a Maxent model:
#' if (require(dismo) && require(rJava) && 
#'     file.exists(system.file('java/maxent.jar', package='dismo'))) {
#'   fnames <- list.files(system.file('ex', package='dismo'), '\\.grd$', 
#'                        full.names=TRUE )
#'   predictors <- stack(fnames)
#'   occurrence <- system.file('ex/bradypus.csv', package='dismo')
#'   occ <- read.table(occurrence, header=TRUE, sep=',')[,-1]
#'   me <- maxent(predictors, occ, args=c('hinge=false', 'threshold=false'),
#'                path=file.path(tempdir(), 'ex1'))
#'   n_features(me)
#'   
#'   # comparing multiple models
#'   me2 <- maxent(predictors, occ, args=c('betamultiplier=3'), 
#'                 path=file.path(tempdir(), 'ex2'))
#'   n_features(list(me, me2))
#'   
#'   # or use lambdas files
#'   ff <- file.path(tempdir(), c('ex1', 'ex2'), 'species.lambdas')
#'   ff
#'   n_features(ff)
#' }
n_features <- function(lambdas) {
  if(length(lambdas)==1) {
    lambdas <- list(parse_lambdas(lambdas)$lambdas)
  } else {
    lambdas <- lapply(lambdas, function(x) parse_lambdas(x)$lambdas)  
  }
  k <- sapply(lambdas, function(x) sum(x$lambda != 0))
  k
}
