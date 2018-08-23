#' Calculate which factors limit predicted habitat suitability
#'
#' Calculate a limiting factor surface corresponding to projections of a Maxent
#' model.
#'
#' @param x A \code{RasterStack} containing predictors grids for predictors used
#'   in the Maxent model specified at \code{me}.
#' @param me A a \code{MaxEnt} fitted model object (from \code{dismo::maxent}).
#' @return A \code{Raster} object indicating the most limiting factor at each
#'   grid cell. This is defined, for each grid cell, as the variable whose value
#'   is most responsible for decreasing suitability at that cell. The decrease 
#'   in suitability is calculated, for each predictor in turn, relative to the 
#'   suitability (logistic prediction) that would be achieved if that predictor
#'   took the value equal to the mean at occurrence sites (median for
#'   categorical variables). The predictor associated with the largest decrease
#'   in suitability is the most limiting factor.
#' @references 
#'   Elith, J., Kearney, M. and Phillips, S. (2010), \href{http://doi.org/10.1111/j.2041-210X.2010.00036.x}{The art of modelling range-shifting species.} \emph{Methods in Ecology and Evolution}, 1: 330-342. doi: 10.1111/j.2041-210X.2010.00036.x
#' @importFrom raster stack stackApply as.factor levels which.max
#' @export
#' @examples
#' # This uses the ?dismo::maxent example. Only run if maxent.jar exists.
#' if(require(dismo) && require(rJava) && 
#'    file.exists(system.file('java/maxent.jar', package='dismo'))) {
#'   fnames <- list.files(system.file('ex', package='dismo'), '\\.grd$', 
#'                        full.names=TRUE)
#'   predictors <- stack(fnames)
#'   occurrence <- system.file('ex/bradypus.csv', package='dismo')
#'   occ <- read.table(occurrence, header=TRUE, sep=',')[,-1]
#'   me <- maxent(predictors, occ, factors='biome')
#'   limiting(predictors, me)
#' }
limiting <- function(x, me) {
  lam <- parse_lambdas(me)$lambdas
  best <- lapply(seq_along(me@presence), function(i) {
    if(any(lam[lam$var==names(me@presence)[i], ]$type=='categorical')) {
      as.numeric(names(which.max(table(me@presence[, i]))), 
                 levels=levels(me@presence[, i]))
    } else {
      mean(me@presence[, i])
    }
  })
  L <- lapply(seq_along(names(me@presence)), function(i) {
    p <- x[[names(me@presence)]]
    p[[i]][] <- best[[i]]
    project(me, p, quiet=TRUE)$prediction_logistic
  })
  pred <- project(me, x, quiet=TRUE)$prediction_logistic
  limiting <- raster::which.max(raster::stack(L)- pred)
  limiting <- raster::as.factor(limiting)
  lev <- raster::levels(limiting)[[1]]
  lev$predictor <- names(me@presence)[lev$ID]
  levels(limiting) <- lev
  limiting
}
