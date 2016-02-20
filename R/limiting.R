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
#'   Elith, J., Kearney, M. and Phillips, S. (2010), \href{http://onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2010.00036.x/abstract}{The art of modelling range-shifting species.} \emph{Methods in Ecology and Evolution}, 1: 330–342. doi: 10.1111/j.2041-210X.2010.00036.x
#' @importFrom raster stack stackApply as.factor levels
#' @export
#' @examples
#' library(dismo)
#' library(rJava)
#' 
#' fnames <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
#'                      pattern='grd', full.names=TRUE )
#' predictors <- stack(fnames)
#' occurence <- paste(system.file(package='dismo'), '/ex/bradypus.csv', sep='')
#' occ <- read.table(occurence, header=TRUE, sep=',')[,-1]
#' me <- maxent(predictors, occ, factors='biome')
#' 
#' limiting(predictors, me)
limiting <- function(x, me) {
  best <- lapply(seq_along(me@presence), function(i) {
    if(any(subset(parse_lambdas(me)$lambdas, 
                  var==names(me@presence)[i])$type=='categorical')) {
      as.numeric(names(which.max(table(me@presence[, i]))), 
                 levels=levels(me@presence[, i]))
    } else {
      mean(me@presence[, i]) # should vectorise with colMeans, but then need to sort out categoricals
    }
  })
  L <- lapply(seq_along(names(me@presence)), function(i) {
    p <- x[[names(me@presence)]]
    p[[i]][] <- best[[i]]
    project_maxent(me, p, quiet=TRUE)$prediction_logistic
  })
  pred <- project_maxent(me, x, quiet=TRUE)$prediction_logistic
  limiting <- which.max(stack(L)- pred)
  limiting <- raster::as.factor(limiting)
  lev <- raster::levels(limiting)[[1]]
  lev$predictor <- names(me@presence)[lev$ID]
  levels(limiting) <- lev
  limiting
}
