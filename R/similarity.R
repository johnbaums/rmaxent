#' Calculate Multivariate Environmental Similarity
#'
#' Calculate Multivariate Environmental Similarity with respect to a reference
#' dataset, for a set of environmental variables.
#'
#' @param x a `Raster*`, `list`, `matrix`, or `data.frame`
#'   where each layer/column/element represents focal values of an environmental
#'   variable.
#' @param ref a `list`, `matrix`, or `data.frame` where each
#'   column/element represents reference values for an environmental variable
#'   (corresponding to those given in `x`).
#' @param full (logical) should similarity values be returned for all variables?
#'   If `FALSE` (the default), then only the minimum similarity scores
#'   across variables will be returned.
#' @return If `x` is a `Raster*` object, this function returns a
#'   `Raster` layer giving the corresponding multivariate environmental
#'   similarity grid. If `full` is `TRUE`, this object will be
#'   returned in a list along with a `Raster*` object giving the
#'   environmental similarity for each layer of `x`. If `x` is a
#'   `list`, `matrix`, or `data.frame`, the function will return
#'   a vector giving multivariate similarity, or, if `full` is `TRUE`,
#'   a `list` with a `matrix` giving the environmental similarity for
#'   each element/column of `x`, as well as the multivariate similarity
#'   vector.
#' @details `similarity` uses the MESS algorithm described in Appendix S3
#'   of Elith et al. 2010.
#' @keywords maxent, mess, similarity, environment
#' @references 
#' * Elith, J., Kearney, M., and Phillips, S. (2010) [The art of modelling range-shifting species](https://doi.org/10.1111/j.2041-210X.2010.00036.x). _Methods in Ecology and Evolution_, 1: 330-342. doi:10.1111/j.2041-210X.2010.00036.x
#' @importFrom raster stack nlayers init as.data.frame raster
#' @importFrom methods is
#' @importFrom stats na.omit
#' @export
#' @examples
#' library(dismo)
#' library(raster)
#' ff <- list.files(system.file('ex', package='dismo'), '\\.grd$', 
#'                  full.names=TRUE )
#' predictors <- stack(grep('biome', ff, value=TRUE, invert=TRUE))
#' occ <- read.csv(system.file('ex/bradypus.csv', package='dismo'))[, -1]
#' ref <- extract(predictors, occ)
#' mess <- similarity(predictors, ref, full=TRUE)
similarity <- function(x, ref, full=FALSE) {
  if(!methods::is(ref, 'data.frame')) {
    ref <- as.data.frame(ref)
  }
  if(is(x, 'Raster')) {
    r <- TRUE
    if(isTRUE(full)) {
      out <- raster::stack(replicate(
        raster::nlayers(x), raster::init(x, function(x) NA)))
    } else {
      out <- raster::init(x, function(x) NA)
    }
  } else r <- FALSE
  ref <- stats::na.omit(ref)
  if(!methods::is(x, 'data.frame')) {
    x <- as.data.frame(x)
  }
  if(is.null(dim(ref))) {
    rng <- as.data.frame(range(ref, na.rm=TRUE))
  } else {
    rng <- as.data.frame(apply(ref, 2, range, na.rm=TRUE))
  }
  pct_less <- mapply(function(x, ref) {
    findInterval(x, sort(ref))/length(ref)
  }, x, ref, SIMPLIFY=FALSE)
  sim <- mapply(function(f, rng, p) {
    ifelse(f==0, (p-rng[1])/diff(rng)*100,
           ifelse(f > 0 & f <= 0.5, f*200,
                  ifelse(f > 0.5 & f < 1, (1-f)*200,
                         (rng[2]-p)/diff(rng)*100)))
  }, pct_less, rng, x)
  min_sim <- if(is.matrix(sim)) apply(sim, 1, min) else(min(sim))
  if(isTRUE(r)) {
    out_min <- raster::raster(out)
    out_min[] <- min_sim
    if(isTRUE(full)) {
      out[] <- sim
      list(similarity=out, similarity_min=out_min)
    } else out_min
  } else {
    if(isTRUE(full)) {
      list(similarity=sim, similarity_min=min_sim)
    } else min_sim
  }
}
