#' Calculate information criteria for Maxent models
#'
#' Calculate AIC, AICc, and BIC for Maxent models as implemented in ENMTools.
#'
#' @param x Either a \code{vector} of one or more file paths to Maxent raw 
#'   prediction raster files, or a \code{Raster} or \code{RasterStack}
#'   containing raw Maxent predictions for which information criteria will be
#'   calculated.
#' @param occ A \code{matrix} or \code{data.frame} containing the coordinates
#'   for occurrence data used to fit the model. This is assumed to be common to 
#'   all models being compared.
#' @param lambdas A \code{vector} of one or more file paths to Maxent .lambdas 
#'   files, or a \code{list} of one or more \code{MaxEnt} fitted model objects.
#'   The length of this vector/list should equal the length (or number of 
#'   layers) of \code{x}.
#' @return An n x 6 matrix, where n is the number of Maxent models for which
#'   information criteria are to be calculated. Columns give \code{n} (the
#'   number of occurrence records used for model training), \code{k} (the number
#'   of features with non-zero weights), \code{ll} (the negative log likelihood
#'   of the model), \code{AIC}, \code{AICc}, and \code{BIC} (as calculated in
#'   ENMTools).
#' @section Warning: 
#' These information criteria should \emph{not} be calculated for models that
#' use hinge or threshold features because the number of predictors is not
#' estimated correctly.
#' @references 
#' \itemize{
#'   \item{Warren, D. L., Glor, R. E. and Turelli, M. 2009. \href{http://doi.org/10.1111/j.1600-0587.2009.06142.x}{ENMTools: a toolbox for comparative studies of environmental niche models.} \emph{Ecography} 33:607-611}
#'   \item{\href{http://enmtools.blogspot.com.au/}{ENMTools}}
#' }
#' @importFrom raster cellFromXY extract nlayers values stack
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
#'                path=tempdir())
#'   r <- project(me, predictors, quiet=TRUE)$prediction_raw
#' 
#'   # passing the raster object to pred.raw and the maxent object to lambdas:
#'   ic(r, occ, me)
#'   
#'   # passing a lambdas file path to lambdas:
#'   ic(r, occ, file.path(tempdir(), 'species.lambdas'))
#'   
#'   # passing a raster file path and lambdas file path to lambdas:
#'   writeRaster(r, f <- tempfile(fileext='.tif'))
#'   ic(f, occ, file.path(tempdir(), 'species.lambdas'))
#'   
#'   # comparing multiple models
#'   me2 <- maxent(predictors, occ, args=c('hinge=false', 'threshold=false',
#'                 'betamultiplier=3'), path=tempdir())
#'   r2 <- project(me2, predictors, quiet=TRUE)$prediction_raw               
#'   ic(stack(r, r2), occ, list(me, me2))
#' }
ic <- function(x, occ, lambdas) {
  x <- raster::stack(x)
  if(length(lambdas)==1) {
    lambdas <- list(parse_lambdas(lambdas)$lambdas)
  } else {
    lambdas <- lapply(lambdas, function(x) parse_lambdas(x)$lambdas)  
  }
  if(any(sapply(lambdas, '[[', 'type') %in% c('threshold', 'hinge'))) 
    stop('Cannot calculate information criteria when threshold or hinge',
         ' features are in use.',
         call.=FALSE)
  k <- sapply(lambdas, function(x) sum(x$lambda != 0))
  occ <- occ[!duplicated(cellFromXY(x, occ)), ]
  n <- nrow(occ)
  if(all(k > n)) 
    stop('Number of parameters exceeds number of occurrence points for.')
  if(any(k > n)) {
    warning('Number of parameters exceeds number of occurrence points for ',
            'model:', toString(which(k > n)), 
            '. Information criteria not calculated for those models.')
    return(c(n=n, k=k, ll=NA, AIC=NA, AICc=NA, BIC=NA))
  }
  out <- t(mapply(function(pred, lam, k) {
    pred <- pred/sum(raster::values(pred), na.rm=TRUE)
    vals <- na.omit(raster::extract(pred, occ))
    n <- length(vals)
    ll <- sum(log(vals))
    AIC <- 2*k - 2*ll
    AICc <- AIC + ((2*k*(k+1))/(n - k - 1))
    BIC <- k*log(n) - 2*ll
    c(n=n, k=k, ll=ll, AIC=AIC, AICc=AICc, BIC=BIC)
  }, unstack(x), lambdas, k))
  row.names(out) <- names(x)
  out
}
