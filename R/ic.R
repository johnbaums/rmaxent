#' Calculate information criteria for Maxent models
#'
#' Calculate AIC, AICc, and BIC for Maxent models as implemented in ENMTools.
#'
#' @param x Either a file path to a Maxent raw prediction raster, or a
#'   \code{Raster} or \code{RasterStack} containing raw Maxent predictions for
#'   which information criteria will be calculated.
#' @param occ A \code{matrix} or \code{data.frame} containing the coordinates
#'   for occurrence data used to fit the model.
#' @param lambdas The path to a Maxent .lambdas file, or a \code{MaxEnt} fitted
#'   model object.
#' @return An n x 6 matrix, where n is the number of Maxent rasters for which
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
#'   \item{Warren, D. L., Glor, R. E. and Turelli, M. 2009. \href{http://onlinelibrary.wiley.com/doi/10.1111/j.1600-0587.2009.06142.x/full}{ENMTools: a toolbox for comparative studies of environmental niche models.} \emph{Ecography} 33:607-611}
#'   \item{\href{http://enmtools.blogspot.com.au/}{ENMTools}}
#' }
#' @importFrom raster raster cellFromXY extract
#' @export
#' @examples
#' library(dismo)
#' fnames <- list.files(path=file.path(system.file(package='dismo'), 'ex'), 
#'                      patt='grd', full=TRUE )
#' predictors <- stack(fnames)
#' occ <- read.csv(file.path(system.file(package="dismo"), 'ex/bradypus.csv'))[, -1]
#' me <- maxent(predictors, occ, args=c('hinge=false', 'threshold=false'),
#'              path=tempdir())
#' r <- project_maxent(me, predictors, quiet=TRUE)$prediction_raw
#' 
#' # passing the raster object to pred.raw and the maxent object to lambdas:
#' maxent_IC(r, occ, me)
#' 
#' # passing a lambdas file path to lambdas:
#' maxent_IC(r, occ, file.path(tempdir(), 'species.lambdas'))
#' 
#' # passing a raster file path and lambdas file path to lambdas:
#' writeRaster(r, f <- tempfile(fileext='.tif'))
#' maxent_IC(f, occ, file.path(tempdir(), 'species.lambdas'))
ic <- function(x, occ, lambdas) {
  if(is.character(x)) x <- raster::raster(x)
  lambdas <- parse_lambdas(lambdas)$lambdas
  if(any(lambdas$type %in% c('threshold', 'hinge'))) 
    stop('Cannot calculate information criteria when threshold or hinge',
         ' features are in use.',
         call.=FALSE)
  k <- sum(lambdas$lambda != 0)
  occ <- occ[!duplicated(cellFromXY(x, occ)), ]
  n <- nrow(occ)
  if(k > n) {
    warning('Number of parameters greater than number of occurrence points for ',
            substitute(x), '. Information criteria not calculated.')
    out <- 
    return(c(n=n, k=k, ll=NA, AIC=NA, AICc=NA, BIC=NA))
  }
  out <- t(sapply(seq_len(nlayers(x)), function(i) {
    x <- x[[i]]
    x <- x/sum(values(x), na.rm=TRUE)
    ll <- sum(log(extract(x, occ)))
    AIC <- 2*k - 2*ll
    AICc <- AIC + ((2*k*(k+1))/(n - k - 1))
    BIC <- k*log(n) - 2*ll
    c(n=n, k=k, ll=ll, AIC=AIC, AICc=AICc, BIC=BIC)
  }))
  row.names(out) <- names(x)
  out
}
