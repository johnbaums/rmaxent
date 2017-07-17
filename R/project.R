#' Project a fitted Maxent model
#'
#' Project a fitted Maxent model by predicting to new environmental data.
#'
#' @param lambdas Either (1) a \code{MaxEnt} fitted model object (fitted with 
#'   the \code{maxent} function in the \code{dismo} package), (2) a file path to
#'   a Maxent .lambdas file, or (3) a \code{lambdas} object returned by 
#'   \code{\link{parse_lambdas}}.
#' @param newdata A \code{RasterStack}, \code{RasterBrick}, \code{list},
#'   \code{data.frame}, \code{data.table}, or \code{matrix} that has
#'   layers/elements/columns whose names correspond to the names of predictors
#'   used to fit the model. These layers/elements/columns must all have the same
#'   length.
#' @param return_lfx Logical. Should \code{Raster} layers be returned giving
#'   lambda*feature values for each feature with a non-zero lambda? Currently
#'   ignored if \code{newdata} is not a \code{Raster*} object.
#' @param mask (Optional; requires that \code{newdata} is a \code{Raster*} 
#'   object.) A \code{Raster} object with \code{NA} values in cells for which 
#'   the model should \emph{not} be projected. These cells will be assigned
#'   \code{NA} in the returned output.
#' @param quiet Logical. Should projection progress be reported?   
#' @return If \code{newdata} is a \code{RasterStack} or \code{RasterBrick}, a 
#'   list with three elements: 
#'  \itemize{
#'   \item{\code{prediction_raw}}{: a \code{Raster} layer giving the raw Maxent
#'   prediction; and}
#'   \item{\code{prediction_logistic}}{: a \code{Raster} layer giving the
#'   logistic Maxent prediction.}
#'   \item{\code{prediction_cloglog}}{: a \code{Raster} layer giving the
#'   cloglog Maxent prediction.}
#' }
#' If \code{newdata} is \emph{not} a \code{RasterStack} or \code{RasterBrick},
#' the raster layers will be replaced with \code{data.table}s in the returned
#' list.
#'
#' Additionally, if \code{newdata} is a \code{RasterStack} or \code{RasterBrick}
#' and \code{return_lfx} is \code{TRUE}, the returned list will include 
#' \code{prediction_lfx} (the logit scores for the linear predictor), and 
#' \code{lfx_all} (the contributions to \code{prediction_lfx} of each feature
#' with a non-zero lambda).
#' @details \code{project} uses feature weights described in a .lambas
#'   file or \code{MaxEnt} object to predict a Maxent model to environmental
#'   data. This function performs the projection entirely in R, without the need
#'   for the Maxent Java software. For tested datasets, it performs the 
#'   projection in roughly one third of the time taken for the same projection
#'   by maxent.jar.
#' @section Warning:
#' This function is still in development, and no guarantee is made for the
#' accuracy of its projections.
#' @keywords maxent, predict, project
#' @references 
#' \itemize{
#'   \item{Wilson, P. W. (2009) \href{http://gsp.humboldt.edu/OLM/GSP_570/Learning Modules/10 BlueSpray_Maxent_Uncertinaty/MaxEnt lambda files.pdf}{\emph{Guidelines for computing MaxEnt model output values from a lambdas file}}.}
#'   \item{\emph{Maxent software for species habitat modeling, version 3.3.3k} help file (software freely available \href{https://www.cs.princeton.edu/~schapire/maxent/}{here}).}
#' }
#' @seealso \code{\link{read_mxe}}
#' @importFrom raster raster mask compareRaster as.data.frame
#' @importFrom data.table data.table as.data.table is.data.table :=
#' @importFrom methods is
#' @importFrom stats complete.cases
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
#'   me <- maxent(predictors, occ, factors='biome')
#' 
#'   # ... and then predict it to the full environmental grids:
#'   pred <- project(me, predictors)
#'   # This is equivalent to using the predict method for MaxEnt objects:
#'   pred2 <- predict(me, predictors, args='outputformat=logistic')
#'   pred3 <- predict(me, predictors, args='outputformat=cloglog')
#'   
#'   all.equal(values(pred$prediction_logistic), values(pred2))
#'   all.equal(values(pred$prediction_cloglog), values(pred3))
#' }
project <- function(lambdas, newdata, return_lfx=FALSE, mask, quiet=FALSE) {
  if(!missing(mask)) {
    if(!methods::is(mask, 'RasterLayer')) {
      stop('mask should be a RasterLayer object')
    } else {
      if(!methods::is(newdata, 'Raster')) {
        stop('If mask is provided, newdata should be a Raster object with the',
             'same dimensions, extent, and CRS.')
      }
      if(!raster::compareRaster(mask, newdata, stopiffalse=FALSE))
        stop('If mask is provided, newdata should be a Raster object with the',
             'same dimensions, extent, and CRS.')
    }
  }
  if(!rmaxent:::is.lambdas(lambdas)) lambdas <- parse_lambdas(lambdas)
  meta <- lambdas[-1]
  lambdas <- lambdas[[1]]
  is_cat <- unique(
    gsub('\\(|==.*\\)', '', lambdas[lambdas$type=='categorical', 'feature']))
  nms <- unique(unlist(strsplit(lambdas$var[lambdas$lambda != 0], ',')))
  clamp_limits <- data.table::data.table(lambdas[lambdas$type=='linear', ])
  lambdas <- lambdas[lambdas$lambda != 0, ]
  lambdas <- split(lambdas, c('other', 'hinge')[
    grepl('hinge', lambdas$type) + 1])
  if (is(newdata, 'RasterStack') | is(newdata, 'RasterBrick') | is(newdata, 'Raster')) {
    pred_raw <- pred_logistic <- pred_cloglog <- pred_lfx <- raster::raster(newdata)
    if(!missing(mask)) {
      newdata <- raster::mask(newdata, mask)
    }
    newdata <- data.table::as.data.table(newdata[])
  }
  if (is.matrix(newdata)) newdata <- data.table::as.data.table(newdata)
  if (is.list(newdata) & !is.data.frame(newdata)) {
    if(length(unique(sapply(newdata, length))) != 1) 
      stop('newdata was provided as a list, but its elements vary in length.')
    newdata <- data.table::as.data.table(newdata)
  }
  if (!is.data.frame(newdata))
    stop('newdata must be a list, data.table, data.frame, matrix, RasterStack,',
         'or RasterBrick.')
  if (!data.table::is.data.table(newdata)) 
    newdata <- data.table::as.data.table(newdata)
  if(!all(nms %in% names(newdata))) {
    stop(sprintf('Variables missing in newdata: %s', 
                 paste(setdiff(nms, colnames(newdata)), collapse=', ')))
  }
  if (any(!names(newdata) %in% nms)) {
    newdata <- newdata[, setdiff(names(newdata), nms) := NULL]
  }
  na <- !complete.cases(newdata)
  newdata <- newdata[!na]
  
  # Clamp newdata
  invisible(lapply(setdiff(names(newdata), is_cat), function(x) {
    clamp_max <- clamp_limits[clamp_limits$feature==x, max]
    clamp_min <- clamp_limits[clamp_limits$feature==x, min]
    newdata[, c(x) := pmax(pmin(get(x), clamp_max), clamp_min)]
  }))
  
  k <- sum(sapply(lambdas, nrow))
  txt <- sprintf('\rCalculating contribution of feature %%%1$dd of %%%1$dd', 
                 nchar(k))
  lfx <- numeric(nrow(newdata))
  lfx_all <- setNames(vector('list', sum(sapply(lambdas, nrow))),
                      unlist(lapply(lambdas[2:1], function(x) x$feature)))
  
  if('other' %in% names(lambdas)) {
    for (i in seq_len(nrow(lambdas$other))) {
      if(!quiet) cat(sprintf(txt, i, k))
      x <- with(newdata, eval(parse(text=lambdas$other$feature[i])))
      # clamp feature
      x <- pmin(pmax(x, lambdas$other$min[i]), lambdas$other$max[i])
      x01 <- (x - lambdas$other$min[i]) / 
        (lambdas$other$max[i] - lambdas$other$min[i])
      lfx_all[[i]] <- lambdas$other$lambda[i] * x01
      lfx <- lfx + lfx_all[[i]]
    }
    rm(x, x01)
  }
  
  if('hinge' %in% names(lambdas)) {
    for (i in seq_len(nrow(lambdas$hinge))) {
      if(!quiet) cat(sprintf(txt, nrow(lambdas$other)+i, k))
      x <- with(newdata, get(sub("'|`", "", lambdas$hinge$feature[i])))
      x01 <- (x - lambdas$hinge$min[i]) / (lambdas$hinge$max[i] - lambdas$hinge$min[i])
      if (lambdas$hinge$type[i]=='reverse_hinge') {
        lfx_all[[nrow(lambdas$other) + i]] <- 
          (x < lambdas$hinge$max[i]) * lambdas$hinge$lambda[i] * (1-x01)
      } else {
        lfx_all[[nrow(lambdas$other) + i]] <- 
          (x >= lambdas$hinge$min[i]) * lambdas$hinge$lambda[i] * x01
      }
      lfx <- lfx + lfx_all[[nrow(lambdas$other) + i]]
    }
    rm(x, x01)
  }
  
  ln_raw <- lfx - meta$linearPredictorNormalizer - log(meta$densityNormalizer)
  raw <- exp(ln_raw)
  logit <- meta$entropy + ln_raw
  cloglog <- 1 - exp(-exp(meta$entropy) * raw)
  logistic <- plogis(logit)
  
  #linpred <- rep(NA_real_, length(na))
  #linpred[!na] <- lfx
  if(exists('pred_raw', inherits=FALSE)) {
    pred_raw[which(!na)] <- raw
    pred_logistic[which(!na)] <- logistic
    pred_cloglog[which(!na)] <- cloglog
    pred_lfx[which(!na)] <- lfx
    lfx_each <- lapply(lfx_all, function(x) {
      r <- raster(pred_raw)
      r[which(!na)] <- x
      r
    })
    out <- list(prediction_raw=pred_raw,
                prediction_logistic=pred_logistic,
                prediction_cloglog=pred_cloglog)
    if(isTRUE(return_lfx)) {
      out <- c(out, 
               list(prediction_lfx=pred_lfx,
                    lfx_all=lfx_each))
    } 
  } else {
    prediction_raw <- prediction_logistic <- prediction_cloglog <-
      rep(NA_real_, length(na))
    prediction_raw[!na] <- raw
    prediction_logistic[!na] <- logistic
    prediction_cloglog[!na] <- cloglog
    #prediction_lfx[!na] <- lfx
    out <- list(prediction_raw=prediction_raw,
                prediction_logistic=prediction_logistic,
                prediction_cloglog=prediction_cloglog)
  }
  return(out)
}
