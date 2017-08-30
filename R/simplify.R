#' Identify an uncorrelated, useful subset of Maxent predictors
#'
#' Given a candidate set of predictor variables, this function identifies a
#' subset that meets specified multicollinearity criteria. Subsequently,
#' backward stepwise variable selection is used to iteratively drop the variable
#' that contributes least to the model, until the contribution of each variable 
#' meets a specified minimum, or until a predetermined minimum number of 
#' predictors remains.
#' 
#' @param occ A \code{data.frame} with predictor values for presence localities,
#'   where columns are predictors, and rows are samples. Additionally, a column
#'   indicating the species being modelled, which should have column name as
#'   specified by argument \code{species_column} (can be a character string, 
#'   numeric ID, etc.). The set of values given in this species indicator column
#'   must be identical to the set given in the corresponding column of
#'   \code{bg}.
#' @param bg A \code{data.frame} with predictor values for background
#'   localities, where columns are predictors, and rows are samples.
#'   Additionally, a column indicating the species being modelled, which should
#'   have column name as specified by argument \code{species_column} (can be a
#'   character string, numeric ID, etc.). The set of values given in this
#'   species indicator column must be identical to the set given in the
#'   corresponding column of \code{occ}.
#' @param path The output path within which output subdirectories will be 
#'   created for each species given in the column of \code{occ} (and \code{bg}) 
#'   specified by \code{species_column}. If missing, a temporary directory will
#'   be used.
#' @param species_column The column of \code{occ} (and \code{bg}) that contains 
#'   values indicating which species the samples belong to (e.g., this might be
#'   species name, or species ID).
#' @param response_curves Logical value indicating whether response curves
#'   should be included in Maxent model html output.
#' @param logistic_format Logical value indicating whether maxentResults.csv 
#'   should report logistic value thresholds (\code{TRUE}) or cloglog value
#'   thresholds (\code{FALSE}). This has no effect for versions of Maxent prior 
#'   to 3.4.
#' @param type The variable contribution metric to use when dropping variables.
#'   This can be \code{'PC'} (percent contribution) or \code{'PI'} (permutation
#'   importance; the default). See the Maxent tutorial for additional details.
#' @param cor_thr The maximum allowable pairwise correlation between predictor
#'   variables (calculated across presence and background localities). 
#' @param pct_thr The minimum allowable percent variable contribution (where 
#'   contribution type is specified by \code{type}). This should be specified as
#'   a value between 0 and 100.
#' @param k_thr The minimum number of variables to be kept in the model.
#' @param quiet Logical value indicating whether progress messages should be 
#'   suppressed (\code{TRUE}) or printed (\code{FALSE}).
#' @return The final fitted \code{MaxEnt} object.
#' @details If \code{path} is provided, subdirectories will be created within 
#'   \code{path}, with names equal to the values provided in the 
#'   \code{species_column} column of \code{occ}. Within these species
#'   subdirectories, two additional directories will be created: "full" contains
#'   the Maxent output corresponding to the model using the full uncorrelated
#'   subset of variables, while "final" contains the Maxent output corresponding
#'   to the model fit with the subset of those variables that each contribute at
#'   least \code{pct_thr}% to the model. Additionally, the \code{MaxEnt} R
#'   objects for the full and final fitted models are saved into these 
#'   directories, each with the name "model.rds". These can be read back into R
#'   with \code{\link{readRDS}}.
#' @keywords maxent, variable selection, correlation
#' @importFrom usdm vifcor
#' @importFrom dismo maxent
#' @export
#' @examples
#' # Below we modify the example given at ?dismo::maxent:
#' if (require(dismo) && require(rJava) &&
#'     file.exists(system.file('java/maxent.jar', package='dismo'))) {
#'   fnames <- list.files(system.file('ex', package="dismo"), '\\.grd$', 
#'                        full.names=TRUE)
#'   fnames <- grep('biome', fnames, value=TRUE, invert=TRUE)
#'   predictors <- scale(stack(fnames))
#'   occurrence <- system.file('ex/bradypus.csv', package='dismo')
#'   occ <- read.table(occurrence, header=TRUE, sep=',')[,-1]
#'   bg <- xyFromCell(predictors, Which(!is.na(sum(predictors)), cells=TRUE))
#'   occ_swd <- data.frame(species='bradypus', extract(predictors, occ))
#'   bg_swd <- data.frame(species='bradypus', extract(predictors, bg))
#'   m <- simplify(occ_swd, bg_swd, cor_thr=0.7, pct_thr=5, k_thr=4, quiet=FALSE)
#' }
simplify <- function(
  occ, bg, path, species_column='species', response_curves=FALSE,
  logistic_format=TRUE, type='PI', cor_thr, pct_thr, k_thr, quiet=TRUE) {
  if(missing(path)) {
    save <- FALSE
    path <- tempdir()
  } else save <- TRUE
  occ_by_species <- split(occ, occ[[species_column]])
  bg_by_species <- split(bg, bg[[species_column]])
  if(!identical(sort(names(occ_by_species)), sort(names(bg_by_species)))) {
    stop('The same set of species names must exist in occ and bg')
  }
  type <- switch(type, 
                 'PI'='permutation.importance',
                 'PC'='contribution',
                 stop('type must be either "PI" or "PC".', call.=FALSE))
  
  args <- c('threshold=false', 'hinge=false')
  if(isTRUE(response_curves)) args <- c(args, 'responsecurves=TRUE')
  if(isTRUE(logistic_format)) args <- c(args, 'outputformat=logistic')
  
  lapply(names(occ_by_species), function(name) {
    if(!quiet) message('\n\nDoing ', name)
    name_ <- gsub(' ', '_', name)
    swd <- rbind(occ_by_species[[name]], bg_by_species[[name]])
    swd <- swd[, -match(species_column, names(swd))]
    if(ncol(swd) < k_thr) stop('Initial number of variables < k_thr')
    round(cor(swd, use='pairwise'), 2)
    pa <- rep(1:0, c(nrow(occ_by_species[[name]]), nrow(bg_by_species[[name]])))
    ok <- as.character(
      usdm::vifcor(swd, maxobservations=nrow(swd), th=cor_thr)@results$Variables)
    swd_uncor <- swd[, ok]
    d <- file.path(path, name_, 'full')
    m <- dismo::maxent(swd_uncor, pa, args=args, path=d)
    if(isTRUE(save)) saveRDS(m, file.path(d, 'model.rds'))
    
    pct <- m@results[grep(type, rownames(m@results)), ]
    pct <- sort(pct[pct > 0])
    names(pct) <- sub(paste0('\\.', type), '', names(pct))
    if(min(pct) >= pct_thr || length(pct) <= k_thr) {
      if(isTRUE(save)) {
        d_out <- file.path(path, name_, 'final')
        dir.create(d_out)
        file.copy(list.files(d, full.names=TRUE), 
                  d_out, recursive=TRUE)
      }
      return(m)
    }
    while(min(pct) < pct_thr && length(pct) > k_thr) {
      message('Dropping ', names(pct)[1])
      swd_uncor <- swd_uncor[, -match(names(pct)[1], names(swd_uncor))]
      tmp <- tempfile()
      if(!quiet) message(
        sprintf('%s variables: %s', ncol(swd_uncor), 
                paste0(colnames(swd_uncor), collapse=', ')))
      m <- dismo::maxent(swd_uncor, pa, args=args, path=tmp)
      pct <- m@results[grep(type, rownames(m@results)), ]
      pct <- sort(pct)
      names(pct) <- sub(paste0('\\.', type), '', names(pct))
    }
    if(isTRUE(save)) {
      d_out <- file.path(path, name_, 'final')
      file.copy(tmp, file.path(path, name_), recursive=TRUE)
      file.rename(file.path(path, name_, basename(tmp)), d_out)
      saveRDS(m, file.path(path, name_, 'final/model.rds')) 
    }
    return(m)
  })
}
