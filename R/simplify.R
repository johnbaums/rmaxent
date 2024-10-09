#' Identify an uncorrelated, useful subset of Maxent predictors
#'
#' Given a candidate set of predictor variables, this function identifies a
#' subset that meets specified multicollinearity criteria. Subsequently,
#' backward stepwise variable selection is used to iteratively drop the variable
#' that contributes least to the model, until the contribution of each variable 
#' meets a specified minimum, or until a predetermined minimum number of 
#' predictors remains.
#' 
#' @param occ A `data.frame` with predictor values for presence localities,
#'   where columns are predictors, and rows are samples. Additionally, a column
#'   indicating the species being modelled, which should have column name as
#'   specified by argument `species_column` (can be a character string, 
#'   numeric ID, etc.). The set of values given in this species indicator column
#'   must be identical to the set given in the corresponding column of
#'   `bg`.
#' @param bg A `data.frame` with predictor values for background
#'   localities, where columns are predictors, and rows are samples.
#'   Additionally, a column indicating the species being modelled, which should
#'   have column name as specified by argument `species_column` (can be a
#'   character string, numeric ID, etc.). The set of values given in this
#'   species indicator column must be identical to the set given in the
#'   corresponding column of `occ`.
#' @param path The output path within which output subdirectories will be 
#'   created for each species given in the column of `occ` (and `bg`) 
#'   specified by `species_column`. If missing, a temporary directory will
#'   be used.
#' @param species_column The column of `occ` (and `bg`) that contains 
#'   values indicating which species the samples belong to (e.g., this might be
#'   species name, or species ID).
#' @param response_curves Logical value indicating whether response curves
#'   should be included in Maxent model html output.
#' @param logistic_format Logical value indicating whether maxentResults.csv 
#'   should report logistic value thresholds (`TRUE`) or cloglog value
#'   thresholds (`FALSE`). This has no effect for versions of Maxent prior 
#'   to 3.4.
#' @param type The variable contribution metric to use when dropping variables.
#'   This can be `'PC'` (percent contribution) or `'PI'` (permutation
#'   importance; the default). See the Maxent tutorial for additional details.
#' @param cor_thr The maximum allowable pairwise correlation between predictor
#'   variables (calculated across presence and background localities). 
#' @param pct_thr The minimum allowable percent variable contribution (where 
#'   contribution type is specified by `type`). This should be specified as
#'   a value between 0 and 100.
#' @param k_thr The minimum number of variables to be kept in the model.
#' @param features Features to include. Specify as a string comprising one or
#'   more of 'l' (linear), 'p' (product), 'q' (quadratic), 't' (threshold), and
#'   'h' (hinge). E.g., `features='lpq'` (equivalently,
#'   `features='plq'`). The default is `'lpq'`.
#' @param replicates The number of cross-validation replicates to perform. When
#'   cross-validation is used, the average (over folds) of the variable
#'   contribution metric is used.
#' @param quiet Logical value indicating whether progress messages should be 
#'   suppressed (`TRUE`) or printed (`FALSE`).
#' @return The final fitted `MaxEnt` object.
#' @details If `path` is provided, subdirectories will be created within 
#'   `path`, with names equal to the values provided in the 
#'   `species_column` column of `occ`. Within these species
#'   subdirectories, two additional directories will be created: "full" contains
#'   the Maxent output corresponding to the model using the full uncorrelated
#'   subset of variables, while "final" contains the Maxent output corresponding
#'   to the model fit with the subset of those variables that each contribute at
#'   least `pct_thr`% to the model. Additionally, the `MaxEnt` R
#'   objects for the full and final fitted models are saved into these 
#'   directories, each with the name "model.rds". These can be read back into R
#'   with [readRDS()].
#' @keywords maxent, variable selection, correlation
#' @importFrom usdm vifcor
#' @importFrom dismo maxent
#' @importFrom methods slot
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
  occ, bg, path, species_column='species', response_curves=TRUE,
  logistic_format=TRUE, type='PI', cor_thr, pct_thr, k_thr,
  features='lpq', replicates=1, quiet=TRUE) {
  if(!species_column %in% colnames(occ))
    stop(species_column, ' is not a column of `occ`', call.=FALSE)
  if(!species_column %in% colnames(bg))
    stop(species_column, ' is not a column of `bg`', call.=FALSE)
  if(missing(path)) {
    save <- FALSE
    path <- tempdir()
  } else save <- TRUE
  features <- unlist(strsplit(gsub('\\s', '', features), ''))
  if(length(setdiff(features, c('l', 'p', 'q', 'h', 't'))) > 1)
    stop("features must be a vector of one or more of ',
         'l', 'p', 'q', 'h', and 't'.")
  off <- setdiff(c('l', 'p', 'q', 't', 'h'), features)
  if(length(off) > 0) {
    off <- c(l='linear=FALSE', p='product=FALSE', q='quadratic=FALSE',
             t='threshold=FALSE', h='hinge=FALSE')[off]
  }
  off <- unname(off)
  occ_by_species <- split(occ, occ[[species_column]])
  bg_by_species <- split(bg, bg[[species_column]])
  if(!identical(sort(names(occ_by_species)), sort(names(bg_by_species)))) {
    stop('The same set of species names must exist in occ and bg')
  }
  type <- switch(type, 
                 'PI'='permutation.importance',
                 'PC'='contribution',
                 stop('type must be either "PI" or "PC".', call.=FALSE))
  
  args <- off
  if(replicates > 1) args <- c(args, paste0('replicates=', replicates))
  if(isTRUE(response_curves)) args <- c(args, 'responsecurves=TRUE')
  if(isTRUE(logistic_format)) args <- c(args, 'outputformat=logistic')
  f <- function(name) {
    if(!quiet) message('\n\nDoing ', name)
    name_ <- gsub(' ', '_', name)
    swd <- rbind(occ_by_species[[name]], bg_by_species[[name]])
    swd <- swd[, -match(species_column, names(swd))]
    if(ncol(swd) < k_thr) stop('Initial number of variables < k_thr', call.=FALSE)
    pa <- rep(1:0, c(nrow(occ_by_species[[name]]), nrow(bg_by_species[[name]])))
    vc <- usdm::vifcor(swd, size=nrow(swd), th=cor_thr)
    vif <- methods::slot(vc, 'results')
    k <- nrow(vif)
    exclude <- methods::slot(vc, 'excluded')
    if(!isTRUE(quiet) & length(exclude) > 0) {
      message('Dropped due to collinearity: ', paste0(exclude, collapse=', '))
    }
    if(k < k_thr)
      stop(sprintf('Number of uncorrelated variables (%s) < k_thr (%s). %s', 
                   k, k_thr, 
                   'Reduce k_thr, increase cor_thr, or find alternative predictors.'),
           call.=FALSE)
    swd_uncor <- swd[, as.character(vif$Variables)]
    d <- file.path(path, name_, if(replicates > 1) 'xval' else 'full')
    m <- dismo::maxent(swd_uncor, pa, args=args, path=d)
    if(isTRUE(save)) saveRDS(m, file.path(d, 'model.rds'))
    pct <- m@results[grep(type, rownames(m@results)), , drop=FALSE]
    pct <- pct[, ncol(pct)]
    pct <- sort(pct)
    names(pct) <- sub(paste0('\\.', type), '', names(pct))
    if(min(pct) >= pct_thr || length(pct) == k_thr) {
      if(replicates > 1) {
        d <- file.path(path, name_, 'full')
        m <- dismo::maxent(
          swd_uncor, pa, args=grep('replicates', args, value=TRUE, invert=TRUE), 
          path=d)
      }
      if(isTRUE(save)) {
        saveRDS(m, file.path(d, 'model.rds'))
      }
      return(m)
    }
    while(min(pct) < pct_thr && length(pct) > k_thr) {
      candidates <- vif[vif$Variables %in% names(pct)[pct==pct[1]], ]
      drop <- as.character(candidates$Variables[which.max(candidates$VIF)])
      message('Dropping ', drop)
      swd_uncor <- swd_uncor[, -match(drop, colnames(swd_uncor))]
      if(!quiet) message(
        sprintf('%s variables: %s', ncol(swd_uncor), 
                paste0(colnames(swd_uncor), collapse=', ')))
      m <- dismo::maxent(swd_uncor, pa, args=args, path=d)
      pct <- m@results[grep(type, rownames(m@results)), , drop=FALSE]
      pct <- pct[, ncol(pct)]
      pct <- sort(pct)
      names(pct) <- sub(paste0('\\.', type), '', names(pct))
    }
    if(replicates > 1) {
      d <- file.path(path, name_, 'full')
      m <- dismo::maxent(
        swd_uncor, pa, args=grep('replicates', args, value=TRUE, invert=TRUE), 
        path=d)
    }
    if(isTRUE(save)) {
      saveRDS(m, file.path(d, 'model.rds'))
    }
    return(m)
  }
  lapply(names(occ_by_species), f)
}
