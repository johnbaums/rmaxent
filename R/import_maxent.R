#' Import a Maxent model into a MaxEnt object
#'
#' Create a MaxEnt object from a previously-fitted Maxent model
#' 
#' @param dir The file path to the directory containing Maxent output.
#' @param prefix The species name, which forms the prefix of the .lambdas file. 
#'   The default is \code{'species'} (and so species.lambdas is expected to
#'   exist in \code{dir}).
#' @return A \code{MaxEnt} object (see \code{\link{MaxEnt-class}}).
#' @importClassesFrom dismo MaxEnt
#' @importFrom methods new
#' @importFrom utils read.csv
#' @export
#' @examples
#' library(dismo)
#' predictors <- stack(list.files(
#'   file.path(system.file(package='dismo'), 'ex'), '\\.grd$', full=TRUE))
#' occ <- read.csv(file.path(system.file(package="dismo"), 'ex/bradypus.csv'))[, -1]
#' d <- file.path(tempdir(), 'demo')
#' dir.create(d)
#' maxent(predictors, occ, factors='biome', path=d)
#' m <- import(d)
import_maxent <- function(dir, prefix='species') {
  m <- methods::new('MaxEnt',
           lambdas=readLines(file.path(dir, paste0(prefix, '.lambdas'))),
           results=t(utils::read.csv(file.path(dir, 'maxentResults.csv'))[, -1]),
           path=dir,
           html=file.path(dir, 'maxent.html'),
           hasabsence=file.exists(file.path(dir, 'absence')))
  if(m@hasabsence) m@absence <- read.csv(file.path(dir, 'absence'))[, -(1:3)]
  if(file.exists(file.path(dir, 'presence')))
    m@presence <- utils::read.csv(file.path(dir, 'presence'))[, -(1:3)]
  m
}
