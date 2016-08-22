#' Import a Maxent model as a MaxEnt object
#'
#' Create a MaxEnt object from a previously-fitted Maxent model
#' 
#' @param dir The file path to the directory containing Maxent output.
#' @param lambdas The name of the .lambdas file representing the fitted model 
#'   (e.g. \code{'species.lambdas'}), excluding the containing path. If not
#'   provided, the first (alphabetically) .lambdas file found in \code{dir} will
#'   be used.
#' @param html The name of the .html file containing reported Maxent results 
#'   (e.g. \code{'maxent.html'}), excluding the containing path. If not 
#'   provided, the first (alphabetically) .html file found in \code{dir} will be
#'   used.
#' @return A \code{MaxEnt} object (see \code{\link{MaxEnt-class}}).
#' @importClassesFrom dismo MaxEnt
#' @importFrom methods new
#' @importFrom utils read.csv
#' @export
#' @examples
#' # Below we use the dismo::maxent example to fit a Maxent model:
#' library(dismo)
#' if (require(dismo) && require(rJava) && 
#'     file.exists(system.file('java/maxent.jar', package='dismo'))) {
#'   predictors <- stack(list.files(
#'     file.path(system.file(package='dismo'), 'ex'), '\\.grd$', full=TRUE))
#'   occ <- read.csv(file.path(system.file(package="dismo"), 'ex/bradypus.csv'))[, -1]
#'   d <- file.path(tempdir(), 'demo')
#'   dir.create(d)
#'   m1 <- maxent(predictors, occ, factors='biome', path=d)
#'   m2 <- import_maxent(d)
#' }
import_maxent <- function(dir, lambdas, html) {
  l <- ifelse(missing(lambdas), 
         list.files(dir, '\\.lambdas$', full.names=TRUE)[1], 
         file.path(dir, basename(lambdas)))
  h <- ifelse(missing(html), 
              list.files(dir, '\\.html$', full.names=TRUE)[1], 
              file.path(dir, basename(html)))
  m <- methods::new('MaxEnt',
           lambdas=readLines(l),
           results=t(utils::read.csv(file.path(dir, 'maxentResults.csv'))[, -1]),
           path=dir,
           html=h,
           hasabsence=file.exists(file.path(dir, 'absence')))
  m@absence <- read.csv(file.path(dir, 'absence'))[, -(1:3)]
  m@presence <- utils::read.csv(file.path(dir, 'presence'))[, -(1:3)]
  m
}
