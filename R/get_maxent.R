#' Download maxent.jar into dismo/java
#'
#' `get_maxent` downloads maxent.jar.
#'
#' This downloads the java file for the requested Maxent version into the
#' appropriate dismo/java folder.
#'
#' @param version Either 'latest' for the latest stable release, or a character
#'   version string indicating the version to download. See
#'   `maxent_versions` for available versions. Default is `'latest'`.
#' @param quiet logical. Should success message be suppressed?
#' @return Invisibly returns `0` for success, and a non-zero integer
#'   otherwise.
#' @seealso [maxent_versions()]
#' @importFrom utils unzip
#' @importFrom httr set_config config GET write_disk stop_for_status
#' @export
#' @details See references for license and citation details.
#' @references 
#' * GitHub repository: [mrmaxent/Maxent](https://github.com/mrmaxent/Maxent)
#' * Steven J. Phillips, Miroslav Dudík, Robert E. Schapire. [Maxent software for modeling species niches and distributions.](https://biodiversityinformatics.amnh.org/open_source/maxent/)
#' @examples
#' \dontrun{
#' get_maxent('3.3.3k')
#' }
get_maxent <- function(version='latest', quiet=FALSE) {
  vv <- maxent_versions(include_beta=TRUE)
  latest <- rev(grep('[Bb]eta', vv, invert=TRUE, value=TRUE))[1]
  if(version=='latest') version <- latest
  if(!version %in% vv)
    stop('Version "', version, '" unavailable.\n',
         'Use maxent_versions() to show accepted version strings',
         call.=FALSE)
  d <- system.file(package='dismo', 'java')
  if(version==latest) {
    cfg <- httr::set_config(httr::config(ssl_verifypeer = 0L))
    on.exit(httr::set_config(httr::config(cfg)))
    u <- 'https://biodiversityinformatics.amnh.org/open_source/maxent/maxent.php?op=download'
    resp <- httr::GET(u, httr::write_disk(f <- tempfile(fileext='.zip')))
    httr::stop_for_status(resp)
    utils::unzip(f, exdir=d, files='maxent.jar', junkpaths=TRUE)
  } else {
    u <- sprintf('%s/blob/master/ArchivedReleases/%s/maxent.jar?raw=true',
                 'https://github.com/mrmaxent/Maxent', version)
    resp <- httr::GET(u, httr::write_disk(file.path(d, 'maxent.jar'), 
                                            overwrite=TRUE))
    httr::stop_for_status(resp)
  }
  if(!quiet) message('Maxent version ', version, ' downloaded.')
  invisible(httr::status_code(resp)!=200)
}
