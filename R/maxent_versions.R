#' Show available versions of Maxent
#'
#' `maxent_versions` returns a vector of available version numbers.
#'
#' This obtains a vector of versions numbers for available versions of Maxent,
#' from https://github.com/mrmaxent/Maxent/tree/master/ArchivedReleases.
#'
#' @param include_beta logical. Should beta versions be included?
#' @return Returns a `character` vector of version numbers.
#' @seealso [get_maxent()]
#' @importFrom xml2 read_html
#' @importFrom rvest html_nodes html_text
#' @importFrom magrittr %>%
#' @importFrom httr set_config config 
#' @export
#' @examples
#' \dontrun{
#' maxent_versions()
#' }
maxent_versions <- function(include_beta=FALSE) {
  u <- 'https://github.com/mrmaxent/Maxent/tree/master/ArchivedReleases'
  v <- xml2::read_html(u) %>%
    rvest::html_nodes(css = '.react-directory-row-name-cell-large-screen .Link--primary') %>%
    rvest::html_text() %>%
    sort()
  if(!include_beta) grep('beta', v, invert = TRUE, value = TRUE) else v
}
