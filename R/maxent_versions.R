#' Show available versions of Maxent
#'
#' `maxent_versions` returns a vector of available version numbers.
#'
#' This obtains a vector of versions numbers for available versions of Maxent,
#' from https://github.com/mrmaxent/Maxent/tree/master/ArchivedReleases, as 
#' well as the latest version number from 
#' http://biodiversityinformatics.amnh.org/open_source/maxent.
#'
#' @param include_beta logical. Should beta versions be included?
#' @return Returns a `character` vector of version numbers.
#' @seealso [get_maxent()]
#' @importFrom xml2 read_html
#' @importFrom rvest html_nodes html_text
#' @importFrom magrittr %>%
#' @importFrom httr set_config config GET content
#' @export
#' @examples
#' \dontrun{
#' maxent_versions()
#' }
maxent_versions <- function(include_beta=FALSE) {
  u <- 'https://github.com/mrmaxent/Maxent/tree/master/ArchivedReleases'
  v <- xml2::read_html(u) %>%
    rvest::html_nodes(xpath = '//div[@role="rowheader"]//a') %>%
    rvest::html_text()
  v <- v[-1]
  cfg <- httr::set_config(httr::config(ssl_verifypeer = 0L))
  on.exit(httr::set_config(httr::config(cfg)))
  u2 <- 'https://biodiversityinformatics.amnh.org/open_source/maxent/index.html'
  v2 <- httr::GET(u2) %>% 
    httr::content(encoding='UTF-8') %>% 
    rvest::html_nodes(xpath='//*[@id="Form"]/h3[1]') %>% 
    rvest::html_text()
  v2 <- gsub('^\\s*Current\\s+version\\s*|\\s*$', '', v2)
  v <- c(v, v2)
  v <- sort(v[v!=''])
  if(!include_beta) grep('beta', v, invert=TRUE, value=TRUE) else v
}
