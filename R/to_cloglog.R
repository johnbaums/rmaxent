#' Transform raw or logistic values to cloglog values
#'
#' Transform Maxent raw or logistic values to cloglog (complementary log-log)
#' values.
#'
#' @param x The values to transform. This can be, e.g., a numeric atomic vector,
#'   numeric matrix or \code{Raster*} object.
#' @param from A character string indicating the type of the source values. This
#'   can be one of 'raw', 'cloglog', or 'logistic'.
#' @param H The entropy of the Maxent distribution. This is only used when
#'   \code{from='raw'}.
#' @keywords maxent, lambdas
#' @seealso \code{\link{to_logistic}}
#' @export
to_cloglog <- function(x, from, H) {
  if(from=='raw' & missing(H)) 
    stop('When from is "raw", H must not be missing.')
  switch(from, 
         'cloglog'=x,
         'raw'=1-exp(-exp(H)*x),
         'logistic'=1-exp(-(x/(1-x))))
}
