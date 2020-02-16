#' Transform raw or cloglog values to logistic values
#'
#' Transform Maxent raw or cloglog (complementary log-log) values to logistic
#' values.
#' 
#' @param x The values to Transform. This can be, e.g., a numeric atomic vector,
#'   numeric matrix or `Raster*` object.
#' @param from A character string indicating the type of the source values. This
#'   can be one of 'raw', 'cloglog', or 'logistic'.
#' @param H The entropy of the Maxent distribution. This is only used when
#'   `from='raw'`.
#' @keywords maxent, lambdas
#' @seealso [to_cloglog()]
#' @export
to_logistic <- function(x, from, H) {
  if(from=='raw' & missing(H)) 
    stop('When from is "raw", H must not be missing.')
  switch(from, 
         'cloglog'=1/(log(1-x) - 1) + 1,
         'raw'=(exp(H)*x)/(1+(exp(H)*x)),
         'logistic'=x)
}
