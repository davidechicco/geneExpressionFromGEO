
#' Function that returns numeric values with 2 decimal numbers.
#' 
#' @param x input numeric value with N decimal numbers.
#' @return a numeric value with 2 decimal numbers.
#' @examples
#' aaa <- dec_two(8.31232)
 dec_two <- function(x) {
  return (format(round(x, 2), nsmall = 2));
}


