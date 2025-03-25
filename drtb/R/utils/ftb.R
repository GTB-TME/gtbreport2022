


#' Number formatter according to GTB rounding rules
#'
#' Formats vectors of numbers (<2 billion)
#'
#' GTB rounding convention:
#'
#' - 0 is written as "0" (output "0")
#'
#' - values under 0.1 are written "<0.1" ("<0.1")
#'
#' - from 0.1 to under 1 are rounded to 2
#'
#'   significant figures (0.NN)
#'
#' - from 1 to under 10 are rounded to 2 significant figures ("N.N")
#'
#' - 10 to under 100 are rounded to 2 significant figures ("NN")
#'
#' - 100 to under 1000 are rounded to 3 significant figures ("NNN")
#'
#' - 1000 upwards are rounded to 3 significant figures ("N NN0 000")
#'
#' - data that are not reported, but could be are represented
#'   as empty cells and should be accompanied by a footnote.
#'
#' - data that cannot be calculated, either because of missing
#'   data, data was not requested, or any other reason are represented
#'   with a dash.
#'
#' When the number represents thousands, show numbers between 0.01 and 0.1:
#'
#' - values under 0.01 are written "<0.01"
#'
#' - values between 0.01 and under 0.1 are rounded to
#'   2 significant figures ("0.0NN")
#'
#' @param x Vector of numbers
#' @examples
#' ftb(348838)
#' ftb(c(0.0359, 0.00036))
#'
#' @export
#'
ftb <- Vectorize(function(x) {
  # formatter according to GTB rounding rules
  # https://docs.google.com/document/d/1cu_syknBiF3scAX9d7hEfN0LZwdG40i8ttN6yua2xTQ/edit
  #' @param x vector of values
  #' @export
  stopifnot(!is.character(x))
  stopifnot(x < 2e9)

  if (!is.na(x)) {
    smallpos <- x > 0 & x < 0.01
    one2ten <- x >= 1 & x < 10
    zero2one <- x >= 0.1 & x < 1

    dg <- ifelse(abs(x) > 0.01 & abs(x) < 100, 2, 3)
    x2 <- signif(x, dg)

    trailing.0 <- x2 == round2(x) & one2ten == TRUE
    trailing0 <- x2 * 10 == round2(x * 10) & zero2one == TRUE & x2 < 1

    x2 <-
      format(
        x2,
        digits = dg,
        nsmall = 0L,
        big.mark = " ",
        justify = 'right',
        drop0trailing = T,
        scientific = F
      )
    if (smallpos)
      x2 <- '<0.01'
    if (trailing.0)
      x2 <- paste0(x2, '.0')
    if (trailing0)
      x2 <- paste0(x2, '0')
  } else
    x2 <- '-'
  return(x2)
}, 'x')



#' always round 0.5 up
#'
round2 <- function(x, digits = 0) sign(x) * trunc(abs(x) * 10^digits + 0.5) / 10^digits
