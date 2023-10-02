#' MLE for geese
#' @param amodel an object of class [geese-class]
#' @param initial Vector of initial parameters
#' @param offset_list Integer vector with the indices of the parameters
#' that will set set to be fixed.
#' @param ... Passed to [stats::optim()]
#' @export
#' @importFrom stats optim
geese_mle <- function(
  amodel,
  initial = rep(0, nterms(amodel)),
  control = list(),
  offset_list = NULL,
  ncores = 1L,
  ...
  ) {

  if (length(control) == 0)
    control$maxit   <- 1e3L
    
  control$fnscale <- -1

  # Normalized Log-likelihood function
  if (length(offset_list)) {

    tmppar <- initial
    not_fixed <- setdiff(seq_along(initial), offset_list)
    fun <- function(p) {

      tmppar[not_fixed] <- p

      ans <- likelihood(p = amodel, par = tmppar, as_log = TRUE, ncores = ncores)

      if (!is.finite(ans))
        return(-.Machine$double.xmax * 1e-100)

      ans

    }

    do.call(
      stats::optim,
      c(
        list(
          par = initial[not_fixed],
          fn  = fun,
          control = control
        ),
        list(...)
      )
    )

  } else {

    fun <- function(p) {

      ans <- likelihood(p = amodel, par = p, as_log = TRUE, ncores = ncores)

      if (!is.finite(ans))
        return(-.Machine$double.xmax * 1e-100)

      ans

    }

    do.call(
      stats::optim,
      c(
        list(
          par = initial,
          fn  = fun,
          control = control
        ),
        list(...)
      )
    )

  }



}
