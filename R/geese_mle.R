#' MLE for geese
#' @param amodel an object of class [geese-class]
#' @param initial Vector of initial parameters
#' @param ... Passed to [stats:optim]
#' @export
#' @importFrom stats optim
geese_mle <- function(
  amodel,
  initial = rep(0, nterms(amodel)),
  control = list(maxit = 1e3L, fnscale = -1),
  ...
  ) {


  # Normalized Log-likelihood function
  fun <- function(p) {

    ans <- likelihood(p = amodel, par = p, as_log = TRUE)

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
