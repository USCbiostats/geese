#' MCMC for aphylo2
#' @param amodel an object of class [aphylo2-class]
#' @param initial Vector of initial parameters
#' @param ... Passed to [fmcmc::MCMC]
#' @export
#' @importFrom fmcmc MCMC
aphylo2_mcmc <- function(
  amodel,
  initial = rep(0, nterms(amodel)),
  ...
  ) {


  # Normalized Log-likelihood function
  fun <- function(p) {

    ans <- -log(likelihood(amodel, p))

    if (!is.finite(ans))
      return(-.Machine$double.xmax * 1e-100)

    ans

  }

  # Calling the FMCMC
  fmcmc::MCMC(
    initial = initial,
    fun     = fun,
    ...
  )

}
