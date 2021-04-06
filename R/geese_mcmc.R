#' MCMC for geese
#' @param amodel an object of class [geese-class]
#' @param initial Vector of initial parameters
#' @param ... Passed to [fmcmc::MCMC]
#' @export
#' @importFrom fmcmc MCMC
geese_mcmc <- function(
  amodel,
  initial = rep(0, nterms(amodel)),
  prior   = function(p) {
    stats::dlogis(p, log = TRUE)
  },
  ...
  ) {


  # Normalized Log-likelihood function
  fun <- function(p) {

    ans <- likelihood(p = amodel, par = p, as_log = TRUE) + sum(prior(p)) * ntrees(amodel)

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