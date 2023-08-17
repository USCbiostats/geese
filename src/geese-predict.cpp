#include <Rcpp.h>
#include "barry/barry.hpp"
#include "barry/models/geese.hpp"
using namespace Rcpp;

#include "geese-utils.h"

//' @export
//' @rdname geese-common
// [[Rcpp::export(rng = false)]]
std::vector< std::vector< double > > predict_geese(
    SEXP p,
    const std::vector< double > & par,
    bool leave_one_out        = false,
    bool use_reduced_sequence = true,
    bool only_annotated       = false
  ) {

  // Preparing output
  std::vector< std::vector< double > > res(0u);

  IF_GEESE(p) {

    Rcpp::XPtr<geese::Geese>ptr(p);
    res = ptr->predict(
      par, nullptr, leave_one_out, only_annotated, use_reduced_sequence
      );

  } IF_FLOCK(p) {

    stop("Use -predict_flock- instead.");

  } IF_NEITHER()

  return res;

}

//' @export
//' @rdname geese-common
// [[Rcpp::export(rng = false)]]
std::vector< std::vector< std::vector< double > > > predict_flock(
    SEXP p,
    const std::vector< double > & par,
    bool leave_one_out        = false,
    bool use_reduced_sequence = true,
    bool only_annotated       = false
) {

  // Preparing output
  std::vector< std::vector< std::vector< double > > > res(0u);

  

  IF_GEESE(p) {

    stop("Use -predict_flock- instead.");

  } IF_FLOCK(p) {

    Rcpp::XPtr<geese::Flock>ptr(p);
    for (auto& d : ptr->dat)
      res.push_back(
        d.predict(
          par, nullptr, leave_one_out, only_annotated, use_reduced_sequence
        )
      );


  } IF_NEITHER()

    return res;

}


//' @export
//' @rdname geese-common
// [[Rcpp::export(rng = true)]]
std::vector< std::vector< double > > predict_geese_simulate(
    SEXP p,
    const std::vector< double > & par,
    size_t nsim,
    bool use_reduced_sequence = true,
    int seed = -1
    ) {

  // Preparing output
  std::vector< std::vector< double > > res(0u);

  if (seed < 0)
    seed = static_cast<size_t>(R::runif(0, 1e6));

  IF_GEESE(p) {

    Rcpp::XPtr<geese::Geese>ptr(p);

    ptr->set_seed(seed);
    res = ptr->predict_sim(
      par, use_reduced_sequence, nsim
    );

  } IF_FLOCK(p) {

    stop("Use -predict_flock_simulate- instead.");

  } IF_NEITHER()

  return res;

}
