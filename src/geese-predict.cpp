#include <random>
#include <iostream>
#include <string>
#include <algorithm>
#include "barry/barry.hpp"
#include "barry/models/geese.hpp"
#include <Rcpp.h>
using namespace Rcpp;

#include "geese-utils.hpp"

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

    Rcpp::XPtr< Geese >ptr(p);
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

    Rcpp::XPtr< Flock >ptr(p);
    for (auto& d : ptr->dat)
      res.push_back(
        d.predict(
          par, nullptr, leave_one_out, only_annotated, use_reduced_sequence
        )
      );


  } IF_NEITHER()

    return res;

}
