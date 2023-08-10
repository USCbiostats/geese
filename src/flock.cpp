#include <Rcpp.h>
#include "barry/barry.hpp"
#include "barry/models/geese.hpp"

using namespace Rcpp;

#include "geese-utils.h"

//' @title Flock of GEESE
//' @description Flocks are pooled models which include multiple phylums.
//' @name flock-class
//' @export
//' @aliases flock
// [[Rcpp::export(rng = false)]]
SEXP new_flock() {

  Rcpp::XPtr<geese::Flock> dat(new geese::Flock());

  dat.attr("class") = "flock";

  return dat;

}

//' @rdname flock-class
//' @param p An object of class `flock`.
//' @param annotations Vector of integer vectors with annotations.
//' @param geneid integer vector with gene ids.
//' @param parent integer vector with parent gene id.
//' @param duplication logical vector indicating the type of event.
//' @export
// [[Rcpp::export(rng = false, invisible = true)]]
int add_geese(
    SEXP p,
    std::vector< std::vector< size_t > > & annotations,
    std::vector< size_t > & geneid,
    std::vector< int > & parent,
    std::vector< bool > & duplication
) {

  CHECK_FLOCK(p)

  Rcpp::XPtr<geese::Flock>ptr(p);
  ptr->add_data(annotations, geneid, parent, duplication);

  return 0;

}

