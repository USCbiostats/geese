#include <Rcpp.h>
#include "barry/barry.hpp"
#include "barry/models/geese.hpp"
using namespace Rcpp;

// Useful macros
#include "geese-utils.h"

//' @title Evolutionary terms
//' @description Model terms for both [geese] and [flock] objects.
//' @export
//' @param p An object of class [geese] or [flock].
//' @param funs Vector of function indices (starting from zero).
//' @param duplication Integer. 0 for speciation, 1 for duplication and 2 for either.
//' @name geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int rule_limit_changes(
    SEXP p,
    int term_pos,
    int lb,
    int ub,
    size_t duplication = 1
) {

  IF_GEESE(p) {

    Rcpp::XPtr<geese::Geese>ptr(p);
    geese::rule_dyn_limit_changes(ptr->get_support_fun(), term_pos, lb, ub, duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr<geese::Flock>ptr(p);
    geese::rule_dyn_limit_changes(ptr->get_support_fun(), term_pos, lb, ub, duplication);

  } IF_NEITHER()

  return 0;

}
