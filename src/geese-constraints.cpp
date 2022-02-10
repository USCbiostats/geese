#include <Rcpp.h>
#include "barry/barry.hpp"
#include "barry/models/geese.hpp"
using namespace Rcpp;

// Useful macros
#include "geese-utils.h"

using namespace barry::counters;

//' @title Evolutionary terms
//' @description Model terms for both [geese] and [flock] objects.
//' @export
//' @param p An object of class [geese] or [flock].
//' @param funs Vector of function indices (starting from zero).
//' @param duplication When `TRUE` indicates that this term is only valid for
//' duplication events.
//' @name geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int rule_limit_changes(
    SEXP p,
    int term_pos,
    int lb,
    int ub,
    bool duplication = true
) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylo::rule_dyn_limit_changes(ptr->get_support_fun(), term_pos, lb, ub, duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylo::rule_dyn_limit_changes(ptr->get_support_fun(), term_pos, lb, ub, duplication);

  } IF_NEITHER()

  return 0;

}
