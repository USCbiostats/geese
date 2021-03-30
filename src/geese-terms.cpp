#include <random>
#include <iostream>
#include <string>
#include <algorithm>
#include "barry/barry.hpp"
#include "barry/models/geese.hpp"
#include <Rcpp.h>
using namespace Rcpp;

/******************************************************************************/

//' @title Evolutionary terms
//' @export
//' @name geese-terms
// [[Rcpp::export(rng = false)]]
int term_gains(
    SEXP p,
    std::vector<unsigned int> & funs,
    bool duplication = true
) {

  Rcpp::XPtr< Geese >ptr(p);
  phylocounters::counter_gains(ptr->counters, funs, duplication);
  return 0;

}

//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false)]]
int term_loss(SEXP p, std::vector<unsigned int> & funs,
              bool duplication = true) {

  Rcpp::XPtr< Geese >ptr(p);
  phylocounters::counter_loss(ptr->counters, funs, duplication);
  return 0;

}

//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false)]]
int term_cogain(SEXP p, unsigned int a, unsigned int b) {

  Rcpp::XPtr< Geese >ptr(p);
  phylocounters::counter_cogain(ptr->counters, a, b);
  return 0;

}

//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false)]]
int term_neofun(SEXP p, unsigned int a, unsigned int b) {

  Rcpp::XPtr< Geese >ptr(p);
  phylocounters::counter_neofun(ptr->counters, a, b);
  return 0;

}

//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false)]]
int term_subfun(SEXP p, unsigned int a, unsigned int b) {

  Rcpp::XPtr< Geese >ptr(p);
  phylocounters::counter_subfun(ptr->counters, a, b);
  return 0;

}

//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false)]]
int term_maxfuns(
    SEXP p, unsigned int lb, unsigned int ub,
    bool duplication = true) {

  Rcpp::XPtr< Geese >ptr(p);
  phylocounters::counter_maxfuns(
    ptr->counters, lb, ub,
    duplication
  );

  return 0;

}


//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false)]]
int term_overall_changes(SEXP p, bool duplication = true) {

  Rcpp::XPtr< Geese >ptr(p);
  phylocounters::counter_overall_changes(ptr->counters, duplication);
  return 0;

}

//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false)]]
int term_kgains(
    SEXP p,
    std::vector<unsigned int> & funs,
    int k = 1,
    bool duplication = true) {

  Rcpp::XPtr< Geese >ptr(p);
  phylocounters::counter_gains_k_offspring(
    ptr->counters, funs,
    static_cast<unsigned int>(k), duplication);
  return 0;

}

//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false)]]
int term_neofun_a2b(
    SEXP p,
    int a,
    int b,
    bool duplication = true
) {

  Rcpp::XPtr< Geese >ptr(p);
  phylocounters::counter_neofun_a2b(
    ptr->counters, a, b, duplication);
  return 0;

}



//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false)]]
int term_genes_changing(
    SEXP p,
    bool duplication = true
) {

  Rcpp::XPtr< Geese >ptr(p);
  phylocounters::counter_genes_changing(ptr->counters, duplication);
  return 0;

}
