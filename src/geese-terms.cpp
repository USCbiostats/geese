#include <random>
#include <iostream>
#include <string>
#include <algorithm>
#include "barry/barry.hpp"
#include "barry/models/geese.hpp"
#include <Rcpp.h>
using namespace Rcpp;

// Useful macros
#ifndef GEESE_FLOCK_CASES
#define GEESE_FLOCK_CASES 1
#define IF_GEESE(a) if (Rf_inherits((a), "geese"))
#define IF_FLOCK(a) else if (Rf_inherits((a), "flock"))
#define IF_NEITHER() else stop("The passed object is neither a 'geese' nor a 'flock'.");
#endif

//' @title Evolutionary terms
//' @export
//' @name geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int term_gains(
    SEXP p,
    std::vector<unsigned int> & funs,
    bool duplication = true
) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylocounters::counter_gains(ptr->counters, funs, duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylocounters::counter_gains(&ptr->support.counters, funs, duplication);

  } IF_NEITHER()

  return 0;

}

//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int term_loss(SEXP p, std::vector<unsigned int> & funs,
              bool duplication = true) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylocounters::counter_loss(ptr->counters, funs, duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylocounters::counter_loss(&ptr->support.counters, funs, duplication);

  } IF_NEITHER()

  return 0;

}

//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int term_cogain(SEXP p, unsigned int a, unsigned int b) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylocounters::counter_cogain(ptr->counters, a, b);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylocounters::counter_cogain(&ptr->support.counters, a, b);

  } IF_NEITHER()

  return 0;

}

//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int term_neofun(SEXP p, unsigned int a, unsigned int b) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylocounters::counter_neofun(ptr->counters, a, b);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylocounters::counter_neofun(&ptr->support.counters, a, b);

  } IF_NEITHER()

  return 0;

}

//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int term_subfun(SEXP p, unsigned int a, unsigned int b) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylocounters::counter_subfun(ptr->counters, a, b);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylocounters::counter_subfun(&ptr->support.counters, a, b);

  } IF_NEITHER()


  return 0;

}

//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int term_maxfuns(
    SEXP p, unsigned int lb, unsigned int ub,
    bool duplication = true) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylocounters::counter_maxfuns(
      ptr->counters, lb, ub, duplication
    );

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylocounters::counter_maxfuns(
      &ptr->support.counters, lb, ub, duplication
    );

  } IF_NEITHER()

  return 0;

}


//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
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
