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
int term_gains(
    SEXP p,
    std::vector<unsigned int> & funs,
    unsigned int duplication = 1
) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylo::counter_gains(ptr->get_counters(), funs, duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylo::counter_gains(ptr->get_counters(), funs, duplication);

  } IF_NEITHER()

  return 0;

}

//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int term_loss(SEXP p, std::vector<unsigned int> & funs,
              unsigned int duplication = 1) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylo::counter_loss(ptr->get_counters(), funs, duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylo::counter_loss(ptr->get_counters(), funs, duplication);

  } IF_NEITHER()

  return 0;

}

//' @export
//' @rdname geese-terms
//' @param a,b Indices of functions (starting from zero)
// [[Rcpp::export(rng = false, invisible = true)]]
int term_cogain(SEXP p, unsigned int a, unsigned int b, unsigned int duplication = 1) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylo::counter_cogain(ptr->get_counters(), a, b, duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylo::counter_cogain(ptr->get_counters(), a, b, duplication);

  } IF_NEITHER()

  return 0;

}

//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int term_neofun(SEXP p, unsigned int a, unsigned int b, unsigned int duplication = 1) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylo::counter_neofun(ptr->get_counters(), a, b, duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylo::counter_neofun(ptr->get_counters(), a, b, duplication);

  } IF_NEITHER()

  return 0;

}

//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int term_subfun(
    SEXP p, unsigned int a, unsigned int b,
    unsigned int duplication = 1
  ) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylo::counter_subfun(ptr->get_counters(), a, b, duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylo::counter_subfun(ptr->get_counters(), a, b, duplication);

  } IF_NEITHER()


  return 0;

}

//' @export
//' @rdname geese-terms
//' @param lb,ub Integers, minimum and maximum number of changes.
// [[Rcpp::export(rng = false, invisible = true)]]
int term_maxfuns(
    SEXP p, unsigned int lb, unsigned int ub,
    unsigned int duplication = 1) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylo::counter_maxfuns(
      ptr->get_counters(), lb, ub, duplication
    );

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylo::counter_maxfuns(
      ptr->get_counters(), lb, ub, duplication
    );

  } IF_NEITHER()

  return 0;

}


//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int term_overall_changes(SEXP p, unsigned int duplication = 1) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylo::counter_overall_changes(ptr->get_counters(), duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylo::counter_overall_changes(ptr->get_counters(), duplication);

  } IF_NEITHER()


  return 0;

}

//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int term_overall_gains(SEXP p, unsigned int duplication = 1) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylo::counter_overall_gains(ptr->get_counters(), duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylo::counter_overall_gains(ptr->get_counters(), duplication);

  } IF_NEITHER()


  return 0;

}

//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int term_overall_loss(SEXP p, unsigned int duplication = 1) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylo::counter_overall_loss(ptr->get_counters(), duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylo::counter_overall_loss(ptr->get_counters(), duplication);

  } IF_NEITHER()


  return 0;

}


//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int term_kgains(
  SEXP p,
  std::vector<unsigned int> & funs,
  int k = 1,
  unsigned int duplication = 1
) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylo::counter_gains_k_offspring(
      ptr->get_counters(), funs,
      static_cast<unsigned int>(k), duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylo::counter_gains_k_offspring(
      ptr->get_counters(), funs,
      static_cast<unsigned int>(k), duplication);

  } IF_NEITHER()

  return 0;

}

//' @export
//' @rdname geese-terms
//' @details In the case of `term_neofun_a2b`, `a` represents the origin function
//' from which `b` is originated.
// [[Rcpp::export(rng = false, invisible = true)]]
int term_neofun_a2b(
    SEXP p,
    int a,
    int b,
    unsigned int duplication = 1
) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylo::counter_neofun_a2b(ptr->get_counters(), a, b, duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylo::counter_neofun_a2b(ptr->get_counters(), a, b, duplication);

  } IF_NEITHER()

  return 0;

}



//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int term_genes_changing(
    SEXP p,
    unsigned int duplication = 1
) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylo::counter_genes_changing(ptr->get_counters(), duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylo::counter_genes_changing(ptr->get_counters(), duplication);

  } IF_NEITHER()

  return 0;

}

//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int term_prop_genes_changing(
    SEXP p,
    unsigned int duplication = 1
) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylo::counter_prop_genes_changing(ptr->get_counters(), duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylo::counter_prop_genes_changing(ptr->get_counters(), duplication);

  } IF_NEITHER()

  return 0;

}



//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int term_coopt(
    SEXP p,
    unsigned int a,
    unsigned int b,
    unsigned int duplication = 1
) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylo::counter_co_opt(ptr->get_counters(), a, b, duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylo::counter_co_opt(ptr->get_counters(), a, b, duplication);

  } IF_NEITHER()

  return 0;

}

//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int term_k_genes_changing(
    SEXP p,
    unsigned int k,
    unsigned int duplication = 1
) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylo::counter_k_genes_changing(ptr->get_counters(), k, duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylo::counter_k_genes_changing(ptr->get_counters(), k, duplication);

  } IF_NEITHER()

  return 0;

}

//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int term_less_than_p_prop_genes_changing(
    SEXP p,
    double prop,
    unsigned int duplication = 1
) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylo::counter_less_than_p_prop_genes_changing(ptr->get_counters(), prop, duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylo::counter_less_than_p_prop_genes_changing(ptr->get_counters(), prop, duplication);

  } IF_NEITHER()

  return 0;

}

//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int term_pairwise_preserving(
    SEXP p,
    unsigned int funA,
    unsigned int funB,
    unsigned int duplication = 1
) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylo::counter_pairwise_preserving(ptr->get_counters(), funA, funB, duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylo::counter_pairwise_preserving(ptr->get_counters(), funA, funB, duplication);

  } IF_NEITHER()

  return 0;

}

//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int term_gains_from_0(
    SEXP p,
    std::vector<unsigned int> & fun,
    unsigned int duplication = 1
) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylo::counter_gains_from_0(ptr->get_counters(), fun, duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylo::counter_gains_from_0(ptr->get_counters(), fun, duplication);

  } IF_NEITHER()

  return 0;

}

//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int term_overall_gains_from_0(
    SEXP p,
    unsigned int duplication = 1
) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylo::counter_overall_gains_from_0(ptr->get_counters(), duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylo::counter_overall_gains_from_0(ptr->get_counters(), duplication);

  } IF_NEITHER()

  return 0;

}

//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int term_pairwise_first_gain(
    SEXP p,
    unsigned int funA,
    unsigned int funB,
    unsigned int duplication = 1
) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylo::counter_pairwise_first_gain(ptr->get_counters(), funA, funB, duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylo::counter_pairwise_first_gain(ptr->get_counters(), funA, funB, duplication);

  } IF_NEITHER()

  return 0;

}


//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int term_preserve_pseudogene(
    SEXP p,
    unsigned int funA,
    unsigned int funB,
    unsigned int duplication = 1
) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylo::counter_preserve_pseudogene(ptr->get_counters(), funA, funB, duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylo::counter_preserve_pseudogene(ptr->get_counters(), funA, funB, duplication);

  } IF_NEITHER()

  return 0;

}


//' @export
//' @rdname geese-terms
// [[Rcpp::export(rng = false, invisible = true)]]
int term_pairwise_overall_change(
    SEXP p,
    unsigned int duplication = 1
) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    phylo::counter_pairwise_overall_change(ptr->get_counters(), duplication);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    phylo::counter_pairwise_overall_change(ptr->get_counters(), duplication);

  } IF_NEITHER()

  return 0;

}



