#include "barry/barry.hpp"
#include "barry/models/geese.hpp"
#include <Rcpp.h>
using namespace Rcpp;

#include "geese-utils.hpp"

//' @title Common functions for `geese` and `flock`.
//' @name geese-common
//' @param p An object of class [geese] or [flock].
//' @export
//' @aliases flock-common
//' @details `init_model` initializes the model. This triggers the calculation
//' of the support using the vector of terms included. Initializing a model can
//' only be done once.
// [[Rcpp::export(rng = false, invisible = true)]]
int init_model(SEXP p) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    ptr->init();

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    ptr->init();

  } IF_NEITHER()

  return 0;

}

//' @rdname geese-common
//' @export
//' @returns `nterms` returns the number of terms included in the model. This
//' is different from the number of parameters as the later includes the number
//' of functions in the data.
// [[Rcpp::export(rng = false)]]
int nterms(SEXP p) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    return static_cast<int>(ptr->nterms());

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    return static_cast<int>(ptr->nterms());

  } IF_NEITHER()

}

//' @rdname geese-common
//' @returns `nnodes` returns the number of nodes in the data, this includes
//' internal nodes. If `p` is a [flock], then it will return a vector of length
//' `ntrees()`.
//' @export
// [[Rcpp::export(rng = false)]]
IntegerVector nnodes(SEXP p) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    return wrap(ptr->nnodes());

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    return wrap(ptr->nnodes());

  } IF_NEITHER()


}

//' @rdname geese-common
//' @returns `ntrees` returns the number of trees in the model. For a geese object
//' this will be equal to one.
//' @export
// [[Rcpp::export(rng = false)]]
int ntrees(SEXP p) {

  IF_GEESE(p) {

    return 1;

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    return ptr->ntrees();

  } IF_NEITHER()


}

//' @rdname geese-common
//' @export
// [[Rcpp::export(rng = false)]]
IntegerVector nleafs(SEXP p) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    return wrap(ptr->nleafs());

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    return wrap(ptr->nleafs());

  } IF_NEITHER()


}

//' @rdname geese-common
//' @export
//' @param par Numeric vector of length `nterms()`.
//' @param trunc_seq When `TRUE` uses the truncated pruning sequence (see details).
//' @details Using the truncated pruning sequence (`trunc_seq = TRUE`) involves
//' traversing the trees throught the induced subtree. This is relevant in the case
//' that not all the leafs are annotated.
// [[Rcpp::export(rng = false)]]
double likelihood(
    SEXP p, const std::vector< double > & par,
    bool as_log    = false,
    bool trunc_seq = true
  ) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    return ptr->likelihood(par, as_log, trunc_seq);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    return ptr->likelihood_joint(par, as_log, trunc_seq);

  } IF_NEITHER()


}

//' @rdname geese-common
//' @export
// [[Rcpp::export(rng = false)]]
NumericMatrix get_probabilities(SEXP p) {

  CHECK_GEESE(p)

  Rcpp::XPtr< Geese >ptr(p);
  unsigned int N = ptr->nodes.size();
  unsigned int M = ptr->get_states().size();
  NumericMatrix m(N, M);
  std::fill(m.begin(), m.end(), 0.0);

  for (auto& i : ptr->sequence) {
    unsigned int k = 0u;
    for (auto& p : ptr->nodes.at(i).subtree_prob)
      m(ptr->nodes.at(i).id, k++) = p;
  }

  return m;
}

//' @rdname geese-common
//' @export
// [[Rcpp::export(rng = false)]]
std::vector< unsigned int > get_sequence(
    SEXP p,
    bool reduced_sequence = true
) {

  CHECK_GEESE(p)

  Rcpp::XPtr< Geese >ptr(p);
  return reduced_sequence ? ptr->reduced_sequence : ptr->sequence;
}

//' @rdname geese-common
//' @export
// [[Rcpp::export(rng = false, invisible = true)]]
int set_seed(SEXP p, unsigned int s) {

  IF_GEESE(p) {

    Rcpp::XPtr< Geese >ptr(p);
    ptr->set_seed(s);

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock >ptr(p);
    ptr->set_seed(s);

  } IF_NEITHER()

  return 0;
}

//' @rdname geese-common
//' @export
// [[Rcpp::export(rng = false)]]
std::vector< std::vector< unsigned int > > sim_geese(
    SEXP p,
    const std::vector<double> & par,
    int seed = -1
  ) {

  CHECK_GEESE(p)

  Rcpp::XPtr< Geese > ptr(p);

  if (seed > 0)
    ptr->set_seed(seed);

  return ptr->simulate(par);

}

//' @rdname geese-common
//' @export
// [[Rcpp::export(rng = false)]]
std::vector< std::vector< double > > observed_counts(
    SEXP p
) {

  CHECK_GEESE(p)

  Rcpp::XPtr< Geese > ptr(p);

  return ptr->observed_counts();

}

//' @rdname geese-common
//' @export
// [[Rcpp::export(rng = false, invisible = true)]]
int print_observed_counts(
    SEXP p
) {

  CHECK_GEESE(p)

  Rcpp::XPtr< Geese > ptr(p);

  ptr->print_observed_counts();
  return 0;

}

//' @rdname geese-common
//' @export
// [[Rcpp::export(rng = false)]]
int support_size(SEXP p) {

  int ans;
  IF_GEESE(p) {

    Rcpp::XPtr< Geese > ptr(p);
    ans = static_cast<int>(ptr->get_support()->get_data().size());

  } IF_FLOCK(p) {

    Rcpp::XPtr< Flock > ptr(p);
    ans = static_cast<int>(ptr->get_support()->get_data().size());

  } IF_NEITHER()

  return ans;

}
