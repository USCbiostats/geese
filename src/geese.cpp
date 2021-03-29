#include <random>
#include <iostream>
#include <string>
#include <algorithm>
#include "barry/barry.hpp"
#include "barry/models/geese.hpp"
#include <Rcpp.h>
using namespace Rcpp;

//' @title GEne Evolutionary model using SufficiEncy (GEESE)
//' @name geese-class
//' @param annotations Vector of integer vectors with annotations.
//' @param geneid integer vector with gene ids
//' @param parent integer vector with parent gene id
//' @export
// [[Rcpp::export(rng = false)]]
SEXP new_model(
    std::vector< std::vector< unsigned int > > & annotations,
    std::vector< unsigned int > & geneid,
    std::vector< int > & parent,
    std::vector< bool > & duplication
) {

  Rcpp::XPtr< Geese > dat(
      new Geese(annotations, geneid, parent, duplication
      ));

  return dat;

}

//' @rdname geese-class
//' @export
// [[Rcpp::export(rng = false)]]
int init(SEXP p) {

  Rcpp::XPtr< Geese >ptr(p);
  ptr->init();
  return 0;

}

//' @rdname geese-class
//' @export
// [[Rcpp::export(rng = false)]]
int nterms(SEXP p) {

  Rcpp::XPtr< Geese >ptr(p);
  return ptr->nterms();

}

//' @rdname geese-class
//' @export
// [[Rcpp::export(rng = false)]]
int nnodes(SEXP p) {

  Rcpp::XPtr< Geese >ptr(p);
  return ptr->nnodes();

}

//' @rdname geese-class
//' @export
// [[Rcpp::export(rng = false)]]
int nleafs(SEXP p) {

  Rcpp::XPtr< Geese >ptr(p);
  return ptr->nleafs();

}

//' @rdname geese-class
//' @export
// [[Rcpp::export(rng = false)]]
double likelihood(SEXP p, const std::vector< double > & par) {

  Rcpp::XPtr< Geese >ptr(p);
  return ptr->likelihood(par);

}

//' @rdname geese-class
//' @export
// [[Rcpp::export(rng = false)]]
NumericMatrix get_probabilities(SEXP p) {
  Rcpp::XPtr< Geese >ptr(p);
  unsigned int N = ptr->nodes.size();
  unsigned int M = ptr->states.size();
  NumericMatrix m(N, M);
  std::fill(m.begin(), m.end(), 0.0);

  for (auto& i : ptr->sequence) {
    unsigned int k = 0u;
    for (auto& p : ptr->nodes.at(i).subtree_prob)
      m(ptr->nodes.at(i).id, k++) = p;
  }

  return m;
}

//' @rdname geese-class
//' @export
// [[Rcpp::export(rng = false)]]
std::vector< unsigned int > get_sequence(SEXP p) {
  Rcpp::XPtr< Geese >ptr(p);
  return ptr->sequence;
}

//' @rdname geese-class
//' @export
// [[Rcpp::export(rng = false)]]
int set_seed(SEXP p, unsigned int s) {
  Rcpp::XPtr< Geese > ptr(p);
  ptr->set_seed(s);
  return 0;
}

//' @rdname geese-class
//' @export
// [[Rcpp::export(rng = false)]]
std::vector< std::vector< unsigned int > > sim_geese(
    SEXP p,
    const std::vector<double> & par,
    int seed = -1
  ) {

  Rcpp::XPtr< Geese > ptr(p);

  if (seed > 0)
    ptr->set_seed(seed);

  return ptr->simulate(par);

}

//' @rdname geese-class
//' @export
// [[Rcpp::export(rng = false)]]
std::vector< std::vector< double > > observed_counts(
    SEXP p
) {

  Rcpp::XPtr< Geese > ptr(p);

  return ptr->observed_counts();

}

//' @rdname geese-class
//' @export
// [[Rcpp::export(rng = false)]]
int print_observed_counts(
    SEXP p
) {

  Rcpp::XPtr< Geese > ptr(p);

  ptr->print_observed_counts();
  return 0;

}

//' @export
//' @rdname geese-class
// [[Rcpp::export(rng = false)]]
std::vector< std::vector< double > > predictions(
    SEXP p,
    const std::vector< double > & par,
    bool leave_one_out = false
  ) {

  Rcpp::XPtr< Geese >ptr(p);

  // Baseline predictions
  std::vector< std::vector< double > > res = ptr->predict(par);

  if (!leave_one_out)
    return res;

  // In the case of leave one out, we need to update the predictions
  // accordingly
  std::vector< unsigned int > default_empty(ptr->nfuns(), 9u);
  for (auto& n : ptr->nodes) {

    if (n.second.is_leaf()) {

      Node & ntmp = n.second;

      // Recording the original annotation
      auto old_ann = ntmp.annotations;

      // Removing the entire gene
      ptr->update_annotations(ntmp.id, default_empty);

      // Making the prediction
      res[ntmp.id] = (ptr->predict(par))[ntmp.id];

      // Restoring the gene
      ptr->update_annotations(ntmp.id, old_ann);


    }

  }

  return res;

}
