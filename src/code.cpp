#include <random>
#include <iostream>
#include <string>
#include <algorithm>
#include "barry/barry.hpp"
#include "barry/models/aphylomodel.hpp"
#include <Rcpp.h>
using namespace Rcpp;

//' @title Aphylo2 model
//' @name aphylo2-class
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

  Rcpp::XPtr< APhyloModel > dat(
      new APhyloModel(annotations, geneid, parent, duplication
      ));

  return dat;

}

//' @rdname aphylo2-class
//' @export
// [[Rcpp::export(rng = false)]]
int init(SEXP p) {

  Rcpp::XPtr< APhyloModel >ptr(p);
  ptr->init();
  return 0;

}

//' @rdname aphylo2-class
//' @export
// [[Rcpp::export(rng = false)]]
int nterms(SEXP p) {

  Rcpp::XPtr< APhyloModel >ptr(p);
  return ptr->nterms();

}

//' @rdname aphylo2-class
//' @export
// [[Rcpp::export(rng = false)]]
int nnodes(SEXP p) {

  Rcpp::XPtr< APhyloModel >ptr(p);
  return ptr->nnodes();

}

//' @rdname aphylo2-class
//' @export
// [[Rcpp::export(rng = false)]]
int nleafs(SEXP p) {

  Rcpp::XPtr< APhyloModel >ptr(p);
  return ptr->nleafs();

}

//' @rdname aphylo2-class
//' @export
// [[Rcpp::export(rng = false)]]
double likelihood(SEXP p, const std::vector< double > & par) {

  Rcpp::XPtr< APhyloModel >ptr(p);
  return ptr->likelihood(par);

}

//' @rdname aphylo2-class
//' @export
// [[Rcpp::export(rng = false)]]
NumericMatrix get_probabilities(SEXP p) {
  Rcpp::XPtr< APhyloModel >ptr(p);
  unsigned int N = ptr->nodes.size();
  unsigned int M = ptr->states.size();
  NumericMatrix m(N, M);
  std::fill(m.begin(), m.end(), 0.0);

  for (auto& i : ptr->sequence) {
    unsigned int k = 0u;
    for (auto& p : ptr->nodes.at(i).probabilities)
      m(ptr->nodes.at(i).id, k++) = p;
  }

  return m;
}

//' @rdname aphylo2-class
//' @export
// [[Rcpp::export(rng = false)]]
std::vector< unsigned int > get_sequence(SEXP p) {
  Rcpp::XPtr< APhyloModel >ptr(p);
  return ptr->sequence;
}

//' @rdname aphylo2-class
//' @export
// [[Rcpp::export(rng = false)]]
int set_seed(SEXP p, unsigned int s) {
  Rcpp::XPtr< APhyloModel > ptr(p);
  ptr->set_seed(s);
  return 0;
}

//' @rdname aphylo2-class
//' @export
// [[Rcpp::export(rng = false)]]
std::vector< std::vector< unsigned int > > sim_aphylo2(
    SEXP p,
    const std::vector<double> & par,
    int seed = -1
  ) {

  Rcpp::XPtr< APhyloModel > ptr(p);

  if (seed > 0)
    ptr->set_seed(seed);

  return ptr->simulate(par);

}

//' @rdname aphylo2-class
//' @export
// [[Rcpp::export(rng = false)]]
std::vector< std::vector< double > > observed_counts(
    SEXP p
) {

  Rcpp::XPtr< APhyloModel > ptr(p);

  return ptr->observed_counts();

}

//' @rdname aphylo2-class
//' @export
// [[Rcpp::export(rng = false)]]
int print_observed_counts(
    SEXP p
) {

  Rcpp::XPtr< APhyloModel > ptr(p);

  ptr->print_observed_counts();
  return 0;

}




/******************************************************************************/

//' @title Evolutionary terms
//' @export
//' @name aphylo2-terms
// [[Rcpp::export(rng = false)]]
int term_gains(
    SEXP p,
    std::vector<unsigned int> & funs,
    bool duplication = true
  ) {

  Rcpp::XPtr< APhyloModel >ptr(p);
  phylocounters::counter_gains(&ptr->counters, funs, duplication);
  return 0;

}

//' @export
//' @rdname aphylo2-terms
// [[Rcpp::export(rng = false)]]
int term_loss(SEXP p, std::vector<unsigned int> & funs,
              bool duplication = true) {

  Rcpp::XPtr< APhyloModel >ptr(p);
  phylocounters::counter_loss(&ptr->counters, funs, duplication);
  return 0;

}

//' @export
//' @rdname aphylo2-terms
// [[Rcpp::export(rng = false)]]
int term_cogain(SEXP p, unsigned int a, unsigned int b) {

  Rcpp::XPtr< APhyloModel >ptr(p);
  phylocounters::counter_cogain(&ptr->counters, a, b);
  return 0;

}

//' @export
//' @rdname aphylo2-terms
// [[Rcpp::export(rng = false)]]
int term_neofun(SEXP p, unsigned int a, unsigned int b) {

  Rcpp::XPtr< APhyloModel >ptr(p);
  phylocounters::counter_neofun(&ptr->counters, a, b);
  return 0;

}

//' @export
//' @rdname aphylo2-terms
// [[Rcpp::export(rng = false)]]
int term_subfun(SEXP p, unsigned int a, unsigned int b) {

  Rcpp::XPtr< APhyloModel >ptr(p);
  phylocounters::counter_subfun(&ptr->counters, a, b);
  return 0;

}

//' @export
//' @rdname aphylo2-terms
// [[Rcpp::export(rng = false)]]
int term_maxfuns(SEXP p, unsigned int lb, unsigned int ub) {

  Rcpp::XPtr< APhyloModel >ptr(p);
  phylocounters::counter_maxfuns(&ptr->counters, lb, ub);
  return 0;

}


//
// std::vector< std::vector< uint > > annotations = {
//   {0, 0, 1, 9, 9, 9},
//   {1, 1, 9, 0, 9, 9},
//   {1, 9, 0, 1, 9, 9}
// };
//
// std::vector< uint > geneid = {0, 1, 2, 3, 4, 5};
// std::vector< uint > parent = {4, 4, 5, 5, 6, 6};
//
//
// // Specifying the terms
// PhyloCounters counters;
// counter_gains(&counters, {0, 1, 2});
// counter_loss(&counters, {0, 1, 2});
// counter_cogain(&counters, 0, 1);
// counter_cogain(&counters, 0, 2);
// counter_cogain(&counters, 1, 2);
//
// APhyloModel dat(
//     annotations, geneid, parent, counters
// );
//
// // Model parameters
// std::vector< double > par = {
//   // Main parameters
//   .1, .1, .1, .1, .1, .1, .1, .1, .1,
//   // Root probabilities
//   .1, .1, .1
// };
//
// // Stargint to measure time
// // auto start = std::chrono::system_clock::now();
// double ans = dat.likelihood(par);
