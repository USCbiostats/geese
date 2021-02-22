#include <Rcpp.h>
#include "barry/barry.hpp"
#include "leaf.hpp"
using namespace Rcpp;

//' @title Aphylo2 model
//' @name aphylo2-class
//' @param annotations Vector of integer vectors with annotations.
//' @param geneid integer vector with gene ids
//' @param parent integer vector with parent gene id
//' @export
// [[Rcpp::export]]
SEXP new_model(
    std::vector< std::vector< unsigned int > > & annotations,
    std::vector< unsigned int > & geneid,
    std::vector< unsigned int > & parent
) {

  Rcpp::XPtr< APhyloModel > dat(
      new APhyloModel(annotations, geneid, parent
      ));

  return dat;

}

//' @rdname aphylo2-class
//' @export
// [[Rcpp::export]]
int init(SEXP p) {

  Rcpp::XPtr< APhyloModel >ptr(p);
  ptr->init();
  return 0;

}

//' @rdname aphylo2-class
//' @export
// [[Rcpp::export]]
double likelihood(SEXP p, const std::vector< double > & par) {

  Rcpp::XPtr< APhyloModel >ptr(p);
  return ptr->likelihood(par);

}

/******************************************************************************/

//' @title Evolutionary terms
//' @export
//' @name aphylo2-terms
// [[Rcpp::export]]
int term_gains(SEXP p, std::vector<unsigned int> & funs) {

  Rcpp::XPtr< APhyloModel >ptr(p);
  counter_gains(&ptr->counters, funs);
  return 0;

}

//' @export
//' @rdname aphylo2-terms
// [[Rcpp::export]]
int term_loss(SEXP p, std::vector<unsigned int> & funs) {

  Rcpp::XPtr< APhyloModel >ptr(p);
  counter_loss(&ptr->counters, funs);
  return 0;

}

//' @export
//' @rdname aphylo2-terms
// [[Rcpp::export]]
int term_cogain(SEXP p, unsigned int a, unsigned int b) {

  Rcpp::XPtr< APhyloModel >ptr(p);
  counter_cogain(&ptr->counters, a, b);
  return 0;

}

//' @export
//' @rdname aphylo2-terms
// [[Rcpp::export]]
int term_neofun(SEXP p, unsigned int a, unsigned int b) {

  Rcpp::XPtr< APhyloModel >ptr(p);
  phylocounters::counter_neofun(&ptr->counters, a, b);
  return 0;

}

//' @export
//' @rdname aphylo2-terms
// [[Rcpp::export]]
int term_subfun(SEXP p, unsigned int a, unsigned int b) {

  Rcpp::XPtr< APhyloModel >ptr(p);
  phylocounters::counter_subfun(&ptr->counters, a, b);
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
