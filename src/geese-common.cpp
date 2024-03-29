#include <Rcpp.h>
#include "barry/barry.hpp"
#include "barry/models/geese.hpp"
using namespace Rcpp;

#include "geese-utils.h"

//' @title Common functions for `geese` and `flock`.
//' @name geese-common
//' @param p An object of class [geese] or [flock].
//' @param verb Integer scalar. When `>1`, it will print a progress bar during
//' the initialization of the process of length `verb`.
//' @export
//' @aliases flock-common
//' @details `init_model` initializes the model. This triggers the calculation
//' of the support using the vector of terms included. Initializing a model can
//' only be done once.
// [[Rcpp::export(rng = false, invisible = true)]]
int init_model(SEXP p, int verb = 80) {

  IF_GEESE(p) {

    Rcpp::XPtr<geese::Geese>ptr(p);
    ptr->init(verb);

  } IF_FLOCK(p) {

    Rcpp::XPtr<geese::Flock>ptr(p);
    ptr->init(verb);

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

    Rcpp::XPtr<geese::Geese>ptr(p);
    return static_cast<int>(ptr->nterms());

  } IF_FLOCK(p) {

    Rcpp::XPtr<geese::Flock>ptr(p);
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

    Rcpp::XPtr<geese::Geese>ptr(p);
    return wrap(ptr->nnodes());

  } IF_FLOCK(p) {

    Rcpp::XPtr<geese::Flock>ptr(p);
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

    Rcpp::XPtr<geese::Flock>ptr(p);
    return ptr->ntrees();

  } IF_NEITHER()


}

//' @rdname geese-common
//' @export
// [[Rcpp::export(rng = false)]]
IntegerVector nleafs(SEXP p) {

  IF_GEESE(p) {

    Rcpp::XPtr<geese::Geese>ptr(p);
    return wrap(ptr->nleafs());

  } IF_FLOCK(p) {

    Rcpp::XPtr<geese::Flock>ptr(p);
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
    bool trunc_seq = true,
    int  ncores    = 1
  ) {

  IF_GEESE(p) {

    Rcpp::XPtr<geese::Geese>ptr(p);
    return ptr->likelihood(par, as_log, trunc_seq, static_cast<size_t>(ncores));

  } IF_FLOCK(p) {

    Rcpp::XPtr<geese::Flock>ptr(p);
    return ptr->likelihood_joint(par, as_log, trunc_seq, static_cast<size_t>(ncores));

  } IF_NEITHER()


}

//' @rdname geese-common
//' @export
// [[Rcpp::export(rng = false)]]
NumericMatrix get_probabilities(SEXP p) {

  CHECK_GEESE(p)

  Rcpp::XPtr<geese::Geese>ptr(p);
  size_t N = ptr->nodes.size();
  size_t M = ptr->get_states().size();
  NumericMatrix m(N, M);
  std::fill(m.begin(), m.end(), 0.0);

  for (auto& i : ptr->sequence) {
    size_t k = 0u;
    for (auto& p : ptr->nodes.at(i).subtree_prob)
      m(ptr->nodes.at(i).id, k++) = p;
  }

  return m;
}

//' @rdname geese-common
//' @export
// [[Rcpp::export(rng = false)]]
std::vector< IntegerVector > get_sequence(
    SEXP p,
    bool reduced_sequence = true
) {

  IF_GEESE(p) {

    Rcpp::XPtr<geese::Geese>ptr(p);
    auto res = reduced_sequence ? ptr->reduced_sequence : ptr->sequence;

    IntegerVector resR(res.size());

    for (auto i = 0u; i < res.size(); ++i)
      resR[i] = res[i];

    return std::vector< IntegerVector >({resR});

  } IF_FLOCK(p) {

    Rcpp::XPtr<geese::Flock>ptr(p);
    std::vector< IntegerVector > res;
    res.reserve(ptr->ntrees());

    size_t i = 0u;
    for (const auto & tree : ptr->dat)
    {
      auto tmp = reduced_sequence ? tree.reduced_sequence: tree.sequence;
      IntegerVector resR(tmp.size());
      for (auto j = 0u; j < tmp.size(); ++j)
        resR[i] = tmp[i];

      res.push_back(resR);

    }

    return res;

  } IF_NEITHER()

  return std::vector< IntegerVector >();

}

//' @rdname geese-common
//' @export
// [[Rcpp::export(rng = false, invisible = true)]]
int set_seed(SEXP p, size_t s) {

  IF_GEESE(p) {

    Rcpp::XPtr<geese::Geese>ptr(p);
    ptr->set_seed(s);

  } IF_FLOCK(p) {

    Rcpp::XPtr<geese::Flock>ptr(p);
    ptr->set_seed(s);

  } IF_NEITHER()

  return 0;
}

//' @rdname geese-common
//' @export
// [[Rcpp::export(rng = false)]]
std::vector< std::vector< size_t > > sim_geese(
    SEXP p,
    const std::vector<double> & par,
    int seed = -1
  ) {

  CHECK_GEESE(p)

  Rcpp::XPtr<geese::Geese> ptr(p);

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

  Rcpp::XPtr<geese::Geese> ptr(p);

  return ptr->observed_counts();

}

//' @rdname geese-common
//' @export
// [[Rcpp::export(rng = false, invisible = true)]]
int print_observed_counts(
    SEXP p
) {

  CHECK_GEESE(p)

  Rcpp::XPtr<geese::Geese> ptr(p);

  ptr->print_observed_counts();
  return 0;

}

//' @rdname geese-common
//' @export
// [[Rcpp::export(rng = false)]]
int support_size(SEXP p) {

  int ans;
  IF_GEESE(p) {

    Rcpp::XPtr<geese::Geese> ptr(p);
    ans = static_cast<int>(ptr->support_size());

  } IF_FLOCK(p) {

    Rcpp::XPtr<geese::Flock> ptr(p);
    ans = static_cast<int>(ptr->get_support_fun()->get_data().size());

  } IF_NEITHER()

  return ans;

}

//' @rdname geese-common
//' @export
// [[Rcpp::export(rng = false)]]
int parse_polytomies(SEXP p, bool verbose = true) {

  std::vector< size_t > vec_ans;
  IF_GEESE(p) {

    Rcpp::XPtr<geese::Geese> ptr(p);
//    (void) ptr->parse_polytomies(verbose, &vec_ans);
    return ptr->parse_polytomies(verbose);

  } IF_FLOCK(p) {

    Rcpp::XPtr<geese::Flock> ptr(p);
    //(void) ptr->parse_polytomies(verbose, &vec_ans);
    return ptr->parse_polytomies(verbose);

  } IF_NEITHER()

  return 0; //wrap(vec_ans);

}

//' @rdname geese-common
//' @export
// [[Rcpp::export(rng = false)]]
int nfuns(SEXP p) {

  int ans = 0;
  IF_GEESE(p) {

    Rcpp::XPtr<geese::Geese> ptr(p);
    ans = static_cast<int>(ptr->nfuns());

  } IF_FLOCK(p) {

    Rcpp::XPtr<geese::Flock> ptr(p);
    ans = static_cast<int>(ptr->nfuns());

  } IF_NEITHER()

  return ans;

}

//' @rdname geese-common
//' @export
// [[Rcpp::export(rng = false, name = "names.geese")]]
std::vector< std::string > names_geese(SEXP p) {

  IF_GEESE(p) {

    Rcpp::XPtr<geese::Geese> ptr(p);
    return ptr->colnames();

  } IF_FLOCK(p){

    Rcpp::XPtr<geese::Flock> ptr(p);
    return ptr->colnames();


  } IF_NEITHER()

}

//' Compute the transition probability
//' @param p An object of class `geese` or `flock`.
//' @param params A vector of model parameters.
//' @param duplication logical scalar. Type of evolutionary event.
//' @param state logical vector. State of the parent.
//' @param array matrix indicating the state of the offspring (rows = function, cols = offspring).
//' @param as_log logical scalar. When `TRUE` returns the log.
//' @export
// [[Rcpp::export(rng = false)]]
double transition_prob(
    SEXP p,
    const std::vector< double > & params,
    bool duplication,
    const std::vector< bool > & state,
    const IntegerMatrix array,
    bool as_log = false
  ) {

  if (state.size() != static_cast<size_t>(array.nrow()))
    stop("The length of -state- does not match the number of functions in -nrow-.");

  geese::PhyloArray A(array.nrow(), array.ncol());
  A.set_data(
    new geese::NodeData(std::vector<double>(1.0, array.ncol()), state, duplication),
    true
  );

  for (int i = 0; i < array.nrow(); ++i)
    for (int j = 0; j < array.ncol(); ++j)
      A(i, j) = array(i, j);

  IF_GEESE(p)
  {

    // Preparing data
    Rcpp::XPtr<geese::Geese> ptr(p);
    return ptr->get_model()->likelihood(params, A, -1, as_log);

  } IF_FLOCK(p) {

    // Preparing data
    Rcpp::XPtr<geese::Flock> ptr(p);
    return ptr->get_model()->likelihood(params, A, -1, as_log);

  } IF_NEITHER()

  return 0.0;

}

//' @export
//' @rdname transition_prob
//' @param i,j Locations (index from zero) of the cell to compute.
//' @details `conditional_prob` computes the so-called "Gibbs sampling"
//' probability, in which the likelihood of observing Y(i,j) = 1 conditional
//' on the rest of the data is computed.
// [[Rcpp::export(rng = false)]]
double conditional_prob(
    SEXP p,
    const std::vector< double > & params,
    bool duplication,
    const std::vector< bool > & state,
    const IntegerMatrix array,
    size_t i, size_t j,
    bool as_log = false
) {

  if (state.size() != static_cast<size_t>(array.nrow()))
    stop("The length of -state- does not match the number of functions in -nrow-.");

  geese::PhyloArray A(array.nrow(), array.ncol());
  A.set_data(
    new geese::NodeData(std::vector<double>(1.0, array.ncol()), state, duplication),
    true
  );

  for (int i = 0; i < array.nrow(); ++i)
    for (int j = 0; j < array.ncol(); ++j)
      A(i, j) = array(i, j);

  IF_GEESE(p)
  {

    // Preparing data
    Rcpp::XPtr<geese::Geese> ptr(p);
    return ptr->get_model()->conditional_prob(A, params, i, j);

  } IF_FLOCK(p) {

    // Preparing data
    Rcpp::XPtr<geese::Flock> ptr(p);
    return ptr->get_model()->conditional_prob(A, params, i, j);

  } IF_NEITHER()

    return 0.0;

}

// [[Rcpp::export(rng = false, invisible = true)]]
int print_geese(SEXP p)
{

  IF_GEESE(p)
  {

    // Preparing data
    Rcpp::XPtr<geese::Geese> ptr(p);
    ptr->print();

  } IF_FLOCK(p) {

    // Preparing data
    Rcpp::XPtr<geese::Flock> ptr(p);
    ptr->print();

  } IF_NEITHER()

  return 0;

}

// [[Rcpp::export(rng = false, invisible = true)]]
int print_nodes_cpp(SEXP p)
{
  
    IF_GEESE(p)
    {
  
      // Preparing data
      Rcpp::XPtr<geese::Geese> ptr(p);
      ptr->print_nodes();
  
    } IF_FLOCK(p) {
  
      // Preparing data
      printf_barry("Not implemented yet for flocks.\n");
  
    } IF_NEITHER()
  
    return 0;
  
}

//' Returns the support of the model
//' @export
// [[Rcpp::export(rng = false)]]
std::vector< NumericMatrix > get_support(SEXP p)
{

  std::vector< double > dat;
  std::vector< size_t > * sizes;
  std::vector< size_t > * sizes_acc;
  size_t n_terms;
  std::vector< std::string > vnames = {"weights"};

  IF_GEESE(p)
  {

    // Preparing data
    Rcpp::XPtr<geese::Geese> ptr(p);
    dat       = *(ptr->get_model()->get_stats_support());
    sizes     = ptr->get_model()->get_stats_support_sizes();
    sizes_acc = ptr->get_model()->get_stats_support_sizes_acc();
    n_terms   = ptr->nterms() - ptr->nfuns() + 1u;

    auto names = ptr->colnames();
    for (auto n : names)
      vnames.push_back(n);

  } IF_FLOCK(p) {

    // Preparing data
    Rcpp::XPtr<geese::Flock> ptr(p);
    dat       = *(ptr->get_model()->get_stats_support());
    sizes     = ptr->get_model()->get_stats_support_sizes();
    sizes_acc = ptr->get_model()->get_stats_support_sizes_acc();
    n_terms   = ptr->nterms() - ptr->nfuns() + 1u;

    auto names = ptr->colnames();
    for (auto n : names)
      vnames.push_back(n);

  } IF_NEITHER()


  std::vector< NumericMatrix > ans;
  for (size_t g = 0u; g < sizes_acc->size(); ++g)
  {
    size_t n_rows = sizes->operator[](g);
    NumericMatrix tmp(n_rows, n_terms);

    for (auto i = 0u; i < n_rows; ++i)
      for (auto j = 0u; j < n_terms; ++j)
        tmp(i, j) = dat[
          sizes_acc->operator[](g) * n_terms +
          i * n_terms + j
        ];

    tmp.attr("dimnames") = List::create(
      R_NilValue,
      wrap(vnames)
    );

    ans.push_back(clone(tmp));

  }

  return ans;


}

//' @export
// [[Rcpp::export(name = "get_annotations", rng = false)]]
IntegerMatrix get_annotations_cpp(SEXP p) {

  Rcpp::XPtr<geese::Geese> ptr(p);
  
  auto annotations = ptr->get_annotations();
  size_t nnodes = ptr->nnodes();
  size_t nfuns  = ptr->nfuns();

  IntegerMatrix ans(nnodes, nfuns);

  for (size_t i = 0u; i < nnodes; ++i)
  {
    for (size_t j = 0u; j < nfuns; ++j)
    {
      int tmp = static_cast<int>(annotations[j * nnodes + i]);
      ans(i, j) = tmp == 9 ? NA_INTEGER : tmp;

    }
  }

  return ans;

}