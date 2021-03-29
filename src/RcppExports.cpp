// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// term_gains
int term_gains(SEXP p, std::vector<unsigned int>& funs, bool duplication);
RcppExport SEXP _geese_term_gains(SEXP pSEXP, SEXP funsSEXP, SEXP duplicationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int>& >::type funs(funsSEXP);
    Rcpp::traits::input_parameter< bool >::type duplication(duplicationSEXP);
    rcpp_result_gen = Rcpp::wrap(term_gains(p, funs, duplication));
    return rcpp_result_gen;
END_RCPP
}
// term_loss
int term_loss(SEXP p, std::vector<unsigned int>& funs, bool duplication);
RcppExport SEXP _geese_term_loss(SEXP pSEXP, SEXP funsSEXP, SEXP duplicationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int>& >::type funs(funsSEXP);
    Rcpp::traits::input_parameter< bool >::type duplication(duplicationSEXP);
    rcpp_result_gen = Rcpp::wrap(term_loss(p, funs, duplication));
    return rcpp_result_gen;
END_RCPP
}
// term_cogain
int term_cogain(SEXP p, unsigned int a, unsigned int b);
RcppExport SEXP _geese_term_cogain(SEXP pSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type a(aSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(term_cogain(p, a, b));
    return rcpp_result_gen;
END_RCPP
}
// term_neofun
int term_neofun(SEXP p, unsigned int a, unsigned int b);
RcppExport SEXP _geese_term_neofun(SEXP pSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type a(aSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(term_neofun(p, a, b));
    return rcpp_result_gen;
END_RCPP
}
// term_subfun
int term_subfun(SEXP p, unsigned int a, unsigned int b);
RcppExport SEXP _geese_term_subfun(SEXP pSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type a(aSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(term_subfun(p, a, b));
    return rcpp_result_gen;
END_RCPP
}
// term_maxfuns
int term_maxfuns(SEXP p, unsigned int lb, unsigned int ub, bool duplication);
RcppExport SEXP _geese_term_maxfuns(SEXP pSEXP, SEXP lbSEXP, SEXP ubSEXP, SEXP duplicationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< bool >::type duplication(duplicationSEXP);
    rcpp_result_gen = Rcpp::wrap(term_maxfuns(p, lb, ub, duplication));
    return rcpp_result_gen;
END_RCPP
}
// term_overall_changes
int term_overall_changes(SEXP p, bool duplication);
RcppExport SEXP _geese_term_overall_changes(SEXP pSEXP, SEXP duplicationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    Rcpp::traits::input_parameter< bool >::type duplication(duplicationSEXP);
    rcpp_result_gen = Rcpp::wrap(term_overall_changes(p, duplication));
    return rcpp_result_gen;
END_RCPP
}
// term_kgains
int term_kgains(SEXP p, std::vector<unsigned int>& funs, int k, bool duplication);
RcppExport SEXP _geese_term_kgains(SEXP pSEXP, SEXP funsSEXP, SEXP kSEXP, SEXP duplicationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int>& >::type funs(funsSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< bool >::type duplication(duplicationSEXP);
    rcpp_result_gen = Rcpp::wrap(term_kgains(p, funs, k, duplication));
    return rcpp_result_gen;
END_RCPP
}
// term_neofun_a2b
int term_neofun_a2b(SEXP p, int a, int b, bool duplication);
RcppExport SEXP _geese_term_neofun_a2b(SEXP pSEXP, SEXP aSEXP, SEXP bSEXP, SEXP duplicationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    Rcpp::traits::input_parameter< bool >::type duplication(duplicationSEXP);
    rcpp_result_gen = Rcpp::wrap(term_neofun_a2b(p, a, b, duplication));
    return rcpp_result_gen;
END_RCPP
}
// term_genes_changing
int term_genes_changing(SEXP p, bool duplication);
RcppExport SEXP _geese_term_genes_changing(SEXP pSEXP, SEXP duplicationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    Rcpp::traits::input_parameter< bool >::type duplication(duplicationSEXP);
    rcpp_result_gen = Rcpp::wrap(term_genes_changing(p, duplication));
    return rcpp_result_gen;
END_RCPP
}
// new_model
SEXP new_model(std::vector< std::vector< unsigned int > >& annotations, std::vector< unsigned int >& geneid, std::vector< int >& parent, std::vector< bool >& duplication);
RcppExport SEXP _geese_new_model(SEXP annotationsSEXP, SEXP geneidSEXP, SEXP parentSEXP, SEXP duplicationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< std::vector< std::vector< unsigned int > >& >::type annotations(annotationsSEXP);
    Rcpp::traits::input_parameter< std::vector< unsigned int >& >::type geneid(geneidSEXP);
    Rcpp::traits::input_parameter< std::vector< int >& >::type parent(parentSEXP);
    Rcpp::traits::input_parameter< std::vector< bool >& >::type duplication(duplicationSEXP);
    rcpp_result_gen = Rcpp::wrap(new_model(annotations, geneid, parent, duplication));
    return rcpp_result_gen;
END_RCPP
}
// init
int init(SEXP p);
RcppExport SEXP _geese_init(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(init(p));
    return rcpp_result_gen;
END_RCPP
}
// nterms
int nterms(SEXP p);
RcppExport SEXP _geese_nterms(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(nterms(p));
    return rcpp_result_gen;
END_RCPP
}
// nnodes
int nnodes(SEXP p);
RcppExport SEXP _geese_nnodes(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(nnodes(p));
    return rcpp_result_gen;
END_RCPP
}
// nleafs
int nleafs(SEXP p);
RcppExport SEXP _geese_nleafs(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(nleafs(p));
    return rcpp_result_gen;
END_RCPP
}
// likelihood
double likelihood(SEXP p, const std::vector< double >& par);
RcppExport SEXP _geese_likelihood(SEXP pSEXP, SEXP parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    Rcpp::traits::input_parameter< const std::vector< double >& >::type par(parSEXP);
    rcpp_result_gen = Rcpp::wrap(likelihood(p, par));
    return rcpp_result_gen;
END_RCPP
}
// get_probabilities
NumericMatrix get_probabilities(SEXP p);
RcppExport SEXP _geese_get_probabilities(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(get_probabilities(p));
    return rcpp_result_gen;
END_RCPP
}
// get_sequence
std::vector< unsigned int > get_sequence(SEXP p);
RcppExport SEXP _geese_get_sequence(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(get_sequence(p));
    return rcpp_result_gen;
END_RCPP
}
// set_seed
int set_seed(SEXP p, unsigned int s);
RcppExport SEXP _geese_set_seed(SEXP pSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(set_seed(p, s));
    return rcpp_result_gen;
END_RCPP
}
// sim_geese
std::vector< std::vector< unsigned int > > sim_geese(SEXP p, const std::vector<double>& par, int seed);
RcppExport SEXP _geese_sim_geese(SEXP pSEXP, SEXP parSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type par(parSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_geese(p, par, seed));
    return rcpp_result_gen;
END_RCPP
}
// observed_counts
std::vector< std::vector< double > > observed_counts(SEXP p);
RcppExport SEXP _geese_observed_counts(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(observed_counts(p));
    return rcpp_result_gen;
END_RCPP
}
// print_observed_counts
int print_observed_counts(SEXP p);
RcppExport SEXP _geese_print_observed_counts(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(print_observed_counts(p));
    return rcpp_result_gen;
END_RCPP
}
// predictions
std::vector< std::vector< double > > predictions(SEXP p, const std::vector< double >& par, bool leave_one_out);
RcppExport SEXP _geese_predictions(SEXP pSEXP, SEXP parSEXP, SEXP leave_one_outSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    Rcpp::traits::input_parameter< const std::vector< double >& >::type par(parSEXP);
    Rcpp::traits::input_parameter< bool >::type leave_one_out(leave_one_outSEXP);
    rcpp_result_gen = Rcpp::wrap(predictions(p, par, leave_one_out));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_geese_term_gains", (DL_FUNC) &_geese_term_gains, 3},
    {"_geese_term_loss", (DL_FUNC) &_geese_term_loss, 3},
    {"_geese_term_cogain", (DL_FUNC) &_geese_term_cogain, 3},
    {"_geese_term_neofun", (DL_FUNC) &_geese_term_neofun, 3},
    {"_geese_term_subfun", (DL_FUNC) &_geese_term_subfun, 3},
    {"_geese_term_maxfuns", (DL_FUNC) &_geese_term_maxfuns, 4},
    {"_geese_term_overall_changes", (DL_FUNC) &_geese_term_overall_changes, 2},
    {"_geese_term_kgains", (DL_FUNC) &_geese_term_kgains, 4},
    {"_geese_term_neofun_a2b", (DL_FUNC) &_geese_term_neofun_a2b, 4},
    {"_geese_term_genes_changing", (DL_FUNC) &_geese_term_genes_changing, 2},
    {"_geese_new_model", (DL_FUNC) &_geese_new_model, 4},
    {"_geese_init", (DL_FUNC) &_geese_init, 1},
    {"_geese_nterms", (DL_FUNC) &_geese_nterms, 1},
    {"_geese_nnodes", (DL_FUNC) &_geese_nnodes, 1},
    {"_geese_nleafs", (DL_FUNC) &_geese_nleafs, 1},
    {"_geese_likelihood", (DL_FUNC) &_geese_likelihood, 2},
    {"_geese_get_probabilities", (DL_FUNC) &_geese_get_probabilities, 1},
    {"_geese_get_sequence", (DL_FUNC) &_geese_get_sequence, 1},
    {"_geese_set_seed", (DL_FUNC) &_geese_set_seed, 2},
    {"_geese_sim_geese", (DL_FUNC) &_geese_sim_geese, 3},
    {"_geese_observed_counts", (DL_FUNC) &_geese_observed_counts, 1},
    {"_geese_print_observed_counts", (DL_FUNC) &_geese_print_observed_counts, 1},
    {"_geese_predictions", (DL_FUNC) &_geese_predictions, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_geese(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
