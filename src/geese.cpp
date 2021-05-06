#include <Rcpp.h>
#include "barry/barry.hpp"
#include "barry/models/geese.hpp"
using namespace Rcpp;

#include "geese-utils.h"

//' @title GEne Evolutionary model using SufficiEncy (GEESE)
//' @name geese-class
//' @param p An object of class `geese`.
//' @param annotations Vector of integer vectors with annotations.
//' @param geneid integer vector with gene ids.
//' @param parent integer vector with parent gene id.
//' @param duplication logical vector indicating the type of events.
//' @export
//' @aliases geese
// [[Rcpp::export(rng = false)]]
SEXP new_geese(
    std::vector< std::vector< unsigned int > > & annotations,
    std::vector< unsigned int > & geneid,
    std::vector< int > & parent,
    std::vector< bool > & duplication
) {

  Rcpp::XPtr< Geese > dat(
      new Geese(annotations, geneid, parent, duplication
      ));

  dat.attr("class") = "geese";

  return dat;

}
