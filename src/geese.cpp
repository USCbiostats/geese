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
    ListOf< IntegerVector > & annotations,
    std::vector< size_t > & geneid,
    std::vector< int > & parent,
    std::vector< bool > & duplication
) {

  // Preprocessing annotations. (checking for data types)
  std::vector< std::vector< size_t > > annotations2;
  for (IntegerVector entry : annotations) {

    std::vector< size_t > tmp;
    for (size_t i = 0u; i < entry.size(); ++i) {
      if (Rcpp::traits::is_na<INTSXP>(entry[i]) || entry[i] == 9)
        tmp.push_back(9);
      else if (entry[i] != 0 && entry[i] != 1)
        stop("Invalid annotation value: %d", entry[i]);
      else
        tmp.push_back(entry[i]);

    }

    // Rcpp::print(wrap(tmp));

    annotations2.push_back(tmp);
  }

  // Rcpp::print(wrap(annotations2));

  Rcpp::XPtr<geese::Geese> dat(
      new geese::Geese(annotations2, geneid, parent, duplication
      ));

  dat.attr("class") = "geese";

  return dat;

}
