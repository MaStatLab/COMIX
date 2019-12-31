#include "RcppArmadillo.h"
#include "perturbedSN_pmc.h"
#include "perturbedSN_helpers.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List perturbedSNcpp( arma::mat Y,
                           arma::uvec C,
                           Rcpp::List prior,
                           Rcpp::List pmc,
                           Rcpp::List state,
                           Rcpp::List initParticles, bool init )
{
  Rcpp::RNGScope scope;  
  PMC H(    Y,
            C,
            prior,
            pmc,
            state,
            initParticles, init
  );
  
  List chain = H.get_chain();
  
  List data = Rcpp::List::create(  
    Rcpp::Named( "Y" ) = Y,
    Rcpp::Named( "C" ) = C
  ) ; 
  
  return Rcpp::List::create(  
    Rcpp::Named( "chain" ) = chain,
    Rcpp::Named( "data" ) = data,
    Rcpp::Named( "prior" ) = prior,
    Rcpp::Named( "pmc" ) = pmc //,
    // Rcpp::Named( "control" ) = control
  ) ;    
}