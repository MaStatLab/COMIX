#ifndef GIBBS_H
#define GIBBS_H

#include "RcppArmadillo.h"
#include <omp.h>

using namespace Rcpp;
using namespace arma;

class PMC
{
private:
  mat Y;            // the data (n,J)
  uvec C;           // the group membership (n,1)
  int K;            // number of mixture components 
  int J;            // the number of groups
  int n;            // number of observations
  int p;            // observation dimension
  int num_particles, num_iter, num_burnin, num_thin, num_display;   // number of iteration, burnin and thinning
  // int seed;         // initial random seed
  
  /* --- hyperparameters --- */
  int length_chain;
  vec tau_a;
  double e0;
  vec b0;
  mat E0, invE0, B0, invB0;
  bool merge_step;
  double merge_par;
  double zeta;
  
  /* --- initial values --- */
  uvec T;

  /* --- storage --- */
  umat saveT;
  mat saveZ;
  cube saveXi, saveXi0, savePsi, saveAlpha, saveW;
  cube saveG, saveOmega, saveE;
  mat saveLog_py, savePerplexity;
  vec saveA0;
  umat saveNResampled;

  /* --- functions --- */
  void main_loop(Rcpp::List state, Rcpp::List prior, 
                 Rcpp::List initParticles, bool init);
  
  double sampleA0(double a0, arma::umat N, double a_par);
    
  arma::mat sampleLogWs(   arma::umat N,
                           double a0);
    
  Rcpp::List initialParticles( uvec T );
  
  Rcpp::List sampleXi( mat Y_k, uvec C_k, uvec N_k, Rcpp::List particles );
  
  Rcpp::List sampleG( mat Y_k, uvec C_k, uvec N_k, 
                     Rcpp::List particles, Rcpp::List prior);
    
  Rcpp::List samplePsi( mat Y_k, uvec C_k, uvec N_k, Rcpp::List particles,
                        Rcpp::List prior);
  
  Rcpp::List sampleZ( mat Y_k, uvec C_k, Rcpp::List particles );
  
  Rcpp::List sampleXi0( mat Y_k, uvec C_k, uvec N_k, Rcpp::List particles );
  
  Rcpp::List sampleE( uvec N_k, Rcpp::List particles, Rcpp::List prior );
    
  arma::vec logPriorDens( Rcpp::List particles,
                          Rcpp::List prior );

  arma::vec logPostDens( mat Y_k, uvec C_k, uvec N_k,
                         Rcpp::List particles, 
                         Rcpp::List prior );
    
  Rcpp::List iter(uvec T, int k, umat N, 
                  Rcpp::List particles, arma::mat log_dQ,
                  Rcpp::List prior);

    
  arma::uvec sampleT( arma::cube xi,
                      arma::cube Omega,
                      arma::mat alpha,
                      arma::mat logW );
    
public:
  // constructor 
  PMC( arma::mat Y,
       arma::uvec C,
       Rcpp::List prior,
       Rcpp::List pmc,
       Rcpp::List state,
       Rcpp::List initParticles, bool init );
  
  Rcpp::List get_chain();
};

#endif