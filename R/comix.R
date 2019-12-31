#'
#' This function generates a sample from the posterior of COMIX.
#'
#' @param Y Matrix of the data. Each row represents an observation.
#' @param C Vector of the group label of each observation. Labels must be integers starting from 1.
#' @param prior A list giving the prior information. If unspecified, a default prior is used.
#' The list includes the following hyparameters:
#' \code{zeta} Coarsening parameter. A number between 0 and 1. zeta = 1: sample from standard posterior;
#' zeta < 1: sample from power posterior. The lower zeta is, the more flexible the kernels become.
#' \code{K} Maximal number of mixture components.
#' \code{tau_a} Parameters for gamma prior for concentration parameter of the stick breaking process
#' prior for the weights.
#' m0 = ncol(Y) + 2,
#' Lambda = stats::cov(Y),
#' b0 = colMeans(Y),
#' B0 = 100 * stats::cov(Y),
#' e0 = ncol(Y) + 2,
#' E0 = 0.1 * stats::cov(Y),
#' \code{merge_step} Introduce step to merge mixture components with small KL divergence. Default 
#' is \code{merge_step = TRUE}.
#' \code{merge_par} Parameter controlling merging radius. Default is \code{merge_par = 0.1}.
#' @param prior A list giving the
#' @param state A list giving the
#' @return 
#' @export
comix = function(Y, C, prior = NULL, pmc = NULL, state = NULL)
{
  Y = as.matrix(Y)
  p = ncol(Y)

  # R wrapper (this code is general):
  if(is.null(prior)) {
    prior = list(    zeta = 1,
                     K = 10,
                     tau_a = c(1, 1),
                     m0 = ncol(Y) + 2,
                     Lambda = stats::cov(Y),
                     b0 = colMeans(Y),
                     B0 = 100 * stats::cov(Y),
                     e0 = ncol(Y) + 2,
                     E0 = 0.1 * stats::cov(Y),
                     merge_step = TRUE,
                     merge_par = 0.1)
  } else {
    if(is.null(prior$zeta)) 
      prior$zeta = 1;
    if(is.null(prior$K)) 
      prior$K = 10;
    if(is.null(prior$tau_a))
      prior$tau_a = c(1, 1);
    if(is.null(prior$Lambda))
      prior$Lambda = stats::cov(Y);
    if(is.null(prior$m0))
      prior$m0 = ncol(Y) + 2;
    if(is.null(prior$b0))
      prior$b0 = colMeans(Y);
    if(is.null(prior$B0))
      prior$B0 = 100*stats::cov(Y);
    if(is.null(prior$e0))
      prior$e0 = ncol(Y) + 2;
    if(is.null(prior$E0))
      prior$E0 = 0.1 * stats::cov(Y);
    if(is.null(prior$merge_step))
      prior$merge_step = TRUE;
    if(is.null(prior$merge_par))
      prior$merge_par = 0.1;
  }
  
  
  if(is.null(pmc)) {
    pmc = list(npart = 10, nburn = 1000, nsave = 1000, nskip = 1, ndisplay = 500)
  } else {
    if(is.null(pmc$npart))
      pmc$npart = 10
    if(is.null(pmc$nburn))
      pmc$nburn = 5000
    if(is.null(pmc$nburn))
      pmc$nburn = 5000
    if(is.null(pmc$nsave))
      pmc$nsave = 1000
    if(is.null(pmc$nskip))
      pmc$nskip = 1
    if(is.null(pmc$ndisplay))
      pmc$ndisplay = 100
  }

  if(is.null(state$t)) {
    state$t = stats::kmeans(Y, prior$K, iter.max = 100)$cluster - 1
  }

  J = length(unique(C))
  if( sum( sort(unique(C)) == 1:J )  != J )
  {
    print("ERROR: unique(C) should look like 1, 2, ...")
    return(0);
  }
  C = C - 1

  ans = perturbedSNcpp(Y, C, prior, pmc, state, initParticles = NULL, init=T)
  colnames(ans$data$Y) = colnames(Y)
  ans$data$C = ans$data$C + 1
  ans$chain$t = ans$chain$t + 1
  class(ans) = "COMIX"
  return(ans)
}
