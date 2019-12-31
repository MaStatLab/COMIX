calibrate <- function(x)
{
  C = x$data$C - 1
  Z = x$chain$t - 1
  if (!is.matrix(Z)) Z=matrix(Z, ncol=1)
  ns  = dim(x$chain$xi0)[3]
  K = dim(x$chain$xi0)[2]

  output = calib(x$data$Y,
                 matrix(C,ncol=1),
                 Z,
                 x$chain$xi, dim(x$chain$xi),
                 x$chain$xi0, dim(x$chain$xi0) )
  colnames(output$Y_cal) = colnames(x$data$Y)
  return(output)

}

calibrateNoDist <- function(x)
{
  C = x$data$C - 1
  Z = x$chain$t - 1
  if (!is.matrix(Z)) Z=matrix(Z, ncol=1)
  ns  = dim(x$chain$xi0)[3]
  K = dim(x$chain$xi0)[2]
  
  output = calibNoDist(x$data$Y,
                 matrix(C,ncol=1),
                 Z,
                 x$chain$xi, dim(x$chain$xi),
                 x$chain$xi0, dim(x$chain$xi0) )
  colnames(output$Y_cal) = colnames(x$data$Y)
  return(output)
  
}

relabelChain = function(res) {
  res$chain$t = res$chain$t - 1
  relabeled_chain = relabel(res)
  res$chain = relabeled_chain 
  res
}

# Recover stochastic representation parameters from distribution parameters:
transform_params = function(Omega, alpha) {
  n = NROW(Omega)
  m = NCOL(Omega)
  if (n!=m) stop("Omega is not sqaure")
  if (length(alpha)!=n) stop("alpha is of wrong length")
  if (n==1) {
    omega = sqrt(Omega)
  } else {
    omega = sqrt(diag(diag(Omega)))
  }
  omega_inv = solve(omega)
  OmegaBar = omega_inv %*% Omega %*% omega_inv
  alpha = matrix(alpha, ncol=1)
  numer = OmegaBar %*% alpha
  denom = as.numeric(sqrt(1 + t(alpha) %*% OmegaBar %*% alpha))
  delta = numer/denom
  if (n==1) {
    psi = sqrt(Omega) * delta
  } else {
    psi = sqrt(diag(Omega)) * delta
  }
  alpha.tilde = delta / sqrt(1-delta^2)
  if (n==1) {
    inv.Delta = 1/c(sqrt(1-delta^2))
  } else {
    inv.Delta = diag(1/c(sqrt(1-delta^2)))
  }
  Omega.epsilon = inv.Delta %*% OmegaBar %*% inv.Delta - alpha.tilde%*%t(alpha.tilde)
  return(list(delta=delta, omega=omega, Omega.epsilon=Omega.epsilon,
              psi=psi, Sigma=(Omega - psi%*%t(psi))))
}


Mode = function(x) {
  ux = unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

summarizeChain = function( res ) {
  chainSummary = list()
  K = res$prior$K
  chain = res$chain
  p = ncol(res$data$Y)
  J = length(unique(res$data$C))
  ns = res$pmc$nsave
  
  xi_raw = matrix(0, nrow=J, ncol=p*K)
  chainSummary$xi0 = matrix(0, nrow=p,ncol=K)
  chainSummary$psi = matrix(0, nrow=p,ncol=K)
  Omega_raw = matrix(0, nrow=p, ncol=p*K)
  Sigma_raw = matrix(0, nrow=p, ncol=p*K)
  E_raw = matrix(0, nrow=p, ncol=p*K)
  chainSummary$alpha = matrix(0, nrow=p,ncol=K)
  chainSummary$W = matrix(0, nrow=J, ncol=K)
  for (i in 1:ns) {
    xi_raw = xi_raw + chain$xi[,,i]
    chainSummary$xi0 = chainSummary$xi0 + chain$xi0[,,i]
    chainSummary$psi = chainSummary$psi + chain$psi[,,i]
    Omega_raw = Omega_raw + chain$Omega[,,i]
    Sigma_raw = Sigma_raw + chain$G[,,i]
    E_raw = E_raw + chain$E[,,i]
    chainSummary$alpha = chainSummary$alpha + chain$alpha[,,i]
    chainSummary$W = chainSummary$W + chain$W[,,i]
  }
  
  chainSummary$W = chainSummary$W/ns
  
  xi_raw = xi_raw/ns
  chainSummary$xi = array(0, dim=c(J,p,K))
  for (k in 1:K) {
    chainSummary$xi[,,k] = xi_raw[,(1+p*(k-1)):(p*k)]
  }
  
  chainSummary$xi0 = chainSummary$xi0/ns
  chainSummary$psi = chainSummary$psi/ns
  
  Omega_raw = Omega_raw/ns
  chainSummary$Omega = array(0, dim=c(p,p,K))
  for (k in 1:K) {
    chainSummary$Omega[,,k] = Omega_raw[,(1+p*(k-1)):(p*k)]
  }
  
  Sigma_raw = Sigma_raw/ns
  chainSummary$Sigma = array(0, dim=c(p,p,K))
  for (k in 1:K) {
    chainSummary$Sigma[,,k] = Sigma_raw[,(1+p*(k-1)):(p*k)]
  }
  
  E_raw = E_raw/ns
  chainSummary$E = array(0, dim=c(p,p,K))
  for (k in 1:K) {
    chainSummary$E[,,k] = E_raw[,(1+p*(k-1)):(p*k)]
  }
  
  chainSummary$alpha = chainSummary$alpha/ns
  
  chainSummary$meanvec = array(0, c(J, p, K))
  chainSummary$meanvec0 = matrix(0, p, K)
  for (k in 1:K) {
    del.om = transform_params(chainSummary$Omega[,,k], chainSummary$alpha[,k])
    # chainSummary$psi[,k] = del.om$psi
    chainSummary$Sigma[,,k] = del.om$Sigma
    for (j in 1:J) {
      chainSummary$meanvec[j,,k] = chainSummary$xi[j,,k] + del.om$omega %*% del.om$delta*sqrt(2/pi)
    }
    chainSummary$meanvec0[,k] = chainSummary$xi0[,k] + del.om$omega %*% del.om$delta*sqrt(2/pi)
  }
  
  chainSummary$t = apply(chain$t,2,Mode)
  chainSummary$a0 = mean(chain$a0)
  
  chainSummary
}