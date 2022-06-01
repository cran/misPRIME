PRIME.impute <- function (ind, Xmat, Xmats, Y, IDX, IDXs, bw, k_type = NULL)
{
  p <- ncol(Xmats)
  if (is.null(k_type))
    k_type = "gaussian"
  ZI <- Xmat[ind, ]
  Pi <- which(!is.na(Xmat[ind, ]))
  #w <- sum(Xmat[ind, Pi] * beta[Pi])
  i.peer <- Reduce(intersect, IDX[Pi])
  Pi.peer <- which(is.na(Xmat[ind, ]))

  for(j in 1:p){
    ij.peer <- intersect(i.peer, IDXs[[j]])
    n_ijp <- length(ij.peer)
    Xx <- matrix(Xmat[ij.peer, Pi], nrow = length(ij.peer),
                 ncol = length(Pi))
    h <- apply((Xx + rnorm(n_ijp, sd = 1e-05)),2,sd) * bw
    ww <- Xmat[ind, Pi]
    uu <- Xmat[ij.peer, Pi] - ww
    Wij <- apply((apply(uu,2,kern)), 1, prod)
    ZI[j] <- sum(Xmats[ij.peer, j] * Wij)/sum(Wij)
  }
  return(ZI)
}
