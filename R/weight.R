cv.weight <- function(y, x, L){
  p <- as.numeric(ncol(x))
  dat.cc <- na.omit(cbind(y, x))
  y.cc <- dat.cc[, 1]
  x.cc <- dat.cc[, -1]
  n1 <- as.numeric(nrow(x.cc))
  e.cv <- matrix(0, nrow = n1, ncol = p)
  for (j in 1:p){
    xs <- cbind(x.cc[,-j], bs(x.cc[,j],L))
    xx_j <- ginv(t(xs) %*% xs)
    xy_j <- t(xs) %*% y.cc
    beta_j <- xx_j %*% xy_j
    mu_j <- xs %*% beta_j
    P_j <- xs %*% xx_j %*% t(xs)
    Q_j <- matrix(0, n1, n1)
    d_j <- 1 / (1 - diag(P_j))
    diag(Q_j) <- d_j
    e.cv[, j] <- Q_j %*% (y.cc - mu_j)
  }
  return(solve.weight(e.cv, rep(0, p), p))
}

cp.weight <- function(y, x, L){
  p <- as.numeric(ncol(x))
  dat.cc <- na.omit(cbind(y, x))
  y.cc <- dat.cc[, 1]
  x.cc <- dat.cc[, -1]
  n1 <- as.numeric(nrow(x.cc))
  tr <- vector()
  e <- matrix(NA,n1,p)
  P <- list()
  Omega <- matrix(NA,n1,p)
  for (j in 1:p){
    xs <- cbind(x.cc[,-j], bs(x.cc[,j],L))
    xx_j <- ginv(t(xs) %*% xs)
    xy_j <- t(xs) %*% y.cc
    beta_j <- xx_j %*% xy_j
    mu_j <- xs %*% beta_j
    P[[j]] <- xs %*% xx_j %*% t(xs)
    e_j <- y.cc - xs %*% beta_j
    e[,j] <- e_j
    Omega[,j] <- as.vector(e_j^2)

  }
  Omega_m <- diag(apply(Omega,1,mean))
  for(j in 1:p){
    tr[j] <- sum(diag(P[[j]] %*% Omega[,j]))
  }
  d <- (-2*tr)
  return(solve.weight(e, d, p))
}
