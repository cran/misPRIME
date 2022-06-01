non.spline <- function(X, model_structure, L){
  if(sum(model_structure) > 0 & sum(model_structure) < ncol(X)){
    # Reorder the variables so that the linear part comes first
    non_index0 <- (model_structure == 1)
    lin_index0 <- (model_structure == 0)
    X <- cbind(X[,lin_index0], X[, non_index0])
    p <- ncol(X)
    lin_num <- p-sum(model_structure)
    non_index <- (lin_num + 1):p
    X0 <- X[,-non_index]
    for(nn in non_index){
      Xs <- cbind(X0, bs(X[,nn],L))
      X0 <- Xs
    }
  }

  if(sum(model_structure) == 0){
    Xs <- X
  }

  if(sum(model_structure) == ncol(X)){
    X0 <- NULL
    for(nn in seq(ncol(X))){
      Xs <- cbind(X0, bs(X[,nn],L))
      X0 <- Xs
    }
  }
  return(list(X = X, Xs = Xs))
}

indx.comp <- function (X)
{
  p <- ncol(X)
  IDX <- NULL
  for (j in 1:p) {
    idxj <- list(which(!is.na(X[, j])))
    IDX <- c(IDX, idxj)
  }
  IDX
}

kern <- function (u, type = "gaussian")
{
  if (!is.element(type, c("epk", "biweight", "triangle",
                          "gaussian", "triweight", "tricube",
                          "cosine", "uniform",'eucle')))
    stop("type must belong to 'epk', 'biweight', 'triangle',\n         'guassian', 'triweight', 'tricube', 'cosine', 'uniform'!")
  if (type == "epk")
    f = 0.75 * (1 - u^2) * (u <= 1 & u >= -1)
  else if (type == "biweight")
    f = (15/16) * (1 - u^2)^2 * (u <= 1 & u >= -1)
  else if (type == "triangle")
    f = (1 - abs(u)) * (u <= 1 & u >= -1)
  else if (type == "gaussian")
    f = dnorm(u)
  else if (type == "triweight")
    f = 35/32 * (1 - u^2)^3 * (u <= 1 & u >= -1)
  else if (type == "tricube")
    f = 70/81 * (1 - abs(u)^3)^3 * (u <= 1 & u >= -1)
  else if (type == "cosine")
    f = pi/4 * cos(pi * u/2)
  else if (type == "uniform")
    f = 1/2 * (u <= 1 & u >= -1)
  return(f)
}

solve.weight <- function(Dm, dv, k){
  Dmat <- 2 * (t(Dm) %*% Dm)
  dvec <- dv

  A1 <- rep(1, k)
  A2 <- A3 <- matrix(0, k, k)
  diag(A2) <- 1
  diag(A3) <- -1
  Amat <- cbind(A1, A2, A3)
  bvec <- c(1, rep(0, k), rep(-1, k))

  sol <- solve.QP(Dmat, dvec, Amat, bvec = bvec, meq = 1)

  return(sol$solution)
}
