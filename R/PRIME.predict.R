PRIME.predict <- function( X_test, coef, Cmat, weight, model_structure, method = c('PRIME', 'PRIME-MA'),
                           intercept = F, L = NULL)
{
  if (is.null(X_test))
    stop("\"X_test\" must be given!")
  if (is.null(n <- nrow(X_test)))
    stop("'X_test' must be a matrix")
  if (n == 0L)
    stop("0 (non-NA) cases")
  if (!is.matrix(X_test))
    X_test <- as.matrix(X_test)
  if (is.null(L))
    L <- 3
  p <- ncol(X_test)
  if (method == "PRIME") {
    if (missing(model_structure))
      stop("Users must supply a nonlinear index of model structure.")
    XX_test  <- non.spline(X_test, model_structure, L)
    Xs_test <- XX_test$Xs
    Y_test <- Xs_test %*% coef
  }

  if (method=='PRIME-MA') {
    if (missing(Cmat)|is.null(Cmat))
      stop("Users must supply a coefficient matrix of candidate models.")
    if (missing(weight)|is.null(weight))
      stop("Users must supply a weight.")

    mu_test <- matrix(NA, nrow(X_test), ncol(X_test))
    for(j in 1:p){
      if(intercept){
        Xts <- cbind(1, X_test[,-j], bs(X_test[,j], L))
        mu_test[,j] <- Xts %*% Cmat[[j]]
      }

      if(!intercept){
        Xts <- cbind(X_test[,-j], bs(X_test[,j], L))
        mu_test[,j] <- Xts %*% Cmat[[j]]
      }
    }
    Y_test <- mu_test %*% weight
  }
  return(list(Y_test = Y_test))
}

