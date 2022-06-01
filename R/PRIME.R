#' @title Partial Replacement IMputation Estimation (PRIME) for Missing Covariates
#' @description partial replacement imputation estimation conducts imputation and regression simultaneously for missing covariates in additive partially linear model.
#'
#' @param Y a numeric vector, the response variable.
#' @param X a numeric matrix that may include NAs (missing), the covariate matrix.
#' @param method Users can choose \code{PRIME} or \code{PRIME-MA}. If \code{method="PRIME"}, users must provide the model structure (nonlinear part index) in the input argument; If \code{method=="PRIME-MA"}, then the program automatically applies model averaging methods to reduce  reduce the loss of misspecification of models without model structure.
#' @param model_structure only available when \code{method="PRIME"}. It is a 0/1 index vector representing whether each variable is linear/nonlinear in the partially linear model. For details see Example section.
#' @param intercept logical. if \code{TRUE}, an intercept is included in the basis; default is \code{FALSE}.
#' @param bw a positive value, specify the bandwidth in estimating missing values, default as \code{NULL}. When \code{bw=NULL}, it is automatically selected by Silverman's rule of thumb method.
#' @param k_type an optional character string, specify the type of kernel used in iterative estimating algorithm and support 'epk', 'biweight', 'triangle', 'gaussian', 'triweight', 'tricube', 'cosine', 'uniform' in current version, default as 'gaussian'.
#' @param weight_type Options for computing weights for \code{PRIME-MA} method. Users can choose among \code{CP} and \code{CV}.
#' @param L an optional positive integer, degree of the piecewise polynomial, default as '3' for cubic splines.
#'
#' @return  an object of class "prime" is a list containing at least the following components:
#' \item{coef}{only available when \code{method="PRIME"}. A vector of coefficients of partially linear model.}
#' \item{beta}{only available when \code{method="PRIME"}. A vector of coefficients of linear parts in partially linear model.}
#' \item{Cmat}{only available when \code{method="PRIME-MA"}. A list of coefficients of candidate partially linear models.}
#' \item{weight}{only available when \code{method="PRIME-MA"}. The weights for candidate models, each candidate model involves one nonlinear part and others are linear parts.}
#'
#' @export
#' @import MASS splines quadprog stats
#' @examples
#' data(PRIME_SimuData)
#' X = PRIME_SimuData[,-1]
#' Y = PRIME_SimuData[,1]
#' model_structure <- c(rep(0,5),1,1,1)
#'
#' # estimation
#' result <- PRIME(Y, X, method = 'PRIME', model_structure, intercept = FALSE, weight_type = 'CV')
#' result$coef
#' result$beta
PRIME <- function (Y, X, method = c('PRIME', 'PRIME-MA'), model_structure = NULL, intercept = FALSE, bw = NULL, k_type = NULL,
                   weight_type = c("CP", "CV"), L = NULL)
{
  if (is.null(Y) || is.null(X))
    stop("\"X\",\"Y\" must be given\n simutaneously!")
  if (is.null(n <- nrow(X)))
    stop("'X' must be a matrix")
  if (n == 0L)
    stop("0 (non-NA) cases")
  if (!is.null(bw) && !is.numeric(bw))
    stop("\"bw\" must be NULL or a positive scalar!")
  if (!is.matrix(X))
    X <- as.matrix(X)
  if (method == "PRIME") {
    if (missing(model_structure))
      stop("Users must supply a nonlinear index of model structure.")
    if (!is.null(model_structure)){
      if (!all(as.numeric(model_structure) %in% c(0, 1)))
        stop("There can only be 0 or 1 in model_structure")
    }
  }
  if (method=='PRIME-MA') {
    if (missing(weight_type))
      stop("Users must supply a weight type.")
    if (!all(weight_type %in% c('CP', 'CV')))
      stop("There can only be 'CP' or 'CV' for weight type")
  }

  p <- ncol(X)
  n <- nrow(X)

  if (method == 'PRIME'){
    XX  <- non.spline(X, model_structure, L)
    X <- XX$X
    Xs <- XX$Xs
    IDX <- indx.comp(X)
    IDXs <- indx.comp(Xs)
  }

  if (is.null(bw))
    bw <- 1.06*n^{-1/5}
  if (is.null(L))
    L <- 3

  NA_ind <- which(apply(is.na(X), 1, sum) != 0)


  if(method == 'PRIME'){

    Z1 <- Xs
    Z2 <- Xs
    for (i in NA_ind) {
      temp_Z1 <- PRIME.impute(ind = i, Xmat = X, Xmats= Xs,
                              Y = Y, IDX = IDX, IDXs=IDXs, bw = bw, k_type = k_type)
      if (is.null(temp_Z1)) {
        return(NULL)
      }
      else {
        Z1[i, ] <- temp_Z1

      }
    }
    Z2[is.na(Xs)] <- Z1[is.na(Xs)]
    if (intercept) {
      coef <- lm(Y ~ Z2)$coefficients
    }
    if (!intercept) {
      coef <- lm(Y ~ Z2 - 1)$coefficients
    }

    beta <- NULL
    if (sum(model_structure) > 0 & sum(model_structure) < p){
      lin_index <- 1:(p-sum(model_structure))
      beta <- coef[lin_index]
      lin_index0 <- which(model_structure==0)
      names(beta) <- paste('X',lin_index0, sep = '')
    }
    if (sum(model_structure) == 0 ){
      beta <- coef
      names(beta) <- paste('X',1:p, sep = '')
    }
    Cmat <- NULL
    weight <- NULL
  }

  if(method == 'PRIME-MA'){

    Cmat <- list()
    for(jj in 1:ncol(X)){
      Xs <- cbind(X[,-jj], bs(X[,jj],3))
      IDX <- indx.comp(X)
      Z1 <- Xs
      Z2 <- Xs
      IDXs <- indx.comp(Xs)
      for (i in NA_ind) {
        temp_Z1 <- PRIME.impute(ind = i, Xmat = X, Xmats= Xs,
                                Y = Y, IDX = IDX, IDXs=IDXs, bw = bw, k_type = k_type)
        if (is.null(temp_Z1)) {
          return(NULL)
        }
        else {
          Z1[i, ] <- temp_Z1
        }
      }
      Z2[is.na(Xs)] <- Z1[is.na(Xs)]
      if (intercept) {
        coef <- lm(Y ~ Z2)$coefficients
      }
      if (!intercept) {
        coef <- lm(Y ~ Z2 - 1)$coefficients
      }
      Cmat[[jj]] <- coef
    }
    if (weight_type == 'CV'){
      weight <- cv.weight(Y, X, L)
    }
    if (weight_type == 'CP'){
      weight <- cv.weight(Y, X, L)
    }
    coef <- NULL
    beta <- NULL
  }
  res <- list(coef = as.vector(coef), beta = beta, Cmat = Cmat, weight = weight)
  return(res)
}
