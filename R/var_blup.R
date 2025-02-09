var_blup <- function(data, model) {
  # Extract variances parameters from the model
  sig2_u <- as.data.frame(lme4::VarCorr(model))$vcov[1]
  sig2_e <- as.data.frame(lme4::VarCorr(model))$vcov[2]
  
  # Extract model matrices
  Z <- lme4::getME(model, "Z")
  Zt <- lme4::getME(model, "Zt")
  X <- lme4::getME(model, "X")
  
  D <- sig2_u * diag(length(unique(data$id)))
  R <- sig2_e * diag(nrow(Z))
  
  V <- Z %*% D %*% Matrix::t(Z) + R
  invV <- Matrix::chol2inv(chol(V))
  
  XWX <- Matrix::t(X) %*% invV %*% X
  inv_XWX <- Matrix::chol2inv(chol(XWX))
  
  vbi  <- D %*% Matrix::t(Z) %*% (invV - invV %*% X %*% inv_XWX %*% Matrix::t(X) %*% invV) %*% Z %*% D
  
  return(data.frame(id = unique(data$id), v_blup = Matrix::diag(D - vbi)))
}

