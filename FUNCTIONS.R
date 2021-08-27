
## FUNCTIONS TO IMPLEMENT THE POWER CACULATION METHOD

var.logistic <- function(I,T,q,beta,theta,trtSeq,n,invR,Z_part)
{
  sumZWZ <- matrix(0,ncol=T+1,nrow=T+1)
  for (s in 1:(I/q)) { 
    for (m in 1:q) {
      u_c <- c(plogis(beta+theta*trtSeq[s,]))
      u <- rep(u_c,each=n)
      A_half <- diag(sqrt(u*(1-u)))
      W <- A_half %*% invR %*% A_half
      Z <- cbind(Z_part,c(rep(0,s*n),rep(1,(T-s)*n)))
      ZWZ <- t(Z) %*% W %*% Z
      sumZWZ <- sumZWZ + ZWZ
    } 
  }
  invsumZWZ <- solve(sumZWZ) 
  var_logistic <- invsumZWZ[T+1,T+1] 
  return(var_logistic)
}

var.logistic.trt.het <- function(I,T,q,beta,theta,trtSeq,n,rho1,rho2,rho3,Z_part)
{
  sumZWZ <- matrix(0,ncol=T+1,nrow=T+1)
  for (s in 1:(I/q)) { 
    for (m in 1:q) {
      u_c <- c(plogis(beta+theta*trtSeq[s,]))
      u <- rep(u_c,each=n)
      A_half <- diag(sqrt(u*(1-u)))
      R <- (1-rho1)*diag(n*T) + 
           (rho1-rho2)*kronecker((1-trtSeq[s,])%*%t(1-trtSeq[s,]),matrix(1,n,n)) + 
           rho2*matrix(1,n*T,n*T) +
           (rho3-rho2)*kronecker(trtSeq[s,]%*%t(trtSeq[s,]),matrix(1,n,n)) + 
           (rho1-rho3)*diag(kronecker(trtSeq[s,],rep(1,n)))  # heterogeneous treatment correlation matrix
      invR <- solve(R)
      W <- A_half %*% invR %*% A_half
      Z <- cbind(Z_part,c(rep(0,s*n),rep(1,(T-s)*n)))
      ZWZ <- t(Z) %*% W %*% Z
      sumZWZ <- sumZWZ + ZWZ
    } 
  }
  invsumZWZ <- solve(sumZWZ) 
  var_logistic <- invsumZWZ[T+1,T+1] 
  return(var_logistic)
}

var.logistic.vary <- function(I,T,q,beta,theta,S.ij,n,rho)
{
  sumZWZ <- matrix(0,ncol=T+1,nrow=T+1)
  for (i in 1:I) { 
    u_c <- c(plogis(beta+theta*S.ij[i,]))
    u <- rep(u_c,each=n[i])
    A_half <- diag(sqrt(u*(1-u)))
    R <- (1-rho)*diag(n[i]*T) + rho*matrix(1,nrow=n[i]*T,ncol=n[i]*T)  # exchangeable correlation matrix
    invR <- solve(R)
    W <- A_half %*% invR %*% A_half
    Z_part <- kronecker(diag(T),rep(1,n[i]))
    Z <- cbind(Z_part,rep(S.ij[i,],each=n[i]))
    ZWZ <- t(Z) %*% W %*% Z
    sumZWZ <- sumZWZ + ZWZ
  }
  invsumZWZ <- solve(sumZWZ) 
  var_logistic <- invsumZWZ[T+1,T+1] 
  return(var_logistic)
}

## FUNCTIONS TO SIMULATE DATA
# Reference: 
# Qaqish, B. F. (2003). A family of multivariate binary distributions for simulating correlated
# binary variables with specified marginal means and correlations. Biometrika 90, 455–463.

create.B <- function(I,T,q,beta,theta,trtSeq,n,R)
{
  B <- NULL
  for (i in 1:(I/q)) {
    u_c <- c(plogis(beta+theta*trtSeq[i,]))
    u <- rep(u_c,each=n)
    A_half <- diag(sqrt(u*(1-u)))
    v <- A_half %*% R %*% A_half # covariance matrix
    b <- v
    for (f in 2:(n*T)) {  
      f1 <- f-1
      gf <- v[1:f1,1:f1]
      sf <- v[1:f1,f]
      bf <- solve(gf,sf)  # b as defined in Qaqish equation 3
      b[1:f1,f] <- bf
    }
    B <- cbind(B,b)
  }
  return(B)
}

create.B.het <- function(I,T,q,beta,theta,trtSeq,n,rho1,rho2,rho3)
{
  B <- NULL
  for (i in 1:(I/q)) {
    u_c <- c(plogis(beta+theta*trtSeq[i,]))
    u <- rep(u_c,each=n)
    A_half <- diag(sqrt(u*(1-u)))
    R <- (1-rho1)*diag(n*T) + 
      (rho1-rho2)*kronecker((1-trtSeq[i,]) %*% t(1-trtSeq[i,]),matrix(1,n,n)) + 
      rho2*matrix(1,n*T,n*T) +
      (rho3-rho2)*kronecker(trtSeq[i,] %*% t(trtSeq[i,]),matrix(1,n,n)) + 
      (rho1-rho3)*diag(kronecker(trtSeq[i,],rep(1,n)))   # correlation matrix
    v <- A_half %*% R %*% A_half      # covariance matrix
    b <- v
    for (f in 2:(n*T)) {              # prepare coeffs
      f1 <- f-1
      gf <- v[1:f1,1:f1]
      sf <- v[1:f1,f]
      bf <- solve(gf,sf)
      b[1:f1,f] <- bf
    }
    B <- cbind(B,b)
  }
  return(B)
}

create.response <- function(I,T,q,beta,theta,trtSeq,n,B)
{
  y_out <- NULL
  for (i in 1:(I/q)) {
    u_c <- c(plogis(beta+theta*trtSeq[i,]))
    u <- rep(u_c,each=n)
    b <- B[,((i-1)*n*T+1):(i*n*T)]
    y_temp <- matrix(-1,n*T,q)
    for(h in 1:q) {                   # simulate data matrix
      y <- rep(-1,(n*T))
      y[1] <- rbinom(1,1,u[1])
      for (l in 2:(n*T)) {
        l1 <- l-1
        res <- y[1:l1] - u[1:l1]          # residuals
        cl <- u[l] + sum(res*b[1:l1,l])   # cond.mean (Qaqish equation 3)
        y[l] <- rbinom(1,1,cl)
      }
      y_temp[,h] <- y
    }
    y_out <- cbind(y_out,y_temp)         
  }
  y <- c(y_out)
  return(y)
}

## FUNCTIONS TO FIT THE EXPONENTIAL DECAY CORRELATION STRUCTURE WITH GEE 1.5
# References:
# Prentice, R. L. (1988). Correlated binary regression with covariates specific to each binary
# observation. Biometrics 44, 1033–1048.
# Mancl, L. A. and DeRouen, T. A. (2001). A covariance estimator for GEE with improved small-sample 
# properties. Biometrics 57, 126–134.
# Li, F., Turner, E. L., and Preisser, J. S. (2018). Sample size determination for GEE analyses of 
# stepped wedge cluster randomized trials. Biometrics 74, 1450–1458.

FITEXPDECAY <- function(y, X, clsize, clpersize, T, maxiter, epsilon){
  require(MASS)
  
  # Number of parameters in the mean model
  p <- ncol(X)
  # Number of correlation parameters
  g <- 2
  
  # Set initial Fisher scoring iteration increment higher than the convergence threshold 
  delta <- rep(2*epsilon,p)
  delta_A <- rep(2*epsilon,g)
  
  # Initial alpha (correlation parameters)
  alpha <- c(0.1,0.1)
  
  # Initial beta estimate (mean model parameters)
  INITRES <- INITBETA(y,X)
  beta <- INITRES$beta
  Ustar <- INITRES$Ustar
  
  # Fisher scoring iterations (Prentice equation 17)
  niter <- 1
  while((niter<=maxiter) & (max(abs(delta_A))>epsilon | max(abs(delta))>epsilon)){
    # estimate alpha 
    SCORE_A <- SCOREALPHA(beta, alpha, y, X, clsize, clpersize, T, g)
    U_A <- SCORE_A$U
    Ustar_A <- SCORE_A$Ustar
    delta_A <- try(solve(Ustar_A,U_A))
    # correlation parameter range check
    ok <- SCORE_A$ok
    if (('try-error' %in% class(delta_A)==1)|(ok %in% c(0,NA))) {
      alpha <- 0.95*alpha
      delta_A <- rep(2*epsilon,g)
    }
    else {
      alpha <- alpha+delta_A
    }
    
    # estimate beta 
    SCORE_RES <- SCORE(beta, alpha, y, X, clsize, clpersize, T, p)
    U <- SCORE_RES$U
    Ustar <- SCORE_RES$Ustar
    UUtran <- SCORE_RES$UUtran
    delta <- solve(Ustar,U)
    beta <- beta+delta
    converge <- (max(abs(delta))<=epsilon*10 & max(abs(delta_A))<=epsilon*10)
    
    # counter
    niter <- niter+1
  }
  
  # Output results when converged
  if(converge==1){
    # Model-based variance estimate
    model_based <- ginv(Ustar)
    # Sandwich variance estimate    
    sandwich <- model_based%*%UUtran%*%t(model_based)
    # Variance estimate from Mancl and DeRouen  
    UUbc <- matrix(0,p,p)
    locx <- BEGINEND(clsize)
    for (i in 1:length(clsize)){
      X_c <- X[locx[i,1]:locx[i,2],,drop=FALSE]
      y_c <- y[locx[i,1]:locx[i,2]]
      mu_c <- 1/(1+exp(c(-X_c%*%beta)))
      
      C <- X_c*(mu_c*(1-mu_c))
      INVB <- INVBMAT(mu_c,alpha,clpersize[i],T)
      H <- C%*%model_based%*%t(C)%*%INVB
      corr_fact <- solve(diag(clsize[i])-H)
      
      U_c <- t(C)%*%INVB%*%corr_fact%*%(y_c-mu_c)
      UUbc_c <- tcrossprod(U_c)
      UUbc <- UUbc + UUbc_c 
    }
    MDvar <- model_based%*%UUbc%*%t(model_based)
  }
  
  coefficients <- cbind(beta, sqrt(diag(model_based)), sqrt(diag(sandwich)), sqrt(diag(MDvar)))
  rownames(coefficients) <- c(rep("time",length(beta)-1),"trt")
  colnames(coefficients) <- c("beta estimates","model-based SE","sand SE","MD SE")
  rownames(alpha) <- c("tau","rho")
  colnames(alpha) <- c("alpha estimates")
  model_fit <- c(as.character(converge),as.character(niter),as.character(ok))
  names(model_fit) <- c("model converged", "number of iterations", "corr parameters in range")
  return(list(coefficients=coefficients,corr_parameters=alpha,model_fit=model_fit))
}

# Generates initial values for beta
# Approximates logistic regression using Newton's method
INITBETA <- function(y,X){
  z <- y+y-1
  beta <- solve(t(X)%*%X,t(X)%*%z)
  
  for (i in 1:2){
    u <- c(X%*%beta)
    u <- 1/(1+exp(-u))
    v <- u*(1-u)
    z <- t(X)%*%(y-u)
    Ustar <- t(X)%*%(X*v)
    d <- solve(Ustar,z)
    beta <- beta+d
  }
  return(list(beta=c(beta),Ustar=Ustar))
}

# Find the begining and end of each cluster
BEGINEND<-function(clsize){
  last<-cumsum(clsize)
  first<-last-clsize+1
  return(cbind(first,last))
}

# Get INVB=INVAhalf%*%INVR%*%INVAhalf
INVBMAT <- function(mu_c,alpha,n,T){
  alpha0 <- alpha[1]
  alpha1 <- alpha[2]
  Fmat <- alpha1^abs(outer(0:(T-1),0:(T-1),"-"))
  R <- (1-alpha0)*diag(n*T) + alpha0*kronecker(Fmat,matrix(1,n,n))
  INVR <- solve(R)
  INVAhalf <- diag(1/sqrt(mu_c*(1-mu_c)))
  INVB <- INVAhalf%*%INVR%*%INVAhalf
  return(INVB)
}

# Get components for the Fisher scoring iteration for beta
SCORE<-function(beta, alpha, y, X, clsize, clpersize, T, p){
  U <- rep(0,p)
  UUtran <- matrix(0,p,p)
  Ustar <- matrix(0,p,p)
  locx <- BEGINEND(clsize)
  
  for(i in 1:length(clsize)){
    X_c <- X[locx[i,1]:locx[i,2],,drop=FALSE]
    y_c <- y[locx[i,1]:locx[i,2]]
    
    U_c <- rep(0,p)
    Ustar_c <- matrix(0,p,p)
    mu_c <- 1/(1+exp(c(-X_c%*%beta)))
    
    C <- X_c*(mu_c*(1-mu_c))
    INVB <- INVBMAT(mu_c,alpha,clpersize[i],T)
    
    U_c <- t(C)%*%INVB%*%(y_c-mu_c)
    UUtran_c <- tcrossprod(U_c)
    Ustar_c <- t(C)%*%INVB%*%C
    
    U <- U+U_c
    UUtran <- UUtran+UUtran_c
    Ustar <- Ustar+Ustar_c
  }
  return(list(U=U,UUtran=UUtran,Ustar=Ustar))
}

# Get components for the Fisher scoring iteration for alpha
SCOREALPHA<-function(beta, alpha, y, X, clsize, clpersize, T, g){
  U <- rep(0,g)
  Ustar <- matrix(0,g,g)
  locx <- BEGINEND(clsize)
  tau <- alpha[1]
  rho <- alpha[2]
  
  for (i in 1:length(clsize)){
    X_c <- X[locx[i,1]:locx[i,2],,drop=FALSE]
    y_c <- y[locx[i,1]:locx[i,2]]
    
    U_c <- rep(0,g)
    Ustar_c <- matrix(0,g,g)
    mu_c <- 1/(1+exp(c(-X_c%*%beta)))
    
    SQVARFUN <- sqrt(mu_c*(1-mu_c))
    INVSQVAR <- 1/SQVARFUN
    RX <- (y_c-mu_c)*INVSQVAR
    Gi <- tcrossprod(RX)
    Z_c <- Gi[upper.tri(Gi,diag=FALSE)]
    Fmat <- rho^abs(outer(0:(T-1),0:(T-1),"-"))
    Lmat <- (1-tau)*diag(clsize[i]) + tau*kronecker(Fmat,matrix(1,clpersize[i],clpersize[i]))
    rho_c <- Lmat[upper.tri(Lmat,diag=FALSE)]
    n_star <- choose(clsize[i],2)
    INVW_c <- diag(n_star)
    DFmat <- abs(outer(0:(T-1),0:(T-1),"-"))*(rho^(abs(outer(0:(T-1),0:(T-1),"-"))-matrix(1,T,T)))
    DLmat_tau <- -diag(clsize[i]) + kronecker(Fmat,matrix(1,clpersize[i],clpersize[i]))
    DLmat_rho <- tau*kronecker(DFmat,matrix(1,clpersize[i],clpersize[i]))
    E_c <- matrix(0,n_star,2)
    E_c[,1] <- DLmat_tau[upper.tri(DLmat_tau,diag=FALSE)]
    E_c[,2] <- DLmat_rho[upper.tri(DLmat_rho,diag=FALSE)]
    
    U_c <- t(E_c)%*%INVW_c%*%(Z_c-rho_c)
    Ustar_c <- t(E_c)%*%INVW_c%*%E_c
    
    U <- U+U_c
    Ustar <- Ustar+Ustar_c
    
    # Range checks (Prentice Discussion)
    mu_star1 <- mu_c/(1-mu_c)
    mu_star2 <- (1-mu_c)/mu_c
    term1 <- -sqrt(tcrossprod(mu_star1,mu_star1))
    term2 <- -sqrt(tcrossprod(mu_star2,mu_star2))
    term3 <- sqrt(tcrossprod(mu_star1,mu_star2))
    term4 <- sqrt(tcrossprod(mu_star2,mu_star1))
    max <- (term1 + term2 + abs(term1-term2))/2
    min <- (term3 + term4 - abs(term3-term4))/2
    max_c <- max[upper.tri(max,diag=FALSE)]
    min_c <- min[upper.tri(min,diag=FALSE)]
    inrange <- ((max_c <= rho_c) & (rho_c <= min_c))
    ok <- (sum(inrange) == n_star)
  }
  return(list(U=U,Ustar=Ustar,ok=ok))
}

## FUNCTIONS TO CALCULATE THE BIAS-CORRECTED MANCL DEROUEN SANDWICH VARIANCE
# Reference:
# Mancl, L. A. and DeRouen, T. A. (2001). A covariance estimator for GEE with improved small-sample 
# properties. Biometrics 57, 126–134.

md.var <- function(fit_pack,simdata_bin,var)
{
  beta_est <- fit_pack$coefficient
  m <- model.frame(y~treatment+factor(period), simdata_bin)
  mat <- as.data.frame(model.matrix(y~treatment+factor(period), m))
  mat$subj <- simdata_bin$cluster
  step11 <- matrix(0, nrow=length(beta_est), ncol=length(beta_est))
  step12 <- matrix(0, nrow=length(beta_est), ncol=length(beta_est))
  for (i in 1:I) {
    y <- as.matrix(simdata_bin$y[simdata_bin$cluster == unique(simdata_bin$cluster)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])],mat$subj == unique(simdata_bin$cluster)[i]))
    D <- matrix(0, ncol=ncol(covariate), nrow=nrow(covariate))
    for (j in 1:ncol(covariate)) {
      D[, j] <- covariate[,j]*exp(covariate %*% beta_est)/((1 + exp(covariate %*% beta_est))^2)
    }
    Vi <- diag(sqrt(c(exp(covariate %*% beta_est)/(1 + exp(covariate %*% beta_est))^2)), n*T) %*% var %*% 
      diag(sqrt(c(exp(covariate %*% beta_est)/(1 + exp(covariate %*% beta_est))^2)),n*T)
    xx <- t(D) %*% solve(Vi) %*% D
    step11 <- step11 + xx
  }
  for (i in 1:I) {
    y <- as.matrix(simdata_bin$y[simdata_bin$cluster == unique(simdata_bin$cluster)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])],mat$subj == unique(simdata_bin$cluster)[i]))
    D <- matrix(0, ncol=ncol(covariate), nrow=nrow(covariate))
    for (j in 1:ncol(covariate)) {
      D[, j] <- covariate[,j]*exp(covariate %*% beta_est)/((1 + exp(covariate %*% beta_est))^2)
    }
    Vi <- diag(sqrt(c(exp(covariate %*% beta_est)/(1 + exp(covariate %*% beta_est))^2)), n*T) %*% var %*% 
      diag(sqrt(c(exp(covariate %*% beta_est)/(1 + exp(covariate %*% beta_est))^2)), n*T)
    xy <- t(D) %*% solve(Vi) %*% solve(diag(n*T) - D %*% solve(step11) %*% t(D) %*% solve(Vi)) %*% 
      (y - exp(covariate %*% beta_est)/(1 + exp(covariate %*% beta_est)))
    step12 <- step12 + xy %*% t(xy)
  }
  cov.beta <- solve(step11) %*% (step12) %*% solve(step11)
  return(cov.beta[2,2]) 
}

md.var.het <- function(fit_pack,simdata_bin,alpha1,alpha2,alpha3,n,T,q,trtSeq)
{
  beta_est <- fit_pack$coefficient
  m <- model.frame(y~treatment+factor(period), simdata_bin)
  mat <- as.data.frame(model.matrix(y~treatment+factor(period), m))
  mat$subj <- simdata_bin$cluster
  step11 <- matrix(0, nrow=length(beta_est), ncol=length(beta_est))
  step12 <- matrix(0, nrow=length(beta_est), ncol=length(beta_est))
  for (i in 1:I) {
    var <- (1-alpha1)*diag(n*T) + 
      (alpha1-alpha2)*kronecker((1-trtSeq[ceiling(i/q),]) %*% t(1-trtSeq[ceiling(i/q),]),matrix(1,n,n)) + 
      alpha2*matrix(1,n*T,n*T) +
      (alpha3-alpha2)*kronecker(trtSeq[ceiling(i/q),] %*% t(trtSeq[ceiling(i/q),]),matrix(1,n,n)) + 
      (alpha1-alpha3)*diag(kronecker(trtSeq[ceiling(I/q),],rep(1,n)))   
    y <- as.matrix(simdata_bin$y[simdata_bin$cluster == unique(simdata_bin$cluster)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])],mat$subj == unique(simdata_bin$cluster)[i]))
    D <- matrix(0, ncol=ncol(covariate), nrow=nrow(covariate))
    for (j in 1:ncol(covariate)) {
      D[, j] <- covariate[,j]*exp(covariate %*% beta_est)/((1 + exp(covariate %*% beta_est))^2)
    }
    Vi <- diag(sqrt(c(exp(covariate %*% beta_est)/(1 + exp(covariate %*% beta_est))^2)), n*T) %*% var %*% 
      diag(sqrt(c(exp(covariate %*% beta_est)/(1 + exp(covariate %*% beta_est))^2)),n*T)
    xx <- t(D) %*% solve(Vi) %*% D
    step11 <- step11 + xx
  }
  for (i in 1:I) {
    var <- (1-alpha1)*diag(n*T) + 
      (alpha1-alpha2)*kronecker((1-trtSeq[ceiling(i/q),]) %*% t(1-trtSeq[ceiling(i/q),]),matrix(1,n,n)) + 
      alpha2*matrix(1,n*T,n*T) +
      (alpha3-alpha2)*kronecker(trtSeq[ceiling(i/q),] %*% t(trtSeq[ceiling(i/q),]),matrix(1,n,n)) + 
      (alpha1-alpha3)*diag(kronecker(trtSeq[ceiling(I/q),],rep(1,n)))   
    y <- as.matrix(simdata_bin$y[simdata_bin$cluster == unique(simdata_bin$cluster)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])],mat$subj == unique(simdata_bin$cluster)[i]))
    D <- matrix(0, ncol=ncol(covariate), nrow=nrow(covariate))
    for (j in 1:ncol(covariate)) {
      D[, j] <- covariate[,j]*exp(covariate %*% beta_est)/((1 + exp(covariate %*% beta_est))^2)
    }
    Vi <- diag(sqrt(c(exp(covariate %*% beta_est)/(1 + exp(covariate %*% beta_est))^2)), n*T) %*% var %*% 
      diag(sqrt(c(exp(covariate %*% beta_est)/(1 + exp(covariate %*% beta_est))^2)), n*T)
    xy <- t(D) %*% solve(Vi) %*% solve(diag(n*T) - D %*% solve(step11) %*% t(D) %*% solve(Vi)) %*% 
      (y - exp(covariate %*% beta_est)/(1 + exp(covariate %*% beta_est)))
    step12 <- step12 + xy %*% t(xy)
  }
  cov.beta <- solve(step11) %*% (step12) %*% solve(step11)
  return(cov.beta[2,2]) 
}
