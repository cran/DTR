###################################################
### Reference:
### Tang X, Wahed AS: Comparison of treatment regimes with adjustment for auxiliary
### variables. Journal of Applied Statistics 38(12):2925-2938, 2011
###################################################

###################################################
### code chunk number 1: chunklibraries
###################################################

#Libraries required
require(survival)

###################################################
### code chunk number 2: chunkCox
###################################################

DTR.Cox.test <- function(fdata # A complete data frame representing the data for two-stage randomization designs
                               # fdata = data frame {X, TR, R, Z, U, delta, V}
                               # V represents covariates
                               # There could be no covariate, one covariate, and more than one covariates
) {
 
  #Retrieve data
  n <- nrow(fdata)
  X <- fdata$X # X=0 for A1, 1 for A2
  TR <- fdata$TR
  R <- fdata$R
  Z <- fdata$Z # Z=0 for B1, 1 for B2
  U <- fdata$U
  delta <- fdata$delta
  
  #Chek for errors
  if (is.null(X)) stop("X can not be empty") 
  if (is.null(TR)) stop("TR can not be empty")  
  if (is.null(R)) stop("R can not be empty")
  if (is.null(Z)) stop("Z can not be empty")  
  if (is.null(U)) stop("U can not be empty")  
  if (is.null(delta)) stop("delta can not be empty") 
  
  #Retrive covariates
  V <- fdata[, ! names(fdata) %in% c("X", "TR", "R", "Z", "U", "delta")]
  
  if(NCOL(V)==0) {
    
    #Format data
    data.model <- data.frame(
      t.start = c(rep(0, n), TR[which(TR<=U)]),
      t.end = c(U[which(TR>U)], TR[which(TR<=U)], U[which(TR<=U)]),
      t.delta = c(delta[which(TR>U)], rep(0, length(which(TR<=U))), delta[which(TR<=U)]),
      V.X = c(X[which(TR>U)], X[which(TR<=U)], X[which(TR<=U)]),
      V.R = c(rep(0, n), rep(1, length(which(TR<=U)))),
      V.XR = c(rep(0, n), X[which(TR<=U)]),
      V.RZ = c(rep(0, n), Z[which(TR<=U)]),
      V.XRZ = c(rep(0, n), X[which(TR<=U)]*Z[which(TR<=U)])
    )
    
    #Get rid of t.start=t.end=0
    if(length(which(data.model$t.end==0))>0) data.model <- data.model[which(data.model$t.end!=0),]
    
    #Fit the Cox proportional hazard model
    fit <- coxph(Surv(t.start, t.end, t.delta)~V.X + V.R + V.XR + V.RZ + V.XRZ, data=data.model)
        
  } 
    
  if(NCOL(V)==1) {
 
    #Format data
    data.model <- data.frame(
      t.start = c(rep(0, n), TR[which(TR<=U)]),
      t.end = c(U[which(TR>U)], TR[which(TR<=U)], U[which(TR<=U)]),
      t.delta = c(delta[which(TR>U)], rep(0, length(which(TR<=U))), delta[which(TR<=U)]), 
      V.X = c(X[which(TR>U)], X[which(TR<=U)], X[which(TR<=U)]),
      V.R = c(rep(0, n), rep(1, length(which(TR<=U)))),
      V.XR = c(rep(0, n), X[which(TR<=U)]),
      V.RZ = c(rep(0, n), Z[which(TR<=U)]),
      V.XRZ = c(rep(0, n), X[which(TR<=U)]*Z[which(TR<=U)]), 
      V.V = c(V[which(TR>U)], V[which(TR<=U)], V[which(TR<=U)])
    )
  
    #Get rid of t.start=t.end=0
    if(length(which(data.model$t.end==0))>0) data.model <- data.model[which(data.model$t.end!=0),]
  
    #Fit the Cox proportional hazard model
    fit <- coxph(Surv(t.start, t.end, t.delta)~V.X + V.R + V.XR + V.RZ + V.XRZ + V.V, data=data.model)
  
  }
  
  if(NCOL(V)>1) {
    
    #Format data
    data.model <- data.frame(
      t.start = c(rep(0, n), TR[which(TR<=U)]),
      t.end = c(U[which(TR>U)], TR[which(TR<=U)], U[which(TR<=U)]),
      t.delta = c(delta[which(TR>U)], rep(0, length(which(TR<=U))), delta[which(TR<=U)]),
      V.X = c(X[which(TR>U)], X[which(TR<=U)], X[which(TR<=U)]),
      V.R = c(rep(0, n), rep(1, length(which(TR<=U)))),
      V.XR = c(rep(0, n), X[which(TR<=U)]),
      V.RZ = c(rep(0, n), Z[which(TR<=U)]),
      V.XRZ = c(rep(0, n), X[which(TR<=U)]*Z[which(TR<=U)]),  
      rbind(V[which(TR>U),], V[which(TR<=U),], V[which(TR<=U),])
    )
    
    #Get rid of t.start=t.end=0
    if(length(which(data.model$t.end==0))>0) data.model <- data.model[which(data.model$t.end!=0),]
    
    #Fit the Cox proportional hazard model
    fit <- coxph(Surv(t.start, t.end, t.delta)~., data=data.model)
    
  }  
  
  #Retrieve coefficient and variance estiamtes
  EST <- matrix(as.numeric(fit$coef[1:5]), nrow=5, ncol=1)
  VAR <- fit$v[1:5,1:5]

  #Test for H0: A1B1=A1B2=A2B1=A2B2
  #beta1=beta3=beta4=beta5=0
  D_overall <- matrix(c(1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1), nrow=4, ncol=5, byrow=T)
  test_overall <- t(D_overall %*% EST) %*% solve(D_overall %*% VAR %*% t(D_overall)) %*% (D_overall %*% EST)
  p_overall <- pchisq(test_overall, df=4, lower.tail=F)
 
  #Test for H0: A1=A2
  #beta1=beta3=beta5=0
  D_A12 <- matrix(c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,1), nrow=3, ncol=5, byrow=T)
  test_A12 <- t(D_A12 %*% EST) %*% solve(D_A12 %*% VAR %*% t(D_A12)) %*% (D_A12 %*% EST)
  p_A12 <- pchisq(test_A12, df=3, lower.tail=F)
  
  #Test for H0: B1=B2
  #beta4=beta5=0
  D_B12 <- matrix(c(0,0,0,1,0,0,0,0,0,1), nrow=2, ncol=5, byrow=T)
  test_B12 <- t(D_B12 %*% EST) %*% solve(D_B12 %*% VAR %*% t(D_B12)) %*% (D_B12 %*% EST)
  p_B12 <- pchisq(test_B12, df=2, lower.tail=F)

  #Test for H0: A1B1=A1B2
  #beta4=0
  D_1112 <- matrix(c(0,0,0,1,0), nrow=1, ncol=5, byrow=T)
  test_1112 <- t(D_1112 %*% EST) %*% solve(D_1112 %*% VAR %*% t(D_1112)) %*% (D_1112 %*% EST)
  p_1112 <- pchisq(test_1112, df=1, lower.tail=F)
  
  #Test for H0: A1B1=A2B1
  #beta1=beta3=0
  D_1121 <- matrix(c(1,0,0,0,0,0,0,1,0,0), nrow=2, ncol=5, byrow=T)
  test_1121 <- t(D_1121 %*% EST) %*% solve(D_1121 %*% VAR %*% t(D_1121)) %*% (D_1121 %*% EST)
  p_1121 <- pchisq(test_1121, df=2, lower.tail=F)
  
  #Test for H0: A1B1=A2B2
  #beta1=0 & beta3+beta4+beta5=0
  D_1122 <- matrix(c(1,0,0,0,0,0,0,1,1,1), nrow=2, ncol=5, byrow=T)
  test_1122 <- t(D_1122 %*% EST) %*% solve(D_1122 %*% VAR %*% t(D_1122)) %*% (D_1122 %*% EST)
  p_1122 <- pchisq(test_1122, df=2, lower.tail=F)
  
  #Test for H0: A1B2=A2B1
  #beta1=0 & beta3=beta4
  D_1221 <- matrix(c(1,0,0,0,0,0,0,1,-1,0), nrow=2, ncol=5, byrow=T)
  test_1221 <- t(D_1221 %*% EST) %*% solve(D_1221 %*% VAR %*% t(D_1221)) %*% (D_1221 %*% EST)
  p_1221 <- pchisq(test_1221, df=2, lower.tail=F)
  
  #Test for H0: A1B2=A2B2
  #beta1=0 & beta3+beta5=0
  D_1222 <- matrix(c(1,0,0,0,0,0,0,1,0,1), nrow=2, ncol=5, byrow=T)
  test_1222 <- t(D_1222 %*% EST) %*% solve(D_1222 %*% VAR %*% t(D_1222)) %*% (D_1222 %*% EST)
  p_1222 <- pchisq(test_1222, df=2, lower.tail=F)  

  #Test for H0: A2B1=A2B2
  #beta4+beta5=0
  D_2122 <- matrix(c(0,0,0,1,1), nrow=1, ncol=5, byrow=T)
  test_2122 <- t(D_2122 %*% EST) %*% solve(D_2122 %*% VAR %*% t(D_2122)) %*% (D_2122 %*% EST)
  p_2122 <- pchisq(test_2122, df=1, lower.tail=F) 
    
  #Return results
  TEST <- data.frame(c("A1B1=A1B2=A2B1=A2B2", "A1=A2", "B1=B2", "A1B1=A1B2", "A1B1=A2B1", "A1B1=A2B2", "A1B2=A2B1", "A1B2=A2B2", "A2B1=A2B2"),
                       round(c(test_overall, test_A12, test_B12, test_1112, test_1121, test_1122, test_1221, test_1222, test_2122),4),
                     format.pval(round(c(p_overall, p_A12, p_B12, p_1112, p_1121, p_1122, p_1221, p_1222, p_2122),4),eps=0.0001))
  names(TEST) <- c("H0", "test statistic", "p")  
  return(TEST)    
    
}  
  
