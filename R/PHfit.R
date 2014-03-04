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

PHfit <- function(data # A complete data frame representing the data for two-stage randomization designs
                        # data = data frame {X, TR, R, Z, U, delta, V}
                        # V represents covariates
                        # There could be no covariate, one covariate, and more than one covariates
) {
 
  #Retrieve data
  n <- nrow(data)
  X <- data$X # X=0 for A1, 1 for A2
  TR <- data$TR
  R <- data$R
  Z <- data$Z # Z=0 for B1, 1 for B2
  U <- data$U
  delta <- data$delta
  
  #Chek for errors
  if (is.null(X)) stop("X can not be empty") 
  if (is.null(TR)) stop("TR can not be empty")  
  if (is.null(R)) stop("R can not be empty")
  if (is.null(Z)) stop("Z can not be empty")  
  if (is.null(U)) stop("U can not be empty")  
  if (is.null(delta)) stop("delta can not be empty") 
  
  #Retrive covariates
  V <- data[, ! names(data) %in% c("X", "TR", "R", "Z", "U", "delta")]
  
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
    
    #Change the names
    names(data.model) <- c("start", "end", "delta", "X", "R", "XR", "RZ", "XRZ")
    
    #Fit the Cox proportional hazard model
    fit <- coxph(Surv(start, end, delta) ~ ., data=data.model)
    
    #Change the argument of object
    fit$call <- "coxph(Surv(U, delta)~ X + R + XR + RZ + XRZ)" 
    
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
    
    #Change the names
    names(data.model) <- c("start", "end", "delta", "X", "R", "XR", "RZ", "XRZ", 
                           names(data)[! names(data) %in% c("X", "TR", "R", "Z", "U", "delta")])
    
    #Fit the Cox proportional hazard model
    fit <- coxph(Surv(start, end, delta)~., data=data.model)
    
    #Change the call
    fit$call <- paste("coxph(Surv(U, delta) ~ X + R + XR + RZ + XRZ + ", 
                      names(data)[! names(data) %in% c("X", "TR", "R", "Z", "U", "delta")], ")", sep="")
  
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
    
    #Change the names
    names(data.model) <- c("start", "end", "delta", "X", "R", "XR", "RZ", "XRZ", 
                           names(data)[! names(data) %in% c("X", "TR", "R", "Z", "U", "delta")])
    
    #Fit the Cox proportional hazard model
    fit <- coxph(Surv(start, end, delta)~., data=data.model)
    
    #Change the call
    fit$call <- paste("coxph(Surv(U, delta) ~ X + R + XR + RZ + XRZ + ", 
                      paste(names(data)[! names(data) %in% c("X", "TR", "R", "Z", "U", "delta")], collapse=" + "), ")", sep= "")
    
  }  
    
  #Return class
  class(fit) <- "coxph"
  return(fit)    
    
}  
  
