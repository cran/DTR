###################################################
### Reference:
### Guo X, Tsiatis AA: A weighted risk set estimator for survival distributions in two-stage 
### randomization designs with censored survival data. Int. J. Biostatistics 1:1-15, 2005
###################################################

###################################################
### code chunk number 1: chunklibraries
###################################################

require(survival)

###################################################
### code chunk number 2: chunkWRSE
###################################################

WRSE.estimator <- function(data, # A data frame representing the data for for one of the first-stage assignments from sequentially randomized designs
                                 # data = data frame {TR, R, Z, U, delta}
                           t # a vector with time points of interest 
                             # For example, t=c(1, 3, 5) for the survival estimates at 1, 3, and 5 years respectively
) {
  
  #Retrieve data
  n <- nrow(data)
  TR <- data$TR
  R <- data$R
  Z <- data$Z # Z=0 for B1, 1 for B2
  U <- data$U
  delta <- data$delta

  #Chek for errors
  if (is.null(TR)) stop("TR can not be empty")  
  if (is.null(R)) stop("R can not be empty")
  if (is.null(Z)) stop("Z can not be empty")  
  if (is.null(U)) stop("U can not be empty")  
  if (is.null(delta)) stop("delta can not be empty") 
  
  if (is.null(t)) stop("t can not be empty") 
  for(i in 1:length(t)) {    
    if (is.character(t[i])) stop("t has to be numeric")
    if (t[i]<0) stop("t must be a non-negative value")
  }
  
  #Calculate probability of Z given R
  pi.z <- sum(R*Z) / sum(R)

  #Calculate time-varying weight function for each observed Ui
  #W1 for B1; W2 for B2
  #Time as rows, individual as columns
  w1 <- w2 <- matrix(0, nrow=n, ncol=n)
  #Calculate summation of Wj(Ui)Yj(Ui)
  s1 <- s2 <- rep(0, n)
  #Calculate Yi(Uk) for second part of variance sigma
  #Each individual time as one row
  kind <- matrix(0, nrow=n, ncol=n)
  
  for(i in 1:n) {
    
    #Calculate I(TR <= Ui)
    rind <- as.numeric(TR <= U[i])
    
    #Calculate individual weights for each time Ui
    w1[i,] <- 1 - R*rind + R*rind*(1-Z)/(1-pi.z) # weighting for A1B1
    w2[i,] <- 1 - R*rind + R*rind*Z/pi.z # weighting for A1B2
    
    #Calculate Yj(ui)
    yind <- as.numeric(U >= U[i])
    
    #Caculate summation Wj(Ui)Yj(Ui) for each time Ui
    s1[i] <- sum(w1[i,]*yind)
    s2[i] <- sum(w2[i,]*yind)
    
    #Calculate Yi(Uk) for second part of variance sigma
    kind[i,] <- as.numeric(U[i] >= U)
    
  }
  
  #Calculate the estimate and standard error for each time point  
  SURV1 <- SURV2 <- rep(NA, length(t))  
  SE1 <- SE2 <- COV12 <- rep(NA, length(t))
  
  for(j in 1:length(t)) {
    
    ind <- as.numeric(U <= t[j])
    
    #Calculate the survival estimates    
    SURV1[j] <- exp(-sum(diag(w1)[which(s1!=0)]*delta[which(s1!=0)]*ind[which(s1!=0)]/s1[which(s1!=0)]))
    SURV2[j] <- exp(-sum(diag(w2)[which(s2!=0)]*delta[which(s2!=0)]*ind[which(s2!=0)]/s2[which(s2!=0)]))
    
    #Calculate first part (pI) of the variance sigma
    #Calculate second part (pII) of the variance sigma
    pI1 <- pI2 <- rep(0, n)
    pII1 <- pII2 <- rep(0, n)
    sdp1 <- sdp2 <- sdp12 <- rep(0, n)
    
    for(k in 1:n) {
      
      #Calculate first part for each individual
      if(s1[k] != 0) pI1[k] <- w1[k,k]*delta[k]*ind[k]/s1[k]
      if(s2[k] != 0) pI2[k] <- w2[k,k]*delta[k]*ind[k]/s2[k]
      
      #Calculate second part for each individual
      pII1[k] <- sum(w1[which(s1!=0),k]*kind[k,][which(s1!=0)]*diag(w1)[which(s1!=0)]*delta[which(s1!=0)]*ind[which(s1!=0)]/s1[which(s1!=0)]^2)
      pII2[k] <- sum(w2[which(s2!=0),k]*kind[k,][which(s2!=0)]*diag(w2)[which(s2!=0)]*delta[which(s2!=0)]*ind[which(s2!=0)]/s2[which(s2!=0)]^2)
    
      #Calcualte the square of the difference between first part and second part
      sdp1[k] <- (pI1[k] - pII1[k])^2
      sdp2[k] <- (pI2[k] - pII2[k])^2
      sdp12[k] <- (pI1[k] - pII1[k])*(pI2[k] - pII2[k])
      
    }
    
    #Calculate the variance and covariace sigma
    sig1 <- n*sum(sdp1)
    sig2 <- n*sum(sdp2)
    sig12 <- n*sum(sdp12)
    
    #Calculate the variance and covariace
    v1 <- SURV1[j]^2*sig1/n; SE1[j] <- sqrt(v1)
    v2 <- SURV2[j]^2*sig2/n; SE2[j] <- sqrt(v2)
    COV12[j] <- SURV1[j]*SURV2[j]*sig12/n
       
  }
  
  #Return results  
  TIME <- t  
  est <- data.frame(TIME, SURV1, SURV2, SE1, SE2, COV12)  
  return(est)  
  
} 


