###################################################
### Reference:
### Lunceford JK, Davidian M, Tsiatis AA: Estimation of survival distributions of treatment 
### policies in two-stage randomization designs in clinical trials. Biometrics 58:48-57, 2002
###################################################

###################################################
### code chunk number 1: chunklibraries
###################################################

#Libraries required
require(survival)

###################################################
### code chunk number 2: chunkmean
###################################################

LDT.mean.estimator <- function(data, # A data frame representing the data for for one of the first-stage assignments from sequentially randomized designs
                                     # data = data frame {R, Z, U, delta}
                               L=.Machine$double.xmax # Optional restricted survival time L
) {

  #Retrieve data
  n <- nrow(data)
  R <- data$R
  Z <- data$Z
  U <- data$U
  delta <- data$delta
  cens <- 1 - delta
  
  #Chek for errors
  if (is.null(R)) stop("R can not be empty")
  if (is.null(Z)) stop("Z can not be empty")  
  if (is.null(U)) stop("U can not be empty")  
  if (is.null(delta)) stop("delta can not be empty")
  
  if (is.null(L)) stop("L can not be empty")
  if (is.character(L)) stop("L has to be numeric")
  if (L<=0) stop("L must be a positive value")  
  
  #Calculate probability of Z given R  
  pi.z <- sum(R*Z) / sum(R)
  
  #Calculate weighting Q  
	Q1 <-(1-R) + R*(1-Z)/(1-pi.z) # weighting for A1B1
  Q2 <-(1-R) + R*Z/(pi.z) # weighting for A1B2
  
  #Kaplan-Meier estimate of the censoring survival curve  
  cfit <- summary(survfit(Surv(U,cens)~1))  
	K <- rep(0, n)
	for (i in 1:n) { if(round(U[i],4) < round(min(cfit$time),4)) { K[i] <- 1 } else { dt <- round(cfit$time,4) - round(U[i],4); K[i] <- cfit$surv[which(dt==max(dt[dt<=0]))[1]] } }  
  
  #Combine two forms of inverse probability weighting  
  w1 <- w2 <- rep(0, n)  
	w1[which(K!=0)] <- delta[which(K!=0)] * Q1[which(K!=0)] / K[which(K!=0)] # weighting for A1B1
  w2[which(K!=0)] <- delta[which(K!=0)] * Q2[which(K!=0)] / K[which(K!=0)] # weighting for A1B2
    
  #Calculate the survival estimates irrespective of the arms
  s <- rep(0, n)  
  for (i in 1:n) { sind <- as.numeric(U <= U[i]); s[i] <- 1 - sum(delta[which(K!=0)] * sind[which(K!=0)] / K[which(K!=0)]) / sum(delta[which(K!=0)] / K[which(K!=0)]) }
  
  #Calculate the estimate and standard error for restricted mean survival  
	MEAN1 <- MEAN2 <- NA  
  MSE1 <- MSE2 <- MCOV12 <- NA  

  #Calculate the restricted mean survival estimates
  MEAN1 <- sum(w1*U) / sum(w1)
  MEAN2 <- sum(w2*U) / sum(w2)
    
  #Calculate the variance estimates
  #Define G1' and G2'
	G1 <- G2 <- rep(0, n)
	#Define E1 and E12
	E1 <- E2 <- E12 <- rep(0, n)
	#Define Y
	Y <- rep(0, n)

	for (k in 1:n) {

		pind <- as.numeric(U >= U[k])

		#Calculate individual G1' and G2'		
  	if (s[k] != 0) G1[k] <- sum(delta[which(K!=0)] * Q1[which(K!=0)] * (U[which(K!=0)]-MEAN1) * pind[which(K!=0)] / K[which(K!=0)])  / (n * s[k])
		if (s[k] != 0) G2[k] <- sum(delta[which(K!=0)] * Q2[which(K!=0)] * (U[which(K!=0)]-MEAN2) * pind[which(K!=0)] / K[which(K!=0)]) / (n * s[k])		
  	   
	  #Calculate individual E1 and E12
    E1[k] <- sum(delta[K!=0] * (Q1[K!=0]*(U[which(K!=0)]-MEAN1) - G1[k])^2 * pind[which(K!=0)] / K[which(K!=0)]) / n
    E2[k] <- sum(delta[K!=0] * (Q2[K!=0]*(U[which(K!=0)]-MEAN2) - G2[k])^2 * pind[which(K!=0)] / K[which(K!=0)]) / n
		E12[k] <- sum(delta[K!=0] * (Q1[K!=0]*(U[which(K!=0)]-MEAN1) - G1[k]) * (Q2[K!=0]*(U[which(K!=0)]-MEAN2) - G2[k]) * pind[which(K!=0)] / K[which(K!=0)]) / n
		
		#Calculate individual Y	
		Y[k] <- sum(pind)

	}  
    
	#Calculate variance and covariance    
  v1 <- sum(delta[which(K!=0)] * Q1[which(K!=0)]^2 * (U[which(K!=0)]-MEAN1)^2 / K[which(K!=0)]) / (n^2) + sum(E1[which(K!=0)] * (1-delta[which(K!=0)]) * (U[which(K!=0)] <= L) / (K[which(K!=0)] * Y[which(K!=0)])) / n; MSE1 <- sqrt(v1)
	v2 <- sum(delta[which(K!=0)] * Q2[which(K!=0)]^2 * (U[which(K!=0)]-MEAN2)^2 / K[which(K!=0)]) / (n^2) + sum(E2[which(K!=0)] * (1-delta[which(K!=0)]) * (U[which(K!=0)] <= L) / (K[which(K!=0)] * Y[which(K!=0)])) / n; MSE2 <- sqrt(v2)	
  MCOV12 <- sum(delta[which(K!=0)] * Q1[which(K!=0)] * Q2[which(K!=0)] * (U[which(K!=0)]-MEAN1) * (U[which(K!=0)]-MEAN2) / K[which(K!=0)]) / (n^2) + sum(E12[which(K!=0)] * (1-delta[which(K!=0)]) * (U[which(K!=0)] <= L) / (K[which(K!=0)] * Y[which(K!=0)])) / n
  
  #Return results  
  est <- data.frame(MEAN1, MEAN2, MSE1, MSE2, MCOV12)
	return(est)

}






