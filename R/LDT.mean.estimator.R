### R code from vignette source 'Lunceford'

###################################################
### code chunk number 1: chunklibraries
###################################################



  #Libraries required
  require(survival)
  require(rms)
 


###################################################
### code chunk number 5: chunkmean
###################################################

LDT.mean.estimator <- function(data, L) {

  n <- nrow(data)
  R <- data$R
  Z <- data$Z
  V <- data$V
  delta <- data$delta
  cens <- 1 - delta
  
  #Chek for errors
  if (is.null(R)) stop("R can not be empty")
  if (is.null(Z)) stop("Z can not be empty")  
  if (is.null(V)) stop("V can not be empty")  
  if (is.null(delta)) stop("delta can not be empty")    
  if (is.character(L)) stop("L has to be numeric")
  if (L<=0) stop("L must be a non-negative value")  
  
  #First step is the calculation on equation 2

  # Obtain probability of Z given R
  
  pi.z <- sum(R*Z) / sum(R)
  
  # Calculate weighting Q
  
	Q1 <-(1-R) + R*(1-Z)/(1-pi.z) 
  Q2 <-(1-R) + R*Z/(pi.z) 
  
  # Kaplan-Meier estimate of the censoring survival curve
  
  cfit <- summary(survfit(Surv(V,cens)~1))
  
	K <- rep(0, n)
	for (i in 1:n) { if(round(V[i],4) < round(min(cfit$time),4)) { K[i] <- 1 } else { dt <- round(cfit$time,4) - round(V[i],4); K[i] <- cfit$surv[which(dt==max(dt[dt<=0]))[1]] } }  
  
  # Calcualte weight
  
  w1 <- w2 <- rep(0, n)  
	w1[which(K!=0)] <- delta[which(K!=0)] * Q1[which(K!=0)] / K[which(K!=0)]
  w2[which(K!=0)] <- delta[which(K!=0)] * Q2[which(K!=0)] / K[which(K!=0)]
    
  # Obtain the survival estimates irrespective of the arms
  s <- rep(0, n)  
  for (i in 1:n) { sind <- as.numeric(V <= V[i]); s[i] <- 1 - sum(delta[which(K!=0)] * sind[which(K!=0)] / K[which(K!=0)]) / sum(delta[which(K!=0)] / K[which(K!=0)]) }
  
  # Calculate the estimate and standard error for each time point
  
	MEAN1 <- MEAN2 <- NA  
  MSE1 <- MSE2 <- NA  

  # Obtain the mean survival estimates

  MEAN1 <- sum(w1*V) / sum(w1)
  MEAN2 <- sum(w2*V) / sum(w2)
    
  # Obtain the variance estimates

  # Define G1' and G2'

	G1 <- G2 <- rep(0, n)

	# Define E1 and E12

	E1 <- E2 <- rep(0, n)

	# Define Y

	Y <- rep(0, n)

	for (k in 1:n) {

		pind <- as.numeric(V >= V[k])

		# Obtain individual G1' and G2'
		
  	if (s[k] != 0) G1[k] <- sum(delta[which(K!=0)] * Q1[which(K!=0)] * (V[which(K!=0)]-MEAN1) * pind[which(K!=0)] / K[which(K!=0)])  / (n * s[k])
		if (s[k] != 0) G2[k] <- sum(delta[which(K!=0)] * Q2[which(K!=0)] * (V[which(K!=0)]-MEAN2) * pind[which(K!=0)] / K[which(K!=0)]) / (n * s[k])		
  	   
	  # Obtain individual E1 and E12

    E1[k] <- sum(delta[K!=0] * (Q1[K!=0]*(V[which(K!=0)]-MEAN1) - G1[k])^2 * pind[which(K!=0)] / K[which(K!=0)]) / n
    E2[k] <- sum(delta[K!=0] * (Q2[K!=0]*(V[which(K!=0)]-MEAN2) - G2[k])^2 * pind[which(K!=0)] / K[which(K!=0)]) / n
        	  
		# Obtain individual Y
	
		Y[k] <- sum(pind)

	}  
    
	# Obtain variance
    
  v1 <- sum(delta[which(K!=0)] * Q1[which(K!=0)]^2 * (V[which(K!=0)]-MEAN1)^2 / K[which(K!=0)]) / (n^2) + sum(E1[which(K!=0)] * (1-delta[which(K!=0)]) * (V[which(K!=0)] <= L) / (K[which(K!=0)] * Y[which(K!=0)])) / n; MSE1 <- sqrt(v1)
	v2 <- sum(delta[which(K!=0)] * Q2[which(K!=0)]^2 * (V[which(K!=0)]-MEAN2)^2 / K[which(K!=0)]) / (n^2) + sum(E2[which(K!=0)] * (1-delta[which(K!=0)]) * (V[which(K!=0)] <= L) / (K[which(K!=0)] * Y[which(K!=0)])) / n; MSE2 <- sqrt(v2)	
		  
  # Combine all the estimates
  
  est <- cbind(MEAN1, MSE1, MEAN2, MSE2)

	return(est)

}






