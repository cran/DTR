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
### code chunk number 2: chunkLDT
###################################################

LDT.estimator <- function(data, # A data frame representing the data for for one of the first-stage assignments from sequentially randomized designs
                                # data = data frame {R, Z, U, delta}
                          t, # a vector with time points of interest. 
                             # For example, t=c(1, 3, 5) for the survival estimates at 1, 3, and 5 years respectively
                          L=.Machine$double.xmax # Optional restricted survival time L
) {

  #Retrieve data
  n <- nrow(data)
  R <- data$R
  Z <- data$Z
	U <- data$U
	delta <- data$delta
  cens <- 1 - delta
  
  #Chek for input errors
  if (is.null(R)) stop("R can not be empty")
  if (is.null(Z)) stop("Z can not be empty")  
  if (is.null(U)) stop("U can not be empty")  
  if (is.null(delta)) stop("delta can not be empty")    

  if (is.null(t)) stop("t can not be empty") 
  for(i in 1:length(t)) {    
    if (is.character(t[i])) stop("t has to be numeric")
    if (t[i]<0) stop("t must be a non-negative value")
  }

  if (is.null(L)) stop("L can not be empty")
  if (is.character(L)) stop("L has to be numeric")
  if (L<=0) stop("L must be a positive value")
  
  #Calculate probability of Z given R  
  pi.z <- sum(R*Z) / sum(R)
  
  #Calculate weighting Q  
	Q1 <-(1-R) + R*(1-Z)/(1-pi.z) # weighting for A1B1
  Q2 <-(1-R) + R*Z/(pi.z) # weighting for A1B2
  
  #Kaplan-Meier estimate of the censoring survival curve (weighting for censoring) 
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
  
  #Calculate the estimate and standard error for each time point  
	SURV1 <- SURV2 <- rep(NA, length(t))  
  SE1 <- SE2 <- COV12 <- rep(NA, length(t))  

	for(j in 1:length(t)) {

    ind <- as.numeric(U <= t[j])
    
    #Calculate the survival estimates
    SURV1[j] <- 1 - sum(w1*ind) / sum(w1)
    SURV2[j] <- 1 - sum(w2*ind) / sum(w2)
    
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
  	  if (s[k] != 0) G1[k] <- sum(delta[which(K!=0)] * Q1[which(K!=0)] * (ind[which(K!=0)]-1+SURV1[j]) * pind[which(K!=0)] / K[which(K!=0)]) / (n * s[k])
		  if (s[k] != 0) G2[k] <- sum(delta[which(K!=0)] * Q2[which(K!=0)] * (ind[which(K!=0)]-1+SURV2[j]) * pind[which(K!=0)] / K[which(K!=0)]) / (n * s[k])		
  	        
		  #Calculate individual E1, E2 and E12
  	  E1[k] <- sum(delta[K!=0] * (Q1[K!=0]*(ind[which(K!=0)]-1+SURV1[j]) - G1[k])^2 * pind[which(K!=0)] / K[which(K!=0)]) / n
		  E2[k] <- sum(delta[K!=0] * (Q2[K!=0]*(ind[which(K!=0)]-1+SURV2[j]) - G2[k])^2 * pind[which(K!=0)] / K[which(K!=0)]) / n
		  E12[k] <- sum(delta[K!=0] * (Q1[K!=0]*(ind[which(K!=0)]-1+SURV1[j]) - G1[k]) * (Q2[K!=0]*(ind[which(K!=0)]-1+SURV2[j]) - G2[k]) * pind[which(K!=0)] / K[which(K!=0)]) / n
		  
		  #Calculate individual Y	
		  Y[k] <- sum(pind)

	  }

	  #Calcualte variance and covariance    
    v1 <- sum(delta[which(K!=0)] * Q1[which(K!=0)]^2 * (ind[which(K!=0)]-1+SURV1[j])^2 / K[which(K!=0)]) / (n^2) + sum(E1[which(K!=0)] * (1-delta[which(K!=0)]) * (U[which(K!=0)] <= L) / (K[which(K!=0)] * Y[which(K!=0)])) / n; SE1[j] <- sqrt(v1)
	  v2 <- sum(delta[which(K!=0)] * Q2[which(K!=0)]^2 * (ind[which(K!=0)]-1+SURV2[j])^2 / K[which(K!=0)]) / (n^2) + sum(E2[which(K!=0)] * (1-delta[which(K!=0)]) * (U[which(K!=0)] <= L) / (K[which(K!=0)] * Y[which(K!=0)])) / n; SE2[j] <- sqrt(v2)		  
    COV12[j] <- sum(delta[which(K!=0)] * Q1[which(K!=0)] * Q2[which(K!=0)] * (ind[which(K!=0)]-1+SURV1[j]) * (ind[which(K!=0)]-1+SURV2[j]) / K[which(K!=0)]) / (n^2) + sum(E12[which(K!=0)] * (1-delta[which(K!=0)]) * (U[which(K!=0)] <= L) / (K[which(K!=0)] * Y[which(K!=0)])) / n
    
	}

  #Return results
  TIME <- t
  est <- data.frame(TIME, SURV1, SURV2, SE1, SE2, COV12)
	return(est)

}

