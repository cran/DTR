###################################################
### Reference:
### Tang X, Wahed AS: Cumulative hazard ratio estimation for treatment regimes in
### sequentially randomized clinical trials. Statistics in Biosciences, [Epub ahead of print]
###################################################

###################################################
### code chunk number 1: chunklibraries
###################################################

#Libraries required
require(survival)

###################################################
### code chunk number 2: chunkCHR
###################################################

CHR.estimator <- function(fdata, # A complete data frame representing the data for two-stage randomization designs
                                 # fdata = data frame {X, R, Z, U, delta, V}
                                 # V represents covariates
                                 # There could be one covariate, or more than one covariates
                                 # The function does not allow the absence of covariates
                          t # a vector with time points of interest 
                            # For example, t=c(1, 3, 5) for the survival estimates at 1, 3, and 5 years respectively
) {
  
  #Define results
  est <- NULL
  
  #Retrieve data
  n <- nrow(fdata)
  X <- fdata$X # X=0 for A1, X=1 for A2
  R <- fdata$R
  Z <- fdata$Z # Z=0 for B1, Z=1 for B2
  U <- fdata$U
  delta <- fdata$delta
  
  #Chek for errors
  if (is.null(X)) stop("R can not be empty")
  if (is.null(R)) stop("R can not be empty")
  if (is.null(Z)) stop("Z can not be empty")  
  if (is.null(U)) stop("V can not be empty")  
  if (is.null(delta)) stop("delta can not be empty") 
  
  if (is.null(t)) stop("t can not be empty") 
  for(i in 1:length(t)) {
    if (is.character(t[i])) stop("t has to be numeric")
    if (t[i]<0) stop("t must be a non-negative value")      
  }
  
  #Estimate probability of being assigned to A2
  pi.x <- sum(X)/n
  
  #Estimate probability of being assigned to B2 (allowing probability to vary across A1 and A2)
  pi.z1 <- sum((1-X)*R*Z) / sum((1-X)*R)
  pi.z2 <- sum(X*R*Z) / sum(X*R)
  
  #Calculate weight for A1B1, A1B2, A2B1, and A2B2
  w11 <- (1-X)*(1-R)/(1-pi.x) + (1-X)*R*(1-Z)/((1-pi.x)*(1-pi.z1))
  w12 <- (1-X)*(1-R)/(1-pi.x) + (1-X)*R*Z/((1-pi.x)*pi.z1)
  w21 <- X*(1-R)/pi.x + X*R*(1-Z)/(pi.x*(1-pi.z1))
  w22 <- X*(1-R)/pi.x + X*R*Z/(pi.x*pi.z1)
  
  #Retrive covariates
  V <- as.matrix(fdata[, ! names(fdata) %in% c("X", "R", "Z", "U", "delta")])
  
  #################################################
  #################################################
  #If no covariates: ERROR
  #################################################
  #################################################
  
  if(NCOL(V)==0) stop("Covariates can not be empty")  
  
  #################################################
  #################################################
  #If only 1 covariates
  #################################################
  #################################################
  if(NCOL(V)==1) { 
    
    #Obtain the inital beta estimates
    beta <- as.numeric(coxph(Surv(U, delta)~., data=fdata[, ! names(fdata) %in% c("X", "R", "Z")])$coef)
    
    #Solve for beta using Newton-Raphson method
    print("Calling for updateBeta() function to solve for coefficients...")

    for (p in 1:1000) {
 
      #Run updateBeta function
      ebeta <- updateBeta(beta, V, U, delta, w11, w12, w21, w22)
    
      #Calculate difference between updated ebeta and beta
      index <- max(abs(ebeta-beta))
    
      if (index <= 10^(-6)) break else { beta <- ebeta; p <- p + 1 }
    
    }
    
    print("Calculating cumulative hazard ratios and variance/covariance...")
    
    #Change V from matrix to numeric
    V <- as.numeric(V)
    
    #Calcualte exp(beta*V)
    e <- exp(ebeta*V)
    
    for(j in 1:length(t)) {
      
      #Define s0 and v1_bar
      s00 <- s01 <- s02 <- s03 <- rep(0, n)
      v10_bar <- v11_bar <- v12_bar <- v13_bar <- rep(0, n)
      
      #Define lambda and omega
      lambda11 <- lambda12 <- lambda21 <- lambda22 <- 0
      omega <- 0
      
      #Define h
      h11 <- h12 <- h21 <- h22 <- 0
      
      for(i in 1:n) {
        
        ind <- as.numeric(U >= U[i])
        
        #Calculate s0j      
        s00[i] <- sum(ind * w11 * e) / n
        s01[i] <- sum(ind * w12 * e) / n
        s02[i] <- sum(ind * w21 * e) / n
        s03[i] <- sum(ind * w22 * e) / n
        
        #Calculate s1j   
        s10 <- sum(ind * w11 * e * V) / n
        s11 <- sum(ind * w12 * e * V) / n
        s12 <- sum(ind * w21 * e * V) / n
        s13 <- sum(ind * w22 * e * V) / n
        
        #Calculate s2j
        s20 <- sum(ind * w11 * e * V * V) / n
        s21 <- sum(ind * w12 * e * V * V) / n
        s22 <- sum(ind * w21 * e * V * V) / n
        s23 <- sum(ind * w22 * e * V * V) / n
        
        #Calculate z1?_bar
        if(s00[i] != 0) v10_bar[i] <- s10/s00[i] else v10_bar[i] <- 0
        if(s01[i] != 0) v11_bar[i] <- s11/s01[i] else v11_bar[i] <- 0
        if(s02[i] != 0) v12_bar[i] <- s12/s02[i] else v12_bar[i] <- 0
        if(s03[i] != 0) v13_bar[i] <- s13/s03[i] else v13_bar[i] <- 0
        
        #calculate z2?_bar
        if(s00[i] != 0) v20_bar <- s20/s00[i] else v20_bar <- 0
        if(s01[i] != 0) v21_bar <- s21/s01[i] else v21_bar <- 0
        if(s02[i] != 0) v22_bar <- s22/s02[i] else v22_bar <- 0
        if(s03[i] != 0) v23_bar <- s23/s03[i] else v23_bar <- 0     
        
        #Calculate each cumulative baseline hazards
        if(s00[i] != 0) lambda11 <- lambda11 + delta[i] * w11[i] * as.numeric(U[i] <= t[j]) / (n*s00[i])
        if(s01[i] != 0) lambda12 <- lambda12 + delta[i] * w12[i] * as.numeric(U[i] <= t[j]) / (n*s01[i])
        if(s02[i] != 0) lambda21 <- lambda21 + delta[i] * w21[i] * as.numeric(U[i] <= t[j]) / (n*s02[i])
        if(s03[i] != 0) lambda22 <- lambda22 + delta[i] * w22[i] * as.numeric(U[i] <= t[j]) / (n*s03[i])
        
        #Calculate each h
        if(s00[i] != 0) h11 <- h11 - delta[i]*w11[i]*as.numeric(U[i] <= t[j])*v10_bar[i] / (n*s00[i])
        if(s01[i] != 0) h12 <- h12 - delta[i]*w12[i]*as.numeric(U[i] <= t[j])*v11_bar[i] / (n*s01[i])
        if(s02[i] != 0) h21 <- h21 - delta[i]*w21[i]*as.numeric(U[i] <= t[j])*v12_bar[i] / (n*s02[i])
        if(s03[i] != 0) h22 <- h22 - delta[i]*w22[i]*as.numeric(U[i] <= t[j])*v13_bar[i] / (n*s03[i])
        
        #Calculate each tao = s2/s0 - v1^2
        tao11i <- v20_bar - v10_bar[i] * v10_bar[i]
        tao12i <- v21_bar - v11_bar[i] * v11_bar[i]
        tao21i <- v22_bar - v12_bar[i] * v12_bar[i]
        tao22i <- v23_bar - v13_bar[i] * v13_bar[i]
        
        #Calculate each omega
        omega <- omega + (delta[i]*w11[i]*tao11i + delta[i]*w12[i]*tao12i + delta[i]*w21[i]*tao21i + delta[i]*w22[i]*tao22i) / n
        
      }    
      
      #Calculate cumulative hazards ratio
      CHR1211 <- lambda12 / lambda11; LogCHR1211 <- log(CHR1211)
      CHR2111 <- lambda21 / lambda11; LogCHR2111 <- log(CHR2111)
      CHR2211 <- lambda22 / lambda11; LogCHR2211 <- log(CHR2211)
      CHR2112 <- lambda21 / lambda12; LogCHR2112 <- log(CHR2112)
      CHR2212 <- lambda22 / lambda12; LogCHR2212 <- log(CHR2212)
      CHR2221 <- lambda22 / lambda21; LogCHR2221 <- log(CHR2221)
      
      # Define xi 
      xi1211 <- xi2111 <- xi2211 <- xi2112 <- xi2212 <- xi2221 <- rep(0, n)
      
      for(k in 1:n) {
        
        #Defind I(Uk >= U)
        uind <- as.numeric(U[k] >= U)
        
        #Calculate each latter part of psi
        psi11i <- w11[k]*uind*e[k]*delta*w11*(V[k]-v10_bar) / (n*s00); psi11i[is.na(psi11i)] <- 0
        psi12i <- w12[k]*uind*e[k]*delta*w12*(V[k]-v11_bar) / (n*s01); psi12i[is.na(psi12i)] <- 0
        psi21i <- w21[k]*uind*e[k]*delta*w21*(V[k]-v12_bar) / (n*s02); psi21i[is.na(psi21i)] <- 0
        psi22i <- w22[k]*uind*e[k]*delta*w22*(V[k]-v13_bar) / (n*s03); psi22i[is.na(psi22i)] <- 0
        
        #Calculate psi
        psi <- delta[k]*w11[k]*(V[k] - v10_bar[k]) + 
          delta[k]*w12[k]*(V[k] - v11_bar[k]) +
          delta[k]*w21[k]*(V[k] - v12_bar[k]) +
          delta[k]*w22[k]*(V[k] - v13_bar[k]) -
          sum(psi11i) - sum(psi12i) - sum(psi21i) - sum(psi22i)
        
        #Calculate each latter part of phi integral         
        intl11i <- w11[k]*uind*e[k]*delta*w11*as.numeric(U<=t[j]) / (n*s00*s00); intl11i[is.na(intl11i)] <- 0
        intl12i <- w12[k]*uind*e[k]*delta*w12*as.numeric(U<=t[j]) / (n*s01*s01); intl12i[is.na(intl12i)] <- 0
        intl21i <- w21[k]*uind*e[k]*delta*w21*as.numeric(U<=t[j]) / (n*s02*s02); intl21i[is.na(intl21i)] <- 0
        intl22i <- w22[k]*uind*e[k]*delta*w22*as.numeric(U<=t[j]) / (n*s03*s03); intl22i[is.na(intl22i)] <- 0
        
        #Calculate each phi_L      
        if(s00[k] != 0) phiL11i <- delta[k]*w11[k]*as.numeric(U[k]<=t[j])/s00[k] - sum(intl11i) else phiL11i <- 0
        if(s01[k] != 0) phiL12i <- delta[k]*w12[k]*as.numeric(U[k]<=t[j])/s01[k] - sum(intl12i) else phiL12i <- 0
        if(s02[k] != 0) phiL21i <- delta[k]*w21[k]*as.numeric(U[k]<=t[j])/s02[k] - sum(intl21i) else phiL21i <- 0
        if(s03[k] != 0) phiL22i <- delta[k]*w22[k]*as.numeric(U[k]<=t[j])/s03[k] - sum(intl22i) else phiL22i <- 0
        
        #Calculate each phi      
        phi11 <- h11 * psi / omega + phiL11i
        phi12 <- h12 * psi / omega + phiL12i
        phi21 <- h21 * psi / omega + phiL21i
        phi22 <- h22 * psi / omega + phiL22i
        
        #Calculate individual xi
        
        xi1211[k] <- phi12 / lambda11 - lambda12 * phi11 / (lambda11)^2
        xi2111[k] <- phi21 / lambda11 - lambda21 * phi11 / (lambda11)^2
        xi2211[k] <- phi22 / lambda11 - lambda22 * phi11 / (lambda11)^2
        xi2112[k] <- phi21 / lambda12 - lambda21 * phi12 / (lambda12)^2
        xi2212[k] <- phi22 / lambda12 - lambda22 * phi12 / (lambda12)^2
        xi2221[k] <- phi22 / lambda21 - lambda22 * phi21 / (lambda21)^2
        
      }
      
      #Save the results
      
      est[[j]] <- list(ebeta, t[j],
                       data.frame(CHR1211, CHR2111, CHR2211, CHR2112, CHR2212, CHR2221),
                       matrix(c(mean(xi1211*xi1211)/n, mean(xi1211*xi2111)/n, mean(xi1211*xi2211)/n,
                                mean(xi1211*xi2112)/n, mean(xi1211*xi2212)/n, mean(xi1211*xi2221)/n,
                                mean(xi2111*xi1211)/n, mean(xi2111*xi2111)/n, mean(xi2111*xi2211)/n,
                                mean(xi2111*xi2112)/n, mean(xi2111*xi2212)/n, mean(xi2111*xi2221)/n,
                                mean(xi2211*xi1211)/n, mean(xi2211*xi2111)/n, mean(xi2211*xi2211)/n,
                                mean(xi2211*xi2112)/n, mean(xi2211*xi2212)/n, mean(xi2211*xi2221)/n,
                                mean(xi2112*xi1211)/n, mean(xi2112*xi2111)/n, mean(xi2112*xi2211)/n,
                                mean(xi2112*xi2112)/n, mean(xi2112*xi2212)/n, mean(xi2112*xi2221)/n,
                                mean(xi2212*xi1211)/n, mean(xi2212*xi2111)/n, mean(xi2212*xi2211)/n,
                                mean(xi2212*xi2112)/n, mean(xi2212*xi2212)/n, mean(xi2212*xi2221)/n,
                                mean(xi2221*xi1211)/n, mean(xi2221*xi2111)/n, mean(xi2221*xi2211)/n,
                                mean(xi2221*xi2112)/n, mean(xi2221*xi2212)/n, mean(xi2221*xi2221)/n), nrow=6, ncol=6),
                       data.frame(LogCHR1211, LogCHR2111, LogCHR2211, LogCHR2112, LogCHR2212, LogCHR2221),
                       matrix(c(mean(xi1211*xi1211)/(n*CHR1211*CHR1211), mean(xi1211*xi2111)/(n*CHR1211*CHR2111), 
                                mean(xi1211*xi2211)/(n*CHR1211*CHR2211), mean(xi1211*xi2112)/(n*CHR1211*CHR2112), 
                                mean(xi1211*xi2212)/(n*CHR1211*CHR2212), mean(xi1211*xi2221)/(n*CHR1211*CHR2221),
                                mean(xi2111*xi1211)/(n*CHR2111*CHR1211), mean(xi2111*xi2111)/(n*CHR2111*CHR2111), 
                                mean(xi2111*xi2211)/(n*CHR2111*CHR2211), mean(xi2111*xi2112)/(n*CHR2111*CHR2112), 
                                mean(xi2111*xi2212)/(n*CHR2111*CHR2212), mean(xi2111*xi2221)/(n*CHR2111*CHR2221),
                                mean(xi2211*xi1211)/(n*CHR2211*CHR1211), mean(xi2211*xi2111)/(n*CHR2211*CHR2111), 
                                mean(xi2211*xi2211)/(n*CHR2211*CHR2211), mean(xi2211*xi2112)/(n*CHR2211*CHR2112), 
                                mean(xi2211*xi2212)/(n*CHR2211*CHR2212), mean(xi2211*xi2221)/(n*CHR2211*CHR2221),
                                mean(xi2112*xi1211)/(n*CHR2112*CHR1211), mean(xi2112*xi2111)/(n*CHR2112*CHR2111), 
                                mean(xi2112*xi2211)/(n*CHR2112*CHR2211), mean(xi2112*xi2112)/(n*CHR2112*CHR2112), 
                                mean(xi2112*xi2212)/(n*CHR2112*CHR2212), mean(xi2112*xi2221)/(n*CHR2112*CHR2221),
                                mean(xi2212*xi1211)/(n*CHR2212*CHR1211), mean(xi2212*xi2111)/(n*CHR2212*CHR2111), 
                                mean(xi2212*xi2211)/(n*CHR2212*CHR2211), mean(xi2212*xi2112)/(n*CHR2212*CHR2112), 
                                mean(xi2212*xi2212)/(n*CHR2212*CHR2212), mean(xi2212*xi2221)/(n*CHR2212*CHR2221),
                                mean(xi2221*xi1211)/(n*CHR2221*CHR1211), mean(xi2221*xi2111)/(n*CHR2221*CHR2111), 
                                mean(xi2221*xi2211)/(n*CHR2221*CHR2211), mean(xi2221*xi2112)/(n*CHR2221*CHR2112), 
                                mean(xi2221*xi2212)/(n*CHR2221*CHR2212), mean(xi2221*xi2221)/(n*CHR2221*CHR2221)), nrow=6, ncol=6))
      names(est[[j]]) <- c("COEF", "TIME", "CHR", "VAR[CHR]", "LOG(CHR)", "VAR[LOG(CHR)]")
      
    }
    
  }
  
  #################################################
  #################################################
  #If more than 1 covariates
  #################################################
  #################################################
  if(NCOL(V)>1) { 
    
    #Obtain the inital beta estimates
    beta <- as.numeric(coxph(Surv(U, delta)~., data=fdata[, ! names(fdata) %in% c("X", "R", "Z")])$coef)
    
    #Solve for beta using Newton-Raphson method
    print("Calling for updateBeta() function to solve for coefficients...")

    for (p in 1:1000) {
      
      #Run updateBeta function
      ebeta <- updateBeta(beta, V, U, delta, w11, w12, w21, w22)
      
      #Calculate difference between updated ebeta and beta
      index <- max(abs(ebeta-beta))
      
      if (index <= 10^(-6)) break else { beta <- ebeta; p <- p + 1 }
      
    }
    
    print("Calculating cumulative hazard ratios and variance/covariance...")
    
    #Calcualte exp(beta*V)
    e <- as.numeric(exp(matrix(ebeta, nrow=1, ncol=NCOL(V)) %*% t(V)))
    
    for(j in 1:length(t)) {
    
      #Define s0 and v1_bar
      s00 <- s01 <- s02 <- s03 <- rep(0, n)
      v10_bar <- v11_bar <- v12_bar <- v13_bar <- matrix(0, nrow=NCOL(V), ncol=n)
    
      #Define lambda and omega
      lambda11 <- lambda12 <- lambda21 <- lambda22 <- 0
      omega <- matrix(0, nrow=NCOL(V), ncol=NCOL(V))
    
      #Define h
      h11 <- h12 <- h21 <- h22 <- matrix(0, nrow=NCOL(V), ncol=1)
    
      for(i in 1:n) {
      
        ind <- as.numeric(U >= U[i])
      
        #Calculate s0j      
        s00[i] <- sum(ind * w11 * e) / n
        s01[i] <- sum(ind * w12 * e) / n
        s02[i] <- sum(ind * w21 * e) / n
        s03[i] <- sum(ind * w22 * e) / n
      
        #Calculate s1j   
        s10 <- t(matrix(ind*w11*e,nrow=1, ncol=n) %*% V / n)
        s11 <- t(matrix(ind*w12*e,nrow=1, ncol=n) %*% V / n)
        s12 <- t(matrix(ind*w21*e,nrow=1, ncol=n) %*% V / n)
        s13 <- t(matrix(ind*w22*e,nrow=1, ncol=n) %*% V / n)
      
        #Calculate s2j
        s20 <- t(V) %*% diag(as.numeric(ind*w11*e)) %*% V / n
        s21 <- t(V) %*% diag(as.numeric(ind*w12*e)) %*% V / n
        s22 <- t(V) %*% diag(as.numeric(ind*w21*e)) %*% V / n
        s23 <- t(V) %*% diag(as.numeric(ind*w22*e)) %*% V / n
      
        #Calculate z1?_bar
        if(s00[i] != 0) v10_bar[,i] <- s10/s00[i] else v10_bar[,i] <- matrix(0, nrow=NCOL(V), ncol=1)
        if(s01[i] != 0) v11_bar[,i] <- s11/s01[i] else v11_bar[,i] <- matrix(0, nrow=NCOL(V), ncol=1)
        if(s02[i] != 0) v12_bar[,i] <- s12/s02[i] else v12_bar[,i] <- matrix(0, nrow=NCOL(V), ncol=1)
        if(s03[i] != 0) v13_bar[,i] <- s13/s03[i] else v13_bar[,i] <- matrix(0, nrow=NCOL(V), ncol=1)
      
        #calculate z2?_bar
        if(s00[i] != 0) v20_bar <- s20/s00[i] else v20_bar <- matrix(0, nrow=NCOL(V), ncol=NCOL(V))
        if(s01[i] != 0) v21_bar <- s21/s01[i] else v21_bar <- matrix(0, nrow=NCOL(V), ncol=NCOL(V))
        if(s02[i] != 0) v22_bar <- s22/s02[i] else v22_bar <- matrix(0, nrow=NCOL(V), ncol=NCOL(V))
        if(s03[i] != 0) v23_bar <- s23/s03[i] else v23_bar <- matrix(0, nrow=NCOL(V), ncol=NCOL(V))      
      
        #Calculate each cumulative baseline hazards
        if(s00[i] != 0) lambda11 <- lambda11 + delta[i] * w11[i] * as.numeric(U[i] <= t[j]) / (n*s00[i])
        if(s01[i] != 0) lambda12 <- lambda12 + delta[i] * w12[i] * as.numeric(U[i] <= t[j]) / (n*s01[i])
        if(s02[i] != 0) lambda21 <- lambda21 + delta[i] * w21[i] * as.numeric(U[i] <= t[j]) / (n*s02[i])
        if(s03[i] != 0) lambda22 <- lambda22 + delta[i] * w22[i] * as.numeric(U[i] <= t[j]) / (n*s03[i])
      
        #Calculate each h
        if(s00[i] != 0) h11 <- h11 - delta[i]*w11[i]*as.numeric(U[i] <= t[j])*v10_bar[,i] / (n*s00[i])
        if(s01[i] != 0) h12 <- h12 - delta[i]*w12[i]*as.numeric(U[i] <= t[j])*v11_bar[,i] / (n*s01[i])
        if(s02[i] != 0) h21 <- h21 - delta[i]*w21[i]*as.numeric(U[i] <= t[j])*v12_bar[,i] / (n*s02[i])
        if(s03[i] != 0) h22 <- h22 - delta[i]*w22[i]*as.numeric(U[i] <= t[j])*v13_bar[,i] / (n*s03[i])
      
        #Calculate each tao = s2/s0 - v1^2
        tao11i <- v20_bar - v10_bar[,i] %*% t(v10_bar[,i])
        tao12i <- v21_bar - v11_bar[,i] %*% t(v11_bar[,i])
        tao21i <- v22_bar - v12_bar[,i] %*% t(v12_bar[,i])
        tao22i <- v23_bar - v13_bar[,i] %*% t(v13_bar[,i])
        
        #Calculate each omega
        omega <- omega + (delta[i]*w11[i]*tao11i + delta[i]*w12[i]*tao12i + delta[i]*w21[i]*tao21i + delta[i]*w22[i]*tao22i) / n
      
      }    
  
      #Calculate cumulative hazards ratio
      CHR1211 <- lambda12 / lambda11; LogCHR1211 <- log(CHR1211)
      CHR2111 <- lambda21 / lambda11; LogCHR2111 <- log(CHR2111)
      CHR2211 <- lambda22 / lambda11; LogCHR2211 <- log(CHR2211)
      CHR2112 <- lambda21 / lambda12; LogCHR2112 <- log(CHR2112)
      CHR2212 <- lambda22 / lambda12; LogCHR2212 <- log(CHR2212)
      CHR2221 <- lambda22 / lambda21; LogCHR2221 <- log(CHR2221)
    
      # Define xi
      xi1211 <- xi2111 <- xi2211 <- xi2112 <- xi2212 <- xi2221 <- rep(0, n)
    
      for(k in 1:n) {
      
        #Defind I(Uk >= U)
        uind <- as.numeric(U[k] >= U)
      
        #Calculate each latter part of psi
        psi11i <- w11[k]*uind*e[k]*delta*w11*(V[k,]-v10_bar) / (n*s00); psi11i[is.na(psi11i)] <- 0
        psi12i <- w12[k]*uind*e[k]*delta*w12*(V[k,]-v11_bar) / (n*s01); psi12i[is.na(psi12i)] <- 0
        psi21i <- w21[k]*uind*e[k]*delta*w21*(V[k,]-v12_bar) / (n*s02); psi21i[is.na(psi21i)] <- 0
        psi22i <- w22[k]*uind*e[k]*delta*w22*(V[k,]-v13_bar) / (n*s03); psi22i[is.na(psi22i)] <- 0

        #Calculate psi
        psi <- delta[k]*w11[k]*(matrix(V[k,], nrow=NCOL(V), ncol=1) - v10_bar[,k]) + 
          delta[k]*w12[k]*(matrix(V[k,], nrow=NCOL(V), ncol=1) - v11_bar[,k]) +
          delta[k]*w21[k]*(matrix(V[k,], nrow=NCOL(V), ncol=1) - v12_bar[,k]) +
          delta[k]*w22[k]*(matrix(V[k,], nrow=NCOL(V), ncol=1) - v13_bar[,k]) -
          sum(psi11i) - sum(psi12i) - sum(psi21i) - sum(psi22i)
            
        #Calculate each latter part of phi integral         
        intl11i <- w11[k]*uind*e[k]*delta*w11*as.numeric(U<=t[j]) / (n*s00*s00); intl11i[is.na(intl11i)] <- 0
        intl12i <- w12[k]*uind*e[k]*delta*w12*as.numeric(U<=t[j]) / (n*s01*s01); intl12i[is.na(intl12i)] <- 0
        intl21i <- w21[k]*uind*e[k]*delta*w21*as.numeric(U<=t[j]) / (n*s02*s02); intl21i[is.na(intl21i)] <- 0
        intl22i <- w22[k]*uind*e[k]*delta*w22*as.numeric(U<=t[j]) / (n*s03*s03); intl22i[is.na(intl22i)] <- 0

        #Calculate each phi_L      
        if(s00[k] != 0) phiL11i <- delta[k]*w11[k]*as.numeric(U[k]<=t[j])/s00[k] - sum(intl11i) else phiL11i <- 0
        if(s01[k] != 0) phiL12i <- delta[k]*w12[k]*as.numeric(U[k]<=t[j])/s01[k] - sum(intl12i) else phiL12i <- 0
        if(s02[k] != 0) phiL21i <- delta[k]*w21[k]*as.numeric(U[k]<=t[j])/s02[k] - sum(intl21i) else phiL21i <- 0
        if(s03[k] != 0) phiL22i <- delta[k]*w22[k]*as.numeric(U[k]<=t[j])/s03[k] - sum(intl22i) else phiL22i <- 0
      
        #Calculate each phi      
        phi11 <- t(h11) %*% solve(omega) %*% psi + phiL11i
        phi12 <- t(h12) %*% solve(omega) %*% psi + phiL12i
        phi21 <- t(h21) %*% solve(omega) %*% psi + phiL21i
        phi22 <- t(h22) %*% solve(omega) %*% psi + phiL22i
      
        #Calculate individual xi
      
        xi1211[k] <- phi12 / lambda11 - lambda12 * phi11 / (lambda11)^2
        xi2111[k] <- phi21 / lambda11 - lambda21 * phi11 / (lambda11)^2
        xi2211[k] <- phi22 / lambda11 - lambda22 * phi11 / (lambda11)^2
        xi2112[k] <- phi21 / lambda12 - lambda21 * phi12 / (lambda12)^2
        xi2212[k] <- phi22 / lambda12 - lambda22 * phi12 / (lambda12)^2
        xi2221[k] <- phi22 / lambda21 - lambda22 * phi21 / (lambda21)^2
      
      }

      #Save the results
    
      est[[j]] <- list(ebeta, t[j],
                        data.frame(CHR1211, CHR2111, CHR2211, CHR2112, CHR2212, CHR2221),
                        matrix(c(mean(xi1211*xi1211)/n, mean(xi1211*xi2111)/n, mean(xi1211*xi2211)/n,
                             mean(xi1211*xi2112)/n, mean(xi1211*xi2212)/n, mean(xi1211*xi2221)/n,
                             mean(xi2111*xi1211)/n, mean(xi2111*xi2111)/n, mean(xi2111*xi2211)/n,
                             mean(xi2111*xi2112)/n, mean(xi2111*xi2212)/n, mean(xi2111*xi2221)/n,
                             mean(xi2211*xi1211)/n, mean(xi2211*xi2111)/n, mean(xi2211*xi2211)/n,
                             mean(xi2211*xi2112)/n, mean(xi2211*xi2212)/n, mean(xi2211*xi2221)/n,
                             mean(xi2112*xi1211)/n, mean(xi2112*xi2111)/n, mean(xi2112*xi2211)/n,
                             mean(xi2112*xi2112)/n, mean(xi2112*xi2212)/n, mean(xi2112*xi2221)/n,
                             mean(xi2212*xi1211)/n, mean(xi2212*xi2111)/n, mean(xi2212*xi2211)/n,
                             mean(xi2212*xi2112)/n, mean(xi2212*xi2212)/n, mean(xi2212*xi2221)/n,
                             mean(xi2221*xi1211)/n, mean(xi2221*xi2111)/n, mean(xi2221*xi2211)/n,
                             mean(xi2221*xi2112)/n, mean(xi2221*xi2212)/n, mean(xi2221*xi2221)/n), nrow=6, ncol=6),
                        data.frame(LogCHR1211, LogCHR2111, LogCHR2211, LogCHR2112, LogCHR2212, LogCHR2221),
                        matrix(c(mean(xi1211*xi1211)/(n*CHR1211*CHR1211), mean(xi1211*xi2111)/(n*CHR1211*CHR2111), 
                                mean(xi1211*xi2211)/(n*CHR1211*CHR2211), mean(xi1211*xi2112)/(n*CHR1211*CHR2112), 
                                mean(xi1211*xi2212)/(n*CHR1211*CHR2212), mean(xi1211*xi2221)/(n*CHR1211*CHR2221),
                                mean(xi2111*xi1211)/(n*CHR2111*CHR1211), mean(xi2111*xi2111)/(n*CHR2111*CHR2111), 
                                mean(xi2111*xi2211)/(n*CHR2111*CHR2211), mean(xi2111*xi2112)/(n*CHR2111*CHR2112), 
                                mean(xi2111*xi2212)/(n*CHR2111*CHR2212), mean(xi2111*xi2221)/(n*CHR2111*CHR2221),
                                mean(xi2211*xi1211)/(n*CHR2211*CHR1211), mean(xi2211*xi2111)/(n*CHR2211*CHR2111), 
                                mean(xi2211*xi2211)/(n*CHR2211*CHR2211), mean(xi2211*xi2112)/(n*CHR2211*CHR2112), 
                                mean(xi2211*xi2212)/(n*CHR2211*CHR2212), mean(xi2211*xi2221)/(n*CHR2211*CHR2221),
                                mean(xi2112*xi1211)/(n*CHR2112*CHR1211), mean(xi2112*xi2111)/(n*CHR2112*CHR2111), 
                                mean(xi2112*xi2211)/(n*CHR2112*CHR2211), mean(xi2112*xi2112)/(n*CHR2112*CHR2112), 
                                mean(xi2112*xi2212)/(n*CHR2112*CHR2212), mean(xi2112*xi2221)/(n*CHR2112*CHR2221),
                                mean(xi2212*xi1211)/(n*CHR2212*CHR1211), mean(xi2212*xi2111)/(n*CHR2212*CHR2111), 
                                mean(xi2212*xi2211)/(n*CHR2212*CHR2211), mean(xi2212*xi2112)/(n*CHR2212*CHR2112), 
                                mean(xi2212*xi2212)/(n*CHR2212*CHR2212), mean(xi2212*xi2221)/(n*CHR2212*CHR2221),
                                mean(xi2221*xi1211)/(n*CHR2221*CHR1211), mean(xi2221*xi2111)/(n*CHR2221*CHR2111), 
                                mean(xi2221*xi2211)/(n*CHR2221*CHR2211), mean(xi2221*xi2112)/(n*CHR2221*CHR2112), 
                                mean(xi2221*xi2212)/(n*CHR2221*CHR2212), mean(xi2221*xi2221)/(n*CHR2221*CHR2221)), nrow=6, ncol=6))
      names(est[[j]]) <- c("COEF", "TIME", "CHR", "VAR[CHR]", "LOG(CHR)", "VAR[LOG(CHR)]")
    
    }
    
  }
  
  #Return results
  return(est)
  
}



