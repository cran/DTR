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
### code chunk number 2: chunkbeta
###################################################

updateBeta <- function(beta, # Current coefficient(s) for covariate(s)
                       V, # Covariate(s)
                       U, # Observed survival time
                       delta, # Censoring indicator
                       w11, # Weights for A1B1
                       w12, # Weights for A1B2
                       w21, # Weights for A2B1
                       w22 # Weights for A2B2
) {

  #Total number of subjects
  n <- length(U)
  
  #If only 1 covariate
  if(NCOL(V)==1) { 
    
    #Change V from matrix to numeric
    V <- as.numeric(V)
    
    #Calcualte exp(beta*V)
    e <- exp(beta*V)
    
    #Define ui and derivative of U: dui
    ui <- 0
    dui <- 0
    
    for(i in 1:length(U)) {
      
      ind <- as.numeric(U >= U[i])
      
      #Calculate s0j      
      s00 <- sum(ind * w11 * e) / n
      s01 <- sum(ind * w12 * e) / n
      s02 <- sum(ind * w21 * e) / n
      s03 <- sum(ind * w22 * e) / n
      
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
      
      #Calculate v1?_bar
      if(s00 != 0) v10_bar <- s10/s00 else v10_bar <- 0
      if(s01 != 0) v11_bar <- s11/s01 else v11_bar <- 0
      if(s02 != 0) v12_bar <- s12/s02 else v12_bar <- 0
      if(s03 != 0) v13_bar <- s13/s03 else v13_bar <- 0
      
      #calculate v2?_bar
      if(s00 != 0) v20_bar <- s20/s00 else v20_bar <- 0
      if(s01 != 0) v21_bar <- s21/s01 else v21_bar <- 0
      if(s02 != 0) v22_bar <- s22/s02 else v22_bar <- 0
      if(s03 != 0) v23_bar <- s23/s03 else v23_bar <- 0     
      
      #Calculate each ui
      ui <- ui + delta[i]*w11[i]*(V[i] - v10_bar) + 
        delta[i]*w12[i]*(V[i] - v11_bar) +
        delta[i]*w21[i]*(V[i] - v12_bar) +
        delta[i]*w22[i]*(V[i] - v13_bar)
      
      #Calculate each dui
      dui <- dui + delta[i]*w11[i]*(v10_bar^2 - v20_bar) +
        delta[i]*w12[i]*(v11_bar^2 - v21_bar) +
        delta[i]*w21[i]*(v12_bar^2 - v22_bar) +
        delta[i]*w22[i]*(v13_bar^2 - v23_bar)
      
    }
      
  }
  
  #If more than 1 covariates
  if(NCOL(V)>1) { 
    
    #Calcualte exp(beta*V)
    e <- as.numeric(exp(matrix(beta, nrow=1, ncol=NCOL(V)) %*% t(V)))
    
    #Define ui and derivative of U: dui
    ui <- matrix(0, nrow=NCOL(V), ncol=1)
    dui <- matrix(0, nrow=NCOL(V), ncol=NCOL(V))
    
    for(i in 1:length(U)) {
      
      ind <- as.numeric(U >= U[i])
      
      #Calculate s0j      
      s00 <- sum(ind * w11 * e) / n
      s01 <- sum(ind * w12 * e) / n
      s02 <- sum(ind * w21 * e) / n
      s03 <- sum(ind * w22 * e) / n
      
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
      if(s00 != 0) v10_bar <- s10/s00 else v10_bar <- matrix(0, nrow=NCOL(V), ncol=1)
      if(s01 != 0) v11_bar <- s11/s01 else v11_bar <- matrix(0, nrow=NCOL(V), ncol=1)
      if(s02 != 0) v12_bar <- s12/s02 else v12_bar <- matrix(0, nrow=NCOL(V), ncol=1)
      if(s03 != 0) v13_bar <- s13/s03 else v13_bar <- matrix(0, nrow=NCOL(V), ncol=1)
      
      #calculate z2?_bar
      if(s00 != 0) v20_bar <- s20/s00 else v20_bar <- matrix(0, nrow=NCOL(V), ncol=NCOL(V))
      if(s01 != 0) v21_bar <- s21/s01 else v21_bar <- matrix(0, nrow=NCOL(V), ncol=NCOL(V))
      if(s02 != 0) v22_bar <- s22/s02 else v22_bar <- matrix(0, nrow=NCOL(V), ncol=NCOL(V))
      if(s03 != 0) v23_bar <- s23/s03 else v23_bar <- matrix(0, nrow=NCOL(V), ncol=NCOL(V))      
      
      #Calculate each ui
      ui <- ui + delta[i]*w11[i]*(matrix(V[i,], nrow=NCOL(V), ncol=1) - v10_bar) + 
        delta[i]*w12[i]*(matrix(V[i,], nrow=NCOL(V), ncol=1) - v11_bar) +
        delta[i]*w21[i]*(matrix(V[i,], nrow=NCOL(V), ncol=1) - v12_bar) +
        delta[i]*w22[i]*(matrix(V[i,], nrow=NCOL(V), ncol=1) - v13_bar)
      
      #Calculate each dui
      dui <- dui + delta[i]*w11[i]*(v10_bar %*% t(v10_bar) - v20_bar) +
        delta[i]*w12[i]*(v11_bar %*% t(v11_bar) - v21_bar) +
        delta[i]*w21[i]*(v12_bar %*% t(v12_bar) - v22_bar) +
        delta[i]*w22[i]*(v13_bar %*% t(v13_bar) - v23_bar)
      
    }
    
  }
    
  #Update beta
    
  beta_up <- beta - solve(dui) %*% ui    
  return(beta_up)
    
}

