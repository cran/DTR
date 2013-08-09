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
### code chunk number 2: chunktesting
###################################################
# fdata: the complete data frame for a two-stage randomization design
# t: a time point of interest

CHR.Wald.test <- function(fdata, # A complete data frame representing the data for two-stage randomization designs
                                 # fdata = data frame {X, R, Z, U, delta, V}
                                 # V represents covariates
                                 # There could be one covariate, or more than one covariates
                                 # The function does not allow the absence of covariates
                          t=quantile(fdata$U, 0.75) # A time of interest for comparison, optional, default is the 75th percentile of U
) {
    
  #Check for errors
  if (is.null(fdata$X)) stop("X can not be empty")
  if (is.null(fdata$R)) stop("R can not be empty")
  if (is.null(fdata$Z)) stop("Z can not be empty")  
  if (is.null(fdata$U)) stop("U can not be empty")  
  if (is.null(fdata$delta)) stop("delta can not be empty") 
  
  if (is.null(t)) stop("t not be empty")  
  if (is.character(t)) stop("t has to be numeric")
  if (t<0) stop("t must be a non-negative value") 

  #Retrive covariates
  V <- as.matrix(fdata[, ! names(fdata) %in% c("X", "R", "Z", "U", "delta")])

  #################################################
  #################################################
  #If no covariates: ERROR
  #################################################
  #################################################
  if(NCOL(V)==0) { stop("Covariates can not be empty") } 

  #################################################
  #################################################
  #If >= 1 covariate
  #################################################
  #################################################
  if(NCOL(V)>=1) {
    
    #Run the CHR.estimator function
    results <- CHR.estimator(fdata, t)[[1]]
    
  }  
  
  #Retrieve theta estimates and variance estiamtes
  EST <- matrix(as.numeric(results$"LOG(CHR)"), nrow=length(results$"LOG(CHR)"), ncol=1)
  VAR <- results$"VAR[LOG(CHR)]"
  
  #Test for H0: A1B1=A1B2=A2B1=A2B2
  #log(theta1211)=log(theta2111)=log(theta2211)=0
  D_overall <- matrix(c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0), nrow=3, ncol=6, byrow=T)
  test_overall <- t(D_overall %*% EST) %*% solve(D_overall %*% VAR %*% t(D_overall)) %*% (D_overall %*% EST)
  p_overall <- pchisq(test_overall, df=3, lower.tail=F)
  
  #Test for H0: A1B1=A1B2
  #log(theta1211)=0
  D_1112 <- matrix(c(1,0,0,0,0,0), nrow=1, ncol=6, byrow=T)
  test_1112 <- t(D_1112 %*% EST) %*% solve(D_1112 %*% VAR %*% t(D_1112)) %*% (D_1112 %*% EST)
  p_1112 <- pchisq(test_1112, df=1, lower.tail=F)
  
  #Test for H0: A1B1=A2B1
  #log(theta2111)=0
  D_1121 <- matrix(c(0,1,0,0,0,0), nrow=1, ncol=6, byrow=T)
  test_1121 <- t(D_1121 %*% EST) %*% solve(D_1121 %*% VAR %*% t(D_1121)) %*% (D_1121 %*% EST)
  p_1121 <- pchisq(test_1121, df=1, lower.tail=F)
  
  #Test for H0: A1B1=A2B2
  #log(theta2211)=0
  D_1122 <- matrix(c(0,0,1,0,0,0), nrow=1, ncol=6, byrow=T)
  test_1122 <- t(D_1122 %*% EST) %*% solve(D_1122 %*% VAR %*% t(D_1122)) %*% (D_1122 %*% EST)
  p_1122 <- pchisq(test_1122, df=1, lower.tail=F)
  
  #Test for H0: A1B2=A2B1
  #log(theta2112)=0
  D_1221 <- matrix(c(0,0,0,1,0,0), nrow=1, ncol=6, byrow=T)
  test_1221 <- t(D_1221 %*% EST) %*% solve(D_1221 %*% VAR %*% t(D_1221)) %*% (D_1221 %*% EST)
  p_1221 <- pchisq(test_1221, df=1, lower.tail=F)
  
  #Test for H0: A1B2=A2B2
  #log(theta2212)=0
  D_1222 <- matrix(c(0,0,0,0,1,0), nrow=1, ncol=6, byrow=T)
  test_1222 <- t(D_1222 %*% EST) %*% solve(D_1222 %*% VAR %*% t(D_1222)) %*% (D_1222 %*% EST)
  p_1222 <- pchisq(test_1222, df=1, lower.tail=F) 
  
  #Test for H0: A2B1=A2B2
  #log(theta2221)=0
  D_2122 <- matrix(c(0,0,0,0,0,1), nrow=1, ncol=6, byrow=T)
  test_2122 <- t(D_2122 %*% EST) %*% solve(D_2122 %*% VAR %*% t(D_2122)) %*% (D_2122 %*% EST)
  p_2122 <- pchisq(test_2122, df=1, lower.tail=F)
  
  #Return results
  TEST <- data.frame(c("A1B1=A1B2=A2B1=A2B2", "A1B1=A1B2", "A1B1=A2B1", "A1B1=A2B2", "A1B2=A2B1", "A1B2=A2B2", "A2B1=A2B2"),
                     round(c(test_overall, test_1112, test_1121, test_1122, test_1221, test_1222, test_2122),4),
                     format.pval(round(c(p_overall, p_1112, p_1121, p_1122, p_1221, p_1222, p_2122),4), eps=0.0001))
  names(TEST) <- c(paste("H0 (t=", round(t,2), ")", sep=""), "test statistic", "p")
  
  return(TEST)    
  
}  

