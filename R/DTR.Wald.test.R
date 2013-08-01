###################################################
### Reference:
### Lunceford JK, Davidian M, Tsiatis AA: Estimation of survival distributions of treatment 
### policies in two-stage randomization designs in clinical trials. Biometrics 58:48-57, 2002
###################################################

###################################################
### Reference:
### Guo X, Tsiatis AA: A weighted risk set estimator for survival distributions in two-stage 
### randomization designs with censored survival data. Int. J. Biostatistics 1:1-15, 2005
###################################################

###################################################
### code chunk number 1: chunklibraries
###################################################

#Libraries required
require(survival)

###################################################
### code chunk number 2: chunktesting
###################################################

DTR.Wald.test <- function(fdata, # A complete data frame representing the data for two-stage randomization designs
                                 # fdata = data frame {X, R, Z, U, delta} for method = "LDT"
                                 # fdata = data frame {X, TR, R, Z, U, delta} for method = "WRSE"
                          t=quantile(fdata$U, 0.75), # A time of interest for comparison, optional, default is the 75th percentile of U
                          method="LDT", # The estimates used for comparisons
                                        # method = "LDT" for Lunceford et al. (2002)
                                        # method = "WRSE" for Guo and Tsiatis (2005)
                          rmean=FALSE, # rmean=FALSE for comparisons based on estimates at t
                                       # rmean=TRUE for comparisons based on restricted mean survivals
                          L=.Machine$double.xmax # # Optional restricted survival time L
) {
  
  #Check for errors
  if (is.null(fdata$X)) stop("X can not be empty")
  if (is.null(fdata$R)) stop("R can not be empty")
  if (is.null(fdata$Z)) stop("Z can not be empty")  
  if (is.null(fdata$U)) stop("U can not be empty")  
  if (is.null(fdata$delta)) stop("delta can not be empty") 
  
  if (is.null(t)) stop("t can not be empty") 
  if (is.character(t)) stop("t has to be numeric")
  if (t<0) stop("t must be a non-negative value") 
  
  if (method != "LDT" & method != "WRSE") stop("method input can not be recognized")
  
  if (rmean != FALSE & rmean != TRUE & rmean != T & rmean != F) stop("rmean input can not be recognized")

  if (is.null(L)) stop("L can not be empty")
  if (is.character(L)) stop("L has to be numeric")
  if (L<=0) stop("L must be a positive value")
  
  #Define all the survival estimates and variance/covariance estimates
  SURV11 <- SURV12 <- SURV21 <- SURV22 <- NA
  SE11 <- SE12 <- SE21 <- SE22 <- NA
  COV1112 <- COV2122 <- NA

  #Comparisons based on LDT estimates at given t
  if(method=="LDT" & rmean==FALSE) {
    
    #Run LDT.estimator to obtain the defined estimates A1 arm
    est1 <- LDT.estimator(data=fdata[which(fdata$X==0),names(fdata) %in% c("R", "Z", "U", "delta")],t=t,L)
    SURV11 <- est1$SURV1
    SURV12 <- est1$SURV2
    SE11 <- est1$SE1
    SE12 <- est1$SE2
    COV1112 <- est1$COV12
    
    #Run LDT.estimator to obtain the defined estimates A2 arm
    est2 <- LDT.estimator(data=fdata[which(fdata$X==1),names(fdata) %in% c("R", "Z", "U", "delta")],t=t,L)
    SURV21 <- est2$SURV1
    SURV22 <- est2$SURV2
    SE21 <- est2$SE1
    SE22 <- est2$SE2
    COV2122 <- est2$COV12
    
  }

  #Comparisons based LDT restricted mean survivals
  if(method=="LDT" & rmean==TRUE) {
    
    #Run LDT.estimator to obtain the defined estimates A1 arm
    est1 <- LDT.mean.estimator(data=fdata[which(fdata$X==0),names(fdata) %in% c("R", "Z", "U", "delta")],L)
    SURV11 <- est1$MEAN1
    SURV12 <- est1$MEAN2
    SE11 <- est1$MSE1
    SE12 <- est1$MSE2
    COV1112 <- est1$MCOV12
    
    #Run LDT.estimator to obtain the defined estimates A2 arm
    est2 <- LDT.mean.estimator(data=fdata[which(fdata$X==1),names(fdata) %in% c("R", "Z", "U", "delta")],L)
    SURV21 <- est2$MEAN1
    SURV22 <- est2$MEAN2
    SE21 <- est2$MSE1
    SE22 <- est2$MSE2
    COV2122 <- est2$MCOV12
    
  }

  #Comparisons based on WRSE estimates at given t
  if(method=="WRSE") {
    
    #Check for errors
    if (is.null(fdata$TR)) stop("TR can not be empty")
    
    #Run LDT.estimator to obtain the defined estimates A1 arm
    est1 <- WRSE.estimator(data=fdata[which(fdata$X==0),names(fdata) %in% c("TR", "R", "Z", "U", "delta")],t=t)
    SURV11 <- est1$SURV1
    SURV12 <- est1$SURV2
    SE11 <- est1$SE1
    SE12 <- est1$SE2
    COV1112 <- est1$COV12
    
    #Run LDT.estimator to obtain the defined estimates A2 arm
    est2 <- WRSE.estimator(data=fdata[which(fdata$X==0),names(fdata) %in% c("TR", "R", "Z", "U", "delta")],t=t)
    SURV21 <- est2$SURV1
    SURV22 <- est2$SURV2
    SE21 <- est2$SE1
    SE22 <- est2$SE2
    COV2122 <- est2$COV12
    
  }
    
  #Combine survival estimates as column vector
  SURV <- matrix(c(SURV11, SURV12, SURV21, SURV22), nrow=4, ncol=1)
  #Combine variance estimates as matrix
  VAR <- matrix(c(SE11^2, COV1112, 0, 0, COV1112, SE12^2, 0, 0,
                  0, 0, SE21^2, COV2122, 0, 0, COV2122, SE22^2), nrow=4, ncol=4)
  
  #Test for H0: S11=S12=S21=S22  
  D_overall <- matrix(c(1,1,1,-1,0,0,0,-1,0,0,0,-1), nrow=3, ncol=4)
  test_overall <- t(D_overall %*% SURV) %*% solve(D_overall %*% VAR %*% t(D_overall)) %*% (D_overall %*% SURV)
  p_overall <- pchisq(test_overall, df=3, lower.tail=F)  
  
  #Test for H0: S11=S12
  D_1112 <- matrix(c(1,-1,0,0), nrow=1, ncol=4)
  test_1112 <- t(D_1112 %*% SURV) %*% solve(D_1112 %*% VAR %*% t(D_1112)) %*% (D_1112 %*% SURV)
  p_1112 <- pchisq(test_1112, df=1, lower.tail=F) 
  
  #Test for H0: S11=S21
  D_1121 <- matrix(c(1,0,-1,0), nrow=1, ncol=4)
  test_1121 <- t(D_1121 %*% SURV) %*% solve(D_1121 %*% VAR %*% t(D_1121)) %*% (D_1121 %*% SURV)
  p_1121 <- pchisq(test_1121, df=1, lower.tail=F) 
  
  #Test for H0: S11=S22
  D_1122 <- matrix(c(1,0,0,-1), nrow=1, ncol=4)
  test_1122 <- t(D_1122 %*% SURV) %*% solve(D_1122 %*% VAR %*% t(D_1122)) %*% (D_1122 %*% SURV)
  p_1122 <- pchisq(test_1122, df=1, lower.tail=F)  
  
  #Test for H0: S12=S21
  D_1221 <- matrix(c(0,1,-1,0), nrow=1, ncol=4)
  test_1221 <- t(D_1221 %*% SURV) %*% solve(D_1221 %*% VAR %*% t(D_1221)) %*% (D_1221 %*% SURV)
  p_1221 <- pchisq(test_1221, df=1, lower.tail=F) 
  
  #Test for H0: S12=S22
  D_1222 <- matrix(c(0,1,0,-1), nrow=1, ncol=4)
  test_1222 <- t(D_1222 %*% SURV) %*% solve(D_1222 %*% VAR %*% t(D_1222)) %*% (D_1222 %*% SURV)
  p_1222 <- pchisq(test_1222, df=1, lower.tail=F) 
  
  #Test for H0: S21=S22
  D_2122 <- matrix(c(0,0,1,-1), nrow=1, ncol=4)
  test_2122 <- t(D_2122 %*% SURV) %*% solve(D_2122 %*% VAR %*% t(D_2122)) %*% (D_2122 %*% SURV)
  p_2122 <- pchisq(test_2122, df=1, lower.tail=F)  

  #Return restuls
  TEST <- data.frame(c("A1B1=A1B2=A2B1=A2B2", "A1B1=A1B2", "A1B1=A2B1", "A1B1=A2B2", "A1B2=A2B1", "A1B2=A2B2", "A2B1=A2B2"),
                     round(c(test_overall, test_1112, test_1121, test_1122, test_1221, test_1222, test_2122),4),
                     format.pval(round(c(p_overall, p_1112, p_1121, p_1122, p_1221, p_1222, p_2122),4),eps=0.0001))
  if(rmean==FALSE) { names(TEST) <- c(paste("H0 (t=", round(t,2), ")", sep=""), "test statistic", "p") }
  if(method=="WRSE" & rmean==TRUE) { names(TEST) <- c(paste("H0 (t=", round(t,2), ")", sep=""), "test statistic", "p") }
  if(method=="LDT" & rmean==TRUE) { names(TEST) <- c("H0 (mean)", "test statistic", "p") }
  return(TEST)    
    
}  
  
