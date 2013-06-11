### R code from vignette source 'Lunceford'

###################################################
### code chunk number 1: chunklibraries
###################################################



  #Libraries required
  require(survival)
  require(rms)
 

###################################################
### code chunk number 2: chunksimulatedata
###################################################
  #----------------------Functions Begins-----------------------------------------
  # Function input definitions
  # n= number of subjects
  # Remission/Consent indicator R is Bernoulli(pi.r)
  # B treatment indicator R=Z is Bernoulli(pi.z)
  # When R=0 (non responders), a survival time T*lambda was drawn from exponential(lambda)
  # When R=1, T*alpha was drawn from exponential(alpha)
  # T.star.11 post remission survival time under B1
#----------------------Changes made by XT 20130123--------------------------------
  # Remove definitions of some variables, e.g. use "length(which(R==1)) instead of defining new variable for the length
  # Use apply function instead of loop for T.star.12
  # Delete the use of L2=2.5 years because it has already been taken care of by the max of C
  # Change the notation R1 to R
  # Change the function name to be "sim.LDT.data"
  # Change the variable name Delta to delta


  sim.LDT.data<-function(n,L,max.c,pi.r,pi.z,lambda,alpha,beta1,beta2)
  { 
    #Functions needed to check errors:
    #Function to test if a number is an integer or not.
    is.wholenumber <-function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    
    #Function to check if a number is a valid probability value or not
    is.probability<-function(x) if(x>=0 && x<=1) TRUE else FALSE
    
    #Step 0: Error Checks:
    #To check the number of data points
    #Check for character values
    if (is.character(n)) stop("n has to be numeric")
    if (n<0 || !is.wholenumber(n)) stop("n must be a non-negative integer")
    # To check for L and max.c
    #Check for character values
    if (is.character(L)) stop("L has to be numeric")
    if (is.character(max.c)) stop("c has to be numeric")
    if (L<=0) stop("L must be a non-negative value")
    if (max.c<=0) stop("max.c must be a non-negative value")    
    # To check the probabilities
    if (is.character(pi.r)) stop("pi.r has to be numeric")
    if (is.character(pi.z)) stop("pi.z has to be numeric")
    if (!is.probability(pi.r)) stop("pi.r must be a valid probability value")    
    if (!is.probability(pi.z)) stop("pi.z must be a valid probability value")
    # To check for a valid lambda
    if (is.character(lambda)) stop("lambda has to be numeric")
    if (lambda<=0) stop("Lambda must be a positive parameter for the exponential distribution")
    # To check for a valid alpha
    if (is.character(alpha)) stop("alpha has to be numeric")
    if (alpha<=0) stop("alpha must be a positive parameter for the exponential distribution")    
    # To check for a valid beta1
    if (is.character(beta1)) stop("beta1 has to be numeric")
    if (beta1<=0) stop("beta1 must be a positive parameter for the exponential distribution")    
    # To check for a valid beta2
    if (is.character(beta2)) stop("beta2 has to be numeric")    
    
    
    
    #Step 1
    #Generate Remission/Consent indicator
      R<-rbinom(n,1,pi.r)
      
      #Generate B treatment indicator - Only for R==1. Because Z is only supposed to be generated 
      #for R=1, those non-responders will not be randomized again the second time, 
      #so for non-responders, Z is automatically equal to 0
      Z<-rep(0,n)
      Z[which(R==1)]<-rbinom(length(which(R==1)),1,pi.z)
      
      #Generate censoring
      C<-runif(n,min=0,max=max.c)
    
    #Step 2
      # When R=0, a survival time T*lambda was drawn from exponential(lambda)
      T.star.lambda<-rep(0,n)
      T.star.lambda[which(R==0)]<-rexp(length(which(R==0)),rate=lambda)
      
      # When R=1, T*alpha was drawn from exponential(alpha)
      T.star.alpha<-rep(0,n)
      T.star.alpha[which(R==1)]<-rexp(length(which(R==1)),rate=alpha)
    
      #T.star.11
      T.star.11<-rep(0,n)
      T.star.11[which(R==1)]<-rexp(length(which(R==1)),rate=exp(beta1))
      
      #T.star.12 
      T.star.12<-rep(0,n)
      T.star.12[which(R==1)]<-apply(as.array(T.star.11[which(R==1)]), 1, function(x) rexp(1,rate=exp(beta1+beta2*x)))
      
    #Step 3 - CAN OPTIMIZE HERE FOR SIMULAIONS
      T.11<-pmin((1-R)*T.star.lambda+R*(T.star.alpha+T.star.11),rep(L,n))
      T.12<-pmin((1-R)*T.star.lambda+R*(T.star.alpha+T.star.12),rep(L,n))
    
    #Step 4
      T<-(1-R)*T.11+R*(1-Z)*T.11+R*Z*T.12
      V<-pmin(T,C)
      delta<-as.numeric(T<C)
      
    #Step 5  
    data<-data.frame(R,Z,V,delta)
    return(data)
  }
  
