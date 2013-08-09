###################################################
### Plot LDT estimates
### Reference:
### Lunceford JK, Davidian M, Tsiatis AA: Estimation of survival distributions of treatment 
### policies in two-stage randomization designs in clinical trials. Biometrics 58:48-57, 2002
###################################################

###################################################
### Plot WRSE estimates
### Reference:
### Guo X, Tsiatis AA: A weighted risk set estimator for survival distributions in two-stage 
### randomization designs with censored survival data. Int. J. Biostatistics 1:1-15, 2005
###################################################

###################################################
### code chunk number 1: chunklibraries
###################################################

#Libraries required
require(survival)
require(ggplot2)
require(grid)

###################################################
### code chunk number 2: chunkplotting
###################################################

DTR.surv.plot <- function(fdata, # A complete data frame representing the data for two-stage randomization designs
                            # fdata = data frame {X, R, Z, U, delta} for method = "LDT"
                            # fdata = data frame {X, TR, R, Z, U, delta} for method = "WRSE"
                          method="LDT", # The estimates used for plotting
                                   # method = "LDT" for Lunceford et al. (2002)
                                   # method = "WRSE" for Guo and Tsiatis (2005)
                          L=.Machine$double.xmax, # Optional restricted survival time L
                          confidence.interval=FALSE, # Plot confidence intreval or not, default no confidence interval
                          xlab="Time", # x axis label
                          ylab="Survival probability", # y axis label
                          line.color=c("black", "grey40", "grey60", "grey80"), # Line colors for A1B1, A1B2, A2B1, A2B2 in order
                          legend.position="right" # Position of the legend
                                          
) {
  
  #Check for errors
  if (is.null(fdata$X)) stop("X can not be empty")
  if (is.null(fdata$R)) stop("R can not be empty")
  if (is.null(fdata$Z)) stop("Z can not be empty")  
  if (is.null(fdata$U)) stop("U can not be empty")  
  if (is.null(fdata$delta)) stop("delta can not be empty") 
  
  if (method != "LDT" & method != "WRSE") stop("method input can not be recognized")
  
  if (is.null(L)) stop("L can not be empty")
  if (is.character(L)) stop("L has to be numeric")
  if (L<=0) stop("L must be a positive value")
  
  if (confidence.interval != FALSE & confidence.interval != TRUE & confidence.interval != T & confidence.interval != F) stop("confidence.interval input can not be recognized")

  if (is.null(xlab)) stop("xlab can not be empty")
  if (is.character(xlab)==FALSE) stop("xlab has to be character")
  
  if (is.null(ylab)) stop("ylab can not be empty")
  if (is.character(ylab)==FALSE) stop("ylab has to be character")
  
  if (is.null(line.color)) stop("line.color can not be empty")

  if (is.null(legend.position)) stop("legend.position can not be empty")
  
  #Obtain the time sequence for plotting
  t <- seq(0, max(fdata$U), max(fdata$U)/100)
  
  #Obtain x scale
  xpool <- c(0.05, 0.1, seq(0.2,2,0.2), seq(0.5, 10, 0.5), seq(10, 50, 5), seq(50, 100, 10), seq(100, 1000, 50))
  xdiff <- (xpool-max(fdata$U)/10)[(xpool-max(fdata$U)/10)>=0]
  x.scale <- max(fdata$U)/10 + unique(xdiff[which(xdiff==min(xdiff))])
  
  #Define all the survival estimates and variance/covariance estimates
  SURV11 <- SURV12 <- SURV21 <- SURV22 <- rep(NA, length(t))
  SE11 <- SE12 <- SE21 <- SE22 <- rep(NA, length(t))
  
  if(method=="LDT") {

    print("Calculating LDT estimates...")
    #Run LDT.estimator to obtain the defined estimates A1 arm
    est1 <- LDT.estimator(data=fdata[which(fdata$X==0),names(fdata) %in% c("R", "Z", "U", "delta")],t=t,L)
    SURV11 <- est1$SURV1
    SURV12 <- est1$SURV2
    SE11 <- est1$SE1
    SE12 <- est1$SE2
    
    #Run LDT.estimator to obtain the defined estimates A2 arm
    est2 <- LDT.estimator(data=fdata[which(fdata$X==1),names(fdata) %in% c("R", "Z", "U", "delta")],t=t,L)
    SURV21 <- est2$SURV1
    SURV22 <- est2$SURV2
    SE21 <- est2$SE1
    SE22 <- est2$SE2
        
  }
  
  if(method=="WRSE") {

    #Check TR
    if (is.null(fdata$TR)) stop("TR can not be empty")
    
    print("Calculating WRSE estimates...")
    #Run LDT.estimator to obtain the defined estimates A1 arm
    est1 <- WRSE.estimator(data=fdata[which(fdata$X==0),names(fdata) %in% c("TR", "R", "Z", "U", "delta")],t=t)
    SURV11 <- est1$SURV1
    SURV12 <- est1$SURV2
    SE11 <- est1$SE1
    SE12 <- est1$SE2
    
    #Run LDT.estimator to obtain the defined estimates A2 arm
    est2 <- WRSE.estimator(data=fdata[which(fdata$X==1),names(fdata) %in% c("TR", "R", "Z", "U", "delta")],t=t)
    SURV21 <- est2$SURV1
    SURV22 <- est2$SURV2
    SE21 <- est2$SE1
    SE22 <- est2$SE2    
        
  }
  
  #Reformat results for ggplot
  group = c(rep("A1B1", length(t)), rep("A1B2", length(t)), rep("A2B1", length(t)), rep("A2B2", length(t)))
  time = rep(t, 4)
  surv = c(SURV11, SURV12, SURV21, SURV22)
  se = c(SE11, SE12, SE21, SE22)
  plot.result <- data.frame(group, time, surv, se)

  #Plot the estimates from 0 to L without confidence interval
  if(confidence.interval==FALSE) {
    
    g <- ggplot(plot.result, aes(x=time, y=surv, color=group)) +
          geom_line(size=1.5) +
          scale_color_manual(values=line.color) +
          scale_x_continuous(breaks=seq(0, x.scale*10, x.scale)) + 
          xlab(xlab) + ylab(ylab) +
          theme_bw() +
          theme(title=element_text(size=18), axis.text=element_text(size=15), 
            legend.position=legend.position, legend.title=element_blank(),
            legend.text=element_text(size=14), legend.key.size=unit(1.3, "cm"))
   
  } else {
    
    g <- ggplot(plot.result, aes(x=time, y=surv, color=group)) +
      geom_line(size=1.5) + 
      geom_ribbon(aes(ymin=surv-1.96*se,ymax=surv+1.96*se),linetype='blank',alpha=0.1)+
      scale_color_manual(values=line.color) +
      scale_x_continuous(breaks=seq(0, x.scale*10, x.scale)) + 
      xlab(xlab) + ylab(ylab) +
      theme_bw() +
      theme(title=element_text(size=18), axis.text=element_text(size=15), 
            legend.position=legend.position, legend.title=element_blank(),
            legend.text=element_text(size=14), legend.key.size=unit(1.3, "cm"))
            
  }
    
  #Plot
  g
  
}
  



