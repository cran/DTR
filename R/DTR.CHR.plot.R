###################################################
### Plot CHR estimates
### Reference:
### Tang X, Wahed AS: Cumulative hazard ratio estimation for treatment regimes in
### sequentially randomized clinical trials. Statistics in Biosciences, 2013 [Epub ahead of print]
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

DTR.CHR.plot <- function(fdata, # A complete data frame representing the data for two-stage randomization designs
                                # fdata = data frame {X, R, Z, U, delta, V}
                        log.CHR=FALSE, # If log.CHR=FALSE, CHR is plotted;
                                       # If log.CHR=TRUE, natural logarithm of CHR is plotted
                        confidence.interval=FALSE, # Plot confidence intreval or not, default no confidence interval
                        xlab="Time", # x axis label
                        line.color=c("black", "grey30", "grey50", "grey60", "grey70", "grey80"), # Line colors
                        legend.position="right" # Position of the legend
                                          
) {
  
  #Check for errors
  if (is.null(fdata$X)) stop("X can not be empty")
  if (is.null(fdata$R)) stop("R can not be empty")
  if (is.null(fdata$Z)) stop("Z can not be empty")  
  if (is.null(fdata$U)) stop("U can not be empty")  
  if (is.null(fdata$delta)) stop("delta can not be empty")  
  
  if (log.CHR != FALSE & log.CHR != TRUE & log.CHR != T & log.CHR != F) stop("log.CHR input can not be recognized")

  if (confidence.interval != FALSE & confidence.interval != TRUE & confidence.interval != T & confidence.interval != F) stop("confidence.interval input can not be recognized")
  
  if (is.null(xlab)) stop("xlab can not be empty")
  if (is.character(xlab)==FALSE) stop("xlab has to be character")
  
  if (is.null(line.color)) stop("line.color can not be empty")
  
  if (is.null(legend.position)) stop("legend.position can not be empty")  

  #Obtain the time sequence for plotting
  t <- seq(0, max(fdata$U), max(fdata$U)/100)
  
  #Obtain x scale
  xpool <- c(0.05, 0.1, seq(0.2,2,0.2), seq(0.5, 10, 0.5), seq(10, 50, 5), seq(50, 100, 10), seq(100, 1000, 50))
  xdiff <- (xpool-max(fdata$U)/10)[(xpool-max(fdata$U)/10)>=0]
  x.scale <- max(fdata$U)/10 + unique(xdiff[which(xdiff==min(xdiff))])
  
  #Retrive covariates
  V <- as.matrix(fdata[, ! names(fdata) %in% c("X", "R", "Z", "U", "delta")])
  
  #Define all the cumulative hazard ratio estimates and variance/covariance estimates
  CHR1211 <- CHR2111 <- CHR2211 <- CHR2112 <- CHR2212 <- CHR2221 <- rep(NA, length(t))
  LCHR1211 <- LCHR2111 <- LCHR2211 <- LCHR2112 <- LCHR2212 <- LCHR2221 <- rep(NA, length(t))
  SE1211 <- SE2111 <- SE2211 <- SE2112 <- SE2212 <- SE2221 <- rep(NA, length(t))
  LSE1211 <- LSE2111 <- LSE2211 <- LSE2112 <- LSE2212 <- LSE2221 <- rep(NA, length(t))
  
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
    results <- CHR.estimator(fdata, t)
    
  }
  
  #Obtain the estimates
  for(i in 1:length(t)) {

    #Results for each time point
    iresults <- results[[i]]

    #Cumulative hazard ratio estimates
    CHR1211[i] <- iresults$CHR$CHR1211; CHR2111[i] <- iresults$CHR$CHR2111; CHR2211[i] <- iresults$CHR$CHR2211;
    CHR2112[i] <- iresults$CHR$CHR2112; CHR2212[i] <- iresults$CHR$CHR2212; CHR2221[i] <- iresults$CHR$CHR2221;

    #Log cumulative hazard ratio estimates
    LCHR1211[i] <- iresults$"LOG(CHR)"$LogCHR1211; LCHR2111[i] <- iresults$"LOG(CHR)"$LogCHR2111; LCHR2211[i] <- iresults$"LOG(CHR)"$LogCHR2211;
    LCHR2112[i] <- iresults$"LOG(CHR)"$LogCHR2112; LCHR2212[i] <- iresults$"LOG(CHR)"$LogCHR2212; LCHR2221[i] <- iresults$"LOG(CHR)"$LogCHR2221;
    
    #Standard error estiamtes for CHR
    SE1211[i] <- sqrt(iresults$"VAR[CHR]"[1,1]); SE2111[i] <- sqrt(iresults$"VAR[CHR]"[2,2]); SE2211[i] <- sqrt(iresults$"VAR[CHR]"[3,3])
    SE2112[i] <- sqrt(iresults$"VAR[CHR]"[4,4]); SE2212[i] <- sqrt(iresults$"VAR[CHR]"[5,5]); SE2221[i] <- sqrt(iresults$"VAR[CHR]"[6,6])

    #Standard error estiamtes for log CHR
    LSE1211[i] <- sqrt(iresults$"VAR[LOG(CHR)]"[1,1]); LSE2111[i] <- sqrt(iresults$"VAR[LOG(CHR)]"[2,2]); LSE2211[i] <- sqrt(iresults$"VAR[LOG(CHR)]"[3,3])
    LSE2112[i] <- sqrt(iresults$"VAR[LOG(CHR)]"[4,4]); LSE2212[i] <- sqrt(iresults$"VAR[LOG(CHR)]"[5,5]); LSE2221[i] <- sqrt(iresults$"VAR[LOG(CHR)]"[6,6])
    
  }
  
  #Reformat results for ggplot
  group = c(rep("A1B2 vs. A1B1", length(t)), rep("A2B1 vs. A1B1", length(t)), rep("A2B2 vs. A1B1", length(t)), 
            rep("A2B1 vs. A1B2", length(t)), rep("A2B2 vs. A1B2", length(t)), rep("A2B2 vs. A2B1", length(t)))
  time = rep(t, 6)
  chr = c(CHR1211, CHR2111, CHR2211, CHR2112, CHR2212, CHR2221)
  lchr = c(LCHR1211, LCHR2111, LCHR2211, LCHR2112, LCHR2212, LCHR2221)
  se = c(SE1211, SE2111, SE2211, SE2112, SE2212, SE2221)
  lse = c(LSE1211, LSE2111, LSE2211, LSE2112, LSE2212, LSE2221)
  plot.result <- data.frame(group, time, chr, lchr, se, lse)
  
  #Get rid of NAs
  plot.result <- plot.result[!is.na(plot.result$chr),]
  
  #Plot the estimates without confidence interval
  if(log.CHR==FALSE & confidence.interval==FALSE) {
    
    g <- ggplot(plot.result, aes(x=time, y=chr, color=group)) +
          geom_line(size=1.5) +
          geom_hline(yintercept=1, linetype="dashed", size=1) +
          scale_color_manual(values=line.color) +
          scale_x_continuous(breaks=seq(0, x.scale*10, x.scale)) + 
          xlab(xlab) + ylab("Cumulative hazard ratio") +
          theme_bw() +
          theme(title=element_text(size=18), axis.text=element_text(size=15), 
            legend.position=legend.position, legend.title=element_blank(),
            legend.text=element_text(size=14), legend.key.size=unit(1.3, "cm"))
   
  } 
  
  #Plot the estimates with confidence interval
  if(log.CHR==FALSE & confidence.interval==TRUE) {
    
    g <- ggplot(plot.result, aes(x=time, y=chr, color=group)) +
      geom_line(size=1.5) + 
      geom_ribbon(aes(ymin=chr-1.96*se,ymax=chr+1.96*se),linetype='blank',alpha=0.1)+
      geom_hline(yintercept=1, linetype="dashed", size=1) +
      scale_color_manual(values=line.color) +
      scale_x_continuous(breaks=seq(0, x.scale*10, x.scale)) + 
      xlab(xlab) + ylab("Cumulative hazard ratio") +
      theme_bw() +
      theme(title=element_text(size=18), axis.text=element_text(size=15), 
            legend.position=legend.position, legend.title=element_blank(),
            legend.text=element_text(size=14), legend.key.size=unit(1.3, "cm"))
            
  }
  
  #Plot the log estimates without confidence interval
  if(log.CHR==TRUE & confidence.interval==FALSE) {
    
    g <- ggplot(plot.result, aes(x=time, y=lchr, color=group)) +
      geom_line(size=1.5) +
      geom_hline(yintercept=0, linetype="dashed", size=1) +
      scale_color_manual(values=line.color) +
      scale_x_continuous(breaks=seq(0, x.scale*10, x.scale)) + 
      xlab(xlab) + ylab("Log cumulative hazard ratio") +
      theme_bw() +
      theme(title=element_text(size=18), axis.text=element_text(size=15), 
            legend.position=legend.position, legend.title=element_blank(),
            legend.text=element_text(size=14), legend.key.size=unit(1.3, "cm"))
    
  } 
  
  #Plot the log estimates with confidence interval
  if(log.CHR==TRUE & confidence.interval==TRUE) {
    
    g <- ggplot(plot.result, aes(x=time, y=lchr, color=group)) +
      geom_line(size=1.5) + 
      geom_ribbon(aes(ymin=lchr-1.96*lse,ymax=lchr+1.96*lse),linetype='blank',alpha=0.1)+
      geom_hline(yintercept=0, linetype="dashed", size=1) +
      scale_color_manual(values=line.color) +
      scale_x_continuous(breaks=seq(0, x.scale*10, x.scale)) + 
      xlab(xlab) + ylab("Log cumulative hazard ratio") +
      theme_bw() +
      theme(title=element_text(size=18), axis.text=element_text(size=15), 
            legend.position=legend.position, legend.title=element_blank(),
            legend.text=element_text(size=14), legend.key.size=unit(1.3, "cm"))
    
  }
  
  #Plot
  g
  
}
  



