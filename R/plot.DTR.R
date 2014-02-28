###################################################
### code chunk number 1: chunklibraries
###################################################

#Libraries required
require(ggplot2)
require(grid)
require(proto)

###################################################
### code chunk number 2: chunkfunction
###################################################

###################
# The functions: stairstepn, stat_stepribbon, and StatStepribbon were downloaded from:
# https://groups.google.com/forum/?fromgroups=#!topic/ggplot2/9cFWHaH1CPs
###################

# Modified version of ggplot:::stairstep function
stairstepn <- function( data, direction="hv", yvars="y" ) {
  direction <- match.arg( direction, c( "hv", "vh" ) )
  data <- as.data.frame( data )[ order( data$x ), ]
  n <- nrow( data )
  
  if ( direction == "vh" ) {
    xs <- rep( 1:n, each = 2 )[ -2 * n ]
    ys <- c( 1, rep( 2:n, each = 2 ) )
  } else {
    ys <- rep( 1:n, each = 2 )[ -2 * n ]
    xs <- c( 1, rep( 2:n, each = 2))
  }
  
  data.frame(
    x = data$x[ xs ]
    , data[ ys, yvars, drop=FALSE ]
    , data[ xs, setdiff( names( data ), c( "x", yvars ) ), drop=FALSE ]
  ) 
}

stat_stepribbon <- function( mapping=NULL, data=NULL, geom="ribbon", position="identity" ) {
  StatStepribbon$new( mapping=mapping, data=data, geom=geom, position=position )
}

StatStepribbon <- proto(ggplot2:::Stat, {
  objname <- "stepribbon"
  desc <- "Stepwise area plot"
  desc_outputs <- list(
    x = "stepped independent variable",
    ymin = "stepped minimum dependent variable",
    ymax = "stepped maximum dependent variable"
  )
  required_aes <- c( "x", "ymin", "ymax" )
  
  default_geom <- function(.) GeomRibbon
  default_aes <- function(.) aes( x=..x.., ymin = ..y.., ymax=Inf )
  
  calculate <- function( ., data, scales, direction = "hv", yvars = c( "ymin", "ymax" ), ...) {
    stairstepn( data = data, direction = direction, yvars = yvars )
  }
  
  examples <- function(.) {
    DF <- data.frame( x = 1:3, ymin = runif( 3 ), ymax=rep( Inf, 3 ) )
    ggplot( DF, aes( x=x, ymin=ymin, ymax=ymax ) ) + stat_stepribbon()
    
  }
  
})


plot.DTR <- function(x, # A complete data frame representing the data for two-stage randomization designs
                    confidence.interval=FALSE, # Plot confidence intreval or not, default no confidence interval
                    xlab="Time", # x axis label
                    ylab="Survival probability", # y axis label
                    line.color=c("black", "grey40", "grey60", "grey80"), # Line colors for A1B1, A1B2, A2B1, A2B2 in order
                    legend.position="right", ... # Position of the legend
) {
  
  #Check for errors
  
  if (confidence.interval != FALSE & confidence.interval != TRUE & 
        confidence.interval != T & confidence.interval != F) stop("confidence.interval input can not be recognized")
  
  if (is.null(xlab)) stop("xlab can not be empty")
  if (is.character(xlab)==FALSE) stop("xlab has to be character")
  
  if (is.null(ylab)) stop("ylab can not be empty")
  if (is.character(ylab)==FALSE) stop("ylab has to be character")
  
  if (is.null(line.color)) stop("line.color can not be empty")
  
  if (is.null(legend.position)) stop("legend.position can not be empty")
  
  #Reformat results for ggplot
  group <- c(rep("A1B1", length(x$time)), rep("A1B2", length(x$time)), 
          rep("A2B1", length(x$time)), rep("A2B2", length(x$time))) 
  time <- rep(x$time, 4)
  surv <- c(x$SURV11, x$SURV12, x$SURV21, x$SURV22)
  se <- c(x$SE11, x$SE12, x$SE21, x$SE22)
  plot.result <- data.frame(group, time, surv, se)
    
  #Obtain x scale
  xpool <- c(0.05, 0.1, seq(0.2,2,0.2), seq(0.5, 10, 0.5), seq(10, 50, 5), seq(50, 100, 10), seq(100, 1000, 50))
  xdiff <- (xpool-max(x$time)/10)[(xpool-max(x$time)/10)>=0]
  x.scale <- max(x$time)/10 + unique(xdiff[which(xdiff==min(xdiff))])
  
  g <- ggplot(data=plot.result) +
    geom_step(aes(x=time, y=surv, color=group), size=1.2) +
    scale_color_manual(values=line.color) + 
    scale_x_continuous(breaks=seq(0, x.scale*10, x.scale)) + 
    xlab(xlab) + ylab(ylab) +
    theme_bw() +
    theme(title=element_text(size=15), axis.text=element_text(size=15), 
          legend.position=legend.position, legend.title=element_blank(),
          legend.text=element_text(size=15), legend.key.size=unit(1, "cm"))
  
  
  #Plot the estimates from 0 to L without confidence interval
  if(confidence.interval==FALSE) {
    g <- g
  } else {
    g <- g +
      geom_ribbon(aes(x = time, ymin = surv-1.96*se, ymax =
                        surv+1.96*se, fill=group), alpha = 0.10, stat='stepribbon',  direction="hv") +
      scale_fill_manual(values=line.color)  +
      geom_step(aes(x=time, y=surv-1.96*se, color=group), linetype='dotted', size=0.5) +
      geom_step(aes(x=time, y=surv+1.96*se, color=group), linetype='dotted', size=0.5)
  }
  
  #Plot
  g
  
}
