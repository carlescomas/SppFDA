#' Function to obtain the Functional Mark Correlation Function
#'
#' This function computes the Functional Mark Correlation Function as described in Comas et al (2009) and Comas et al. (2013)
#' We consider as a test function for the functional point pattern, L2 distances between functions
#'
#' @examples 
#'   To be implemented
#'  
#' @param xy, object of class ppp (point pattern)
#' @param r, sequence of values to evaluate the Functional mark correlation function
#' @param rf, sequence of values for the function associated to each point
#' @param funpp is a matrix that contains de functions for each point, funct[,1] is the function for point 1
#' @param kernel, kernel function for both the pair correlation function and the Functional mark correlation function
#' @param bw bandwidth for the kernel function for the pair correlation function and Functional mark correlation function
#' @param lambda, point intensity
#' @param stoyan value for the epanechnikov kernel, default stoyan=0.15
#' @param correction for edge-effect for both the pair correlation function and the Functional mark correlation function
#' @param simul, number of permutations for the Monte Carlo test, default simul=999
#' @param alpha, signification level of the Monte Carlo Test
#' @param verbose, show the simulation progression, default verbose=TRUE
#' @param Figure, display resulting FMCF and related max and min evelops. default Figure=TRUE
#' @param miny, maxy, plot parameters, ylim=c(miny,maxy)
#' @return 
#' @author Carles Comas \email{carles.comas@udl.cat}
#' @import spatstat
#' @export

EnvFunMarkCor<-function(xy, r, rf, funpp, stoyan=0.15, kernel = "epanechnikov", bw, lambda, correction = "isotropic",simul=299, alpha=0.05,
                        verbose=TRUE, Figure=TRUE, miny, maxy){

        
         out<-FunMarkCor(xy,rf=rf,funpp=funpp)
         n <- xy$n
         il <- seq(1, simul)

         Envel <- sapply(il, function(il) resample(il,n,verbose), simplify = "array")
         MaxEnv<-c(length(out$r))
         MinEnv<-c(length(out$r))
         k<-round(alpha*(n+1)/2)

          for(i in 1:length(out$r))   MaxEnv[i]<-sort(Envel[i,])[simul-k+1]
          for(i in 1:length(out$r))   MinEnv[i]<-sort(Envel[i,])[k]


       #  for(i in 1:length(out$r)) MaxEnv[i]<-max(Envel[i,])
       #  for(i in 1:length(out$r)) MinEnv[i]<-min(Envel[i,])
        
         
           if (missing(miny)) miny<-min(MinEnv)-min(MinEnv)*0.1
           if (missing(maxy)) maxy<-max(MaxEnv)+max(MaxEnv)*0.1

         if(Figure) {
         plot(out$r, out$FMCF, type="l", ylab=expression(paste(FMCF)),xlab=expression(paste(r)),lwd=2.0, col="black", cex.axis=1.20, cex.lab=1.40,
         ylim=c(miny,maxy), main="Functional Mark Correlation Function")
         lines(out$r, MaxEnv, col="red",lwd=2.0)
         lines(out$r, MinEnv, col="red",lwd=2.0)
         abline(h=1.0, col="green",lwd=1.5)
         }

        invisible(return(list(FMCF = out$FMCF, r = out$r, kernel = out$kernel, lambda = out$lambda, MaxEnv = MaxEnv, MinEnv = MinEnv)))
}

resample<-function(i1,n,verbose){
 
         res<-sample(seq(1,n,1),n)
         funppsam=array(NA,c(length(rf),n))
         dim(funppsam)=c(length(rf),n)

         for(i in 1:n) funppsam[,i]<-funpp[,res[i]]
         out<-FunMarkCor(xy,rf=rf,funpp=funppsam)
         out=out$FMCF
        
       if (verbose) {
         cat("simulation", paste(i1))
         cat("\n")
        flush.console()
       }
      
     return(out = out)
}









