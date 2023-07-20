#' Function to obtain the Auto- and cross-function mark summary characteristics for the mark variogram and mark correlation envelopes
#'
#' This function computes an estimator of this function as described in Eckardt, Comas and Mateu (2023)
#' This function admits a parallel version, where the total number of cores -1 are used
#'
#' @examples 
#'   To be implemented
#'  
#' @param xy, object of class ppp (point pattern)
#' @param r, sequence of values to evaluate the  Auto- and cross-function mark summary characteristics
#' @param rf, sequence of values for the function associated to each point
#' @param f1 is a matrix that contains de functions for each point, funct[,1] is the function for point 1 (first set of functions)
#' @param f2 is a matrix that contains de functions for each point, funct[,1] is the function for point 1 (second set of functions)
#' @param kernel, kernel function for both the pair correlation function and the Functional mark correlation function
#' @param bw bandwidth for the kernel function for the pair correlation function and Functional mark correlation function
#' @param lambda, point intensity
#' @param TestFunc, Test function. This can be TestFunc="variogram" or TestFunc="correlation"
#' @param stoyan value for the epanechnikov kernel, default stoyan=0.15
#' @param correction for edge-effect for both the pair correlation function and the Functional mark correlation function
#' @param nsim, number of permutations for the Monte Carlo test, default nsim=299
#' @param nrank,Integer. Rank of the envelope value amongst the ‘nsim’ simulated values. A rank of 1 means that the minimum and maximum simulated values will be used
#' @param alpha, signification level of the Monte Carlo Test
#' @param Parallel, If 'TRUE' the algorithm will be considered assuming a parallel approach using the parallel package and assuming the total number of cores of your computer - 1.
#' @param verbose, show the simulation progression, only for noParallel codification, default verbose=TRUE
#' @param Figure, display resulting FMCF and related max and min evelops. default Figure=TRUE
#' @param miny, maxy, plot parameters, ylim=c(miny,maxy)
#' @return 
#' @author Carles Comas \email{carles.comas@udl.cat}
#' @import spatstat
#' @import parallel
#' @export

EnvFunCrossMark<-function(xy, r, rf, f1, f2, stoyan=0.15, 
                        TestFunc = "variogram", kernel = "epanechnikov", bw, lambda, correction = "isotropic", nsim=199, nrank = 5, 
                        verbose=TRUE, Figure=TRUE, miny, maxy, Parallel=TRUE){


         win <- xy$window 
         n <- xy$n
         kernel <- kernel
         areaW <- area.owin(win)
         if(missing(lambda)) lambda <- n/areaW
         lambda<-lambda
             
         if (missing(r)){
         rect <- as.rectangle(win)
         maxd <- min(diff(rect$xrange), diff(rect$yrange))/4
         r <- seq(maxd*0.1, maxd, len = 301)[-1]
         r <- sort(r)
         }

         if(r[1] == 0) r <- r[-1]
         if (!inherits(r, "numeric")) stop("r should be a numeric value")
         if (length(r)<=1) stop("r should be a sequence of values")


         if (missing(f2)) f2<-f1
        
         noParallel=TRUE
         if (Parallel) noParallel=FALSE

         out<-FunCrossMark(xy,r,rf=rf,f1=f1,f2=f2,kernel=kernel,TestFunc=TestFunc,lambda=lambda)
                          
         if (noParallel){
         n <- xy$n
         il <- seq(1, nsim)
         Envel <- sapply(il, function(il) resample(il,n,f2,r,verbose,kernel=kernel,TestFunc=TestFunc,lambda=lambda), simplify = "array")
         }

         if(Parallel){
         npp <- xy$n
         n<-nsim 
         no_cores <- detectCores() - 1
         cl <- makeCluster(no_cores,type="FORK")  
         Envel <- parSapply(cl,1:n, function(x) ParResample(n=x,npp,f2,r,kernel=kernel,TestFunc=TestFunc,lambda=lambda))
         stopCluster(cl)
         }


         MaxEnv<-c(length(out$r))
         MinEnv<-c(length(out$r))
         

          alpha<-2*nrank/(nsim+1)
          paste("alpha value:",paste(alpha,collapse=","))


          for(i in 1:length(out$r))   MaxEnv[i]<-sort(Envel[i,])[nsim-nrank+1]
          for(i in 1:length(out$r))   MinEnv[i]<-sort(Envel[i,])[nrank]
         
           minval<-min(min(out$FMCF),min(MinEnv))
           maxval<-max(max(out$FMCF),max(MinEnv))
           if (missing(miny)) miny<-minval-minval*0.05
           if (missing(maxy)) maxy<-maxval+maxval*0.05


         if(Figure) {
         plot(out$r, out$FMCF, type="l",ylab="", xlab=expression(paste(r)),lwd=2.0, col="black", cex.axis=1.20, cex.lab=1.40,
         ylim=c(miny,maxy))
         lines(out$r, MaxEnv, col="red",lwd=2.0)
         lines(out$r, MinEnv, col="red",lwd=2.0)
         abline(h=1.0, col="green",lwd=1.5)
         }

        invisible(return(list(FMCF = out$FMCF, r = out$r, kernel = out$kernel, lambda = out$lambda, MaxEnv = MaxEnv, MinEnv = MinEnv)))
}

resample<-function(i1,n,f2,r,verbose,kernel=kernel,TestFunc=TestFunc,lambda=lambda){
        
         res<-sample(seq(1,n,1),n) 
         f1sam=array(NA,c(length(rf),n))
         dim(f1sam)=c(length(rf),n)
         f2sam<-f1sam

         for(i in 1:n) {
            f1sam[,i]<-f1[,res[i]]
            f2sam[,i]<-f2[,res[i]]
         }
         out<-FunCrossMark(xy,r,rf=rf,f1=f1sam,f2=f2sam,kernel=kernel,TestFunc=TestFunc,lambda=lambda)
        out=out$FMCF
         
       if (verbose) {
         cat("simulation", paste(i1))
         cat("\n")
        flush.console()
       }
      
     return(out = out)
}

ParResample<-function(n,npp,f2,r,kernel=kernel,TestFunc=TestFunc,lambda=lambda){
  i<-n 
  set.seed(NULL)
         res<-sample(seq(1,npp,1),npp)
         f1sam=array(NA,c(length(rf),npp))
         dim(f1sam)=c(length(rf),npp)
         f2sam<-f1sam

         for(i in 1:npp) {
           f1sam[,i]<-f1[,res[i]]
           f2sam[,i]<-f2[,res[i]]
         }

          out<-FunCrossMark(xy,r,rf=rf,f1=f1sam,f2=f2sam,kernel=kernel,TestFunc=TestFunc,lambda=lambda)
          out=out$FMCF
         
      
     return(out = out)
}









