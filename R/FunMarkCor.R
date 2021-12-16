#' Function to obtain the Functional Mark Correlation Function
#'
#' This function computes the Functional Mark Correlation Function as described in Comas et al (2009) and Comas et al. (2013)
#' We consider as a test function for the functional point pattern, L2 distances between functions
#'
#' @examples 
#'  using Data_Duke
#'  Functiona.mark.cor(
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
#' @return 
#' @author Carles Comas \email{carles.comas@udl.cat}
#' @useDynLib SppFDA, .registration=TRUE
#' @import spatstat
#' @export

FunMarkCor <- function(xy, r, rf, funpp, stoyan=0.15, kernel = "epanechnikov", bw, lambda, correction = "isotropic"){
        
         verifyclass(xy, "ppp")
         correc <- c("none", "isotropic")

         if (missing(funpp)) stop("the funpp matrix should be given")
         if (!inherits(funpp, "matrix")) stop("funct should be a matrix")

         if (missing(rf)) stop("the rf sequence should be given")
         if (!inherits(rf, "numeric")) stop("rf should be a numeric value")
         if (length(rf)<7) stop("rf should be a sequence of values larger than 6")

         indc <- match(correction, correc, nomatch = NA)
         if (any(cm <- is.na(indc))){
          nocm <- paste("unrecognised correction method:",paste(dQuote(correction[cm]),collapse=","))
          stop(nocm,call. = FALSE)
         }
         
         indc <- unique(indc)	
         edge <- rep(0, 2)
         edge[indc] <- 1	

         kern <- c("rectangular", "epanechnikov", "biweight")
         indk <- match(kernel, kern, nomatch = NA)
         if (any(km <- is.na(indk))){
         nokm <- paste("unrecognised kernel function:",paste(dQuote(kernel[km]),collapse=","))
         stop(nokm,call. = FALSE)
         }

         indk <- unique(indk)
         ker2 <- rep(0, 3)
         ker2[indk] <- 1

         win <- xy$window 
         n <- xy$n
         areaW <- area.owin(win)

         if(missing(lambda)) lambda <- n/areaW

         if (missing(bw) && (kernel == "epanechnikov") && length(lambda)== 1) {
          bw <- stoyan/sqrt(lambda)
         }
 
         if (missing(bw)){
            if((kernel != "epanechnikov") | lenght(lambda)>1) {
           bw <- bw.pcf(X = xy, kernel = kernel)[1] 
            }
         }
         if (!inherits(bw, "numeric")) stop("bw should be a numeric value")

        

         if (missing(r)){
         rect <- as.rectangle(win)
         maxd <- min(diff(rect$xrange), diff(rect$yrange))/4
         r <- seq(maxd*0.1, maxd, len = 301)[-1]
         r <- sort(r)
         }

         if(r[1] == 0) r <- r[-1]
         if (!inherits(r, "numeric")) stop("r should be a numeric value")
         if (length(r)<=1) stop("r should be a sequence of values")

         kernel <- c(kernel = kernel, bw = bw)
         
         n <- xy$n
         il <- seq(1, n)
         nr <- length(r)
         areaW <- area.owin(win)
         Ar<-areaW

       
         if (length(lambda) == 1) lambda <- rep(lambda, n)
    
         d <- pairdist(xy)
         if(correction == "isotropic") cor <- edge.Ripley(xy, d)
         if(correction == "none") cor <- 1

       
          nrf<-length(rf)
          minrf<-rf[1]

          na1<-c()
          for(i in 1:n){ 
            a1<-funpp[,i]
            na1[i]=length(a1[!is.na(a1)])
          }

          
          maxrf<-array(NA,c(n,n))
          dim(maxrf)=c(n,n)

          for(i in 1:n){
           for(j in 1:n){
            na<-min(na1[i],na1[j])
            maxrf[i,j]<-rf[na]
           }
          }


         FMCF<-rep(10,nr)
         storage.mode(FMCF) <- "double"
         
out<-.Fortran("FunMarkCorf1", funpp=as.double(funpp),d = as.double(d),
                       n = as.integer(n), r = as.double(r),
                       nr = as.integer(nr), nrf=as.integer(nrf), minrf=as.double(minrf),
                       maxrf=as.double(maxrf), ker2 = as.integer(ker2),
                       bw = as.double(bw), lambda = as.double(lambda),
                       cor = as.double(cor), edge = as.integer(edge),
                       Ar=as.double(Ar), out=(FMCF), PACKAGE = "SppFDA")
                      

invisible(return(list(FMCF = out$out, r = r, kernel = kernel, lambda = lambda)))
}











