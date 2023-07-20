#' Function to obtain the Auto- and cross-function mark summary characteristics for the mark variogram and mark correlation
#'
#' This function computes an estimator of this function as described in Eckardt, Comas and Mateu (2023)
#'
#' @examples 
#'  using Data_Duke
#'  Functiona.mark.cor(
#'  
#' @param xy, object of class ppp (point pattern)
#' @param r, sequence of values to evaluate the Auto- and cross-function mark summary characteristics
#' @param rf, sequence of values for the function associated to each point
#' @param f1 is a matrix that contains the first set of functions for each point, funct1[,1] is the function of set 1 for point 1
#' @param f2 is a matrix that contains the second set of functions for each point, funct1[,1] is the function of set 2 for point 1
#' @param kernel, kernel function for both the pair correlation function and the function mark summary characteristic
#' @param TestFunc, Test function. This can be TestFunc="variogram" or TestFunc="correlation"
#' @param bw bandwidth for the kernel function for the pair correlation function and Functional mark correlation function
#' @param lambda, point intensity
#' @param stoyan value for the epanechnikov kernel, default stoyan=0.15
#' @param correction for edge-effect for both the pair correlation function and the Functional mark correlation function
#' @return 
#' @author Carles Comas \email{carles.comas@udl.cat}
#' @useDynLib SppFDA, .registration=TRUE
#' @import spatstat
#' @export

FunCrossMark <- function(xy, r, rf, f1, f2, stoyan=0.15, TestFunc = "variogram", kernel = "epanechnikov", bw, lambda, correction = "isotropic"){
       
         verifyclass(xy, "ppp")
         correc <- c("none", "isotropic")
        
         if (missing(f1)) stop("the f1 matrix should be given")
         if (!inherits(f1, "matrix")) stop("funct1 should be a matrix")
         if (missing(f2)) f2<-f1
     
        

         if (missing(rf)) stop("the rf sequence should be given")
         if (!inherits(rf, "numeric")) stop("rf should be a numeric valuess")
         if (length(rf)<7) stop("rf should be a sequence of values larger than 6")

         indc <- match(correction, correc, nomatch = NA)
         if (any(cm <- is.na(indc))){
          nocm <- paste("unrecognised correction method:",paste(dQuote(correction[cm]),collapse=","))
          stop(nocm,call. = FALSE)
         }
         
         indc <- unique(indc)	
         edge <- rep(0, 2)
         edge[indc] <- 1	
         
         corr <- c("variogram", "correlation")
         indr <- match(TestFunc, corr, nomatch = NA)
         if (any(kmr <- is.na(indr))){
         nokmr <- paste("unrecognised test function:",paste(dQuote(TestFunc[kmr]),collapse=","))
         stop(nokmr,call. = FALSE)
         }
         indr <- unique(indr)
         corr2 <- rep(0, 2)
         corr2[indr] <- 1

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
            if((kernel != "epanechnikov") | length(lambda)>1) {
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

         kern <- c(kernel = kernel, bw = bw)
         
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

         na11<-c()
          for(i in 1:n){ 
            a11<-f1[,i]
            na11[i]=(nrf+1)-length(a11[!is.infinite(a11)])
          }

          na12<-c()
          for(i in 1:n){ 
            a12<-f2[,i]
            na12[i]=(nrf+1)-length(a12[!is.infinite(a12)])
          }


          minrf<-array(NA,c(n,n))
          dim(minrf)=c(n,n)

          for(i in 1:n){
           for(j in 1:n){
            na<-max(na11[i],na12[j])
            minrf[i,j]<-na
           }
          }

       
            na21<-c()
          for(i in 1:n){ 
            a21<-f1[,i]
            na21[i]=length(a21[!is.na(a21)])
          }

          na22<-c()
          for(i in 1:n){ 
            a22<-f2[,i]
            na22[i]=length(a22[!is.na(a22)])
          }

          
          maxrf<-array(NA,c(n,n))
          dim(maxrf)=c(n,n)

          for(i in 1:n){
           for(j in 1:n){
            na3<-min(na21[i],na22[j])
             maxrf[i,j]<-na3
           }
          }  


        f1[is.na(f1)] <- 1.0
        f1[is.infinite(f1)] <- 1.0
        f2[is.na(f2)] <- 1.0
        f2[is.infinite(f2)] <- 1.0
         h<-rf[2]-rf[1]
       
         FMCF<-rep(10,nr)
         storage.mode(FMCF) <- "double"
                   
out<-.Fortran("funcrossmarksub", f1=as.double(f1), f2 = as.double(f2),
                       d = as.double(d), n = as.integer(n), r = as.double(r),
                       nr = as.integer(nr), nrf=as.integer(nrf), minrf=as.double(minrf),
                       maxrf=as.double(maxrf), ker2 = as.integer(ker2),corr2 = as.integer(corr2),
                       bw = as.double(bw), h = as.double(h),lambda = as.double(lambda),
                       cor = as.double(cor), edge = as.integer(edge),
                       Ar=as.double(Ar), out=(FMCF), PACKAGE = "SppFDA")
                      
invisible(return(list(FMCF = out$out, r = r, kernel = kern, lambda = lambda)))
}











