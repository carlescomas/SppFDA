C
C     Author: Carles Comas
C     Date: 15th desember 2021 
C
C     Function  Functional.mark.cor.f to obtain the Functional Mark Correlation Function (FMCF)
C     
C

       SUBROUTINE FunMarkCorf1(funpp,d,n,r,nr,nrf,minrf,maxrf,
     &                 kr,bw,lamd,cor,edge,Ar,FMCF)


       IMPLICIT real*8 (a-h,o-z)

C      PARAMETERS
C      i, j, i1 internal counting variables, integer
C      TestF, matrix with Test function combining functions of the point pattern, test function: L2 distances
C      d, pairwise distance. double precision d(n,n)
C      n, number of points of the point pattern, interger
C      r, sequence of values, for the discret set of distances where the FMCF is evaluated, double precision
C      nr, the length of r, integer
C      nrf, the legnth of rf, integer
C      kr, value of the kernel function, integer value
C      bw, bandwidth value, real*8
C      lamd, point intensity, double precision
C      cor, value of the isotropic correction for point i wrt point j
C      edge, value for the edge-effect correcion, interger 
C      Ar, area of the window of observation, double precision
C      FMCF, Functional Mark Correlation Function, output

       INTEGER  i, j, i1, n, nr, kr, edge, nrf  
       DOUBLE PRECISION d(n,n), r(nr), lamd(n), cor(n,n)
       DOUBLE PRECISION TestF(n,n), FMCF(nr), PCF(nr)
       DOUBLE PRECISION fxs(nrf), funpp(nrf,n),maxrf(n,n)
       REAL*8 Ar, bw, dij, kerns, ExpF, a, b, minrf
       DIMENSION kr(3), edge(2)
     
       
     

C      PAIR CORRELATION FUNCTION

       DO i1=1, nr
        PCF(i1)=0.0
         DO i=1,n
          DO j=1,n
           IF (j.ne.i) THEN
            dij=d(i,j)
              IF (kr(1).eq.1) THEN
c             kerns=1
               kerns=boxkernel((r(i1)-dij)/bw,bw)
              ELSE IF (kr(2).eq.1) THEN
               kerns=ekernel((r(i1)-dij)/bw,bw)
              ELSE IF (kr(3).eq.1) THEN
              kerns=qkernel((r(i1)-dij)/bw,bw)
              END IF

              IF (kerns.ne.0) THEN
C     none   
                IF (edge(1).eq.1) THEN
                 wij=kerns/(lamd(i)*lamd(j)*r(i1))
                 PCF(i1)=PCF(i1)+wij 
                END IF                          
C    isotropic
                IF (edge(2).eq.1) THEN                  
                 wij=(kerns*cor(i,j))/(lamd(i)*lamd(j)*r(i1))
                 PCF(i1)=PCF(i1)+wij
                END IF

              END IF
           END IF
         END DO
        END DO
        PCF(i1)=PCF(i1)/(2*3.141592654*Ar)
       END DO

C   TEST FUNCTION

        DO i=1, n
         DO j=1, n
           IF (j.ne.i) THEN 
             DO i1=1, nrf
             fxs(i1)=(funpp(i1,i)-funpp(i1,j))**2
             ENDDO
             TestF(i,j)=sqrt(SIMP(fxs, minrf, maxrf(i,j),nrf))
           ENDIF
         ENDDO
        ENDDO


C   FUNCTIONAL MARK CORRELATION FUNCTION

      ExpF=ETestF(n,TestF)

      DO i1=1, nr
        FMCF(i1)=0.0
         DO i=1,n
          DO j=1,n
           IF (j.ne.i) THEN
            dij=d(i,j)
              IF (kr(1).eq.1) THEN
               kerns=boxkernel((r(i1)-dij)/bw,bw)
              ELSE IF (kr(2).eq.1) THEN
               kerns=ekernel((r(i1)-dij)/bw,bw)
              ELSE IF (kr(3).eq.1) THEN
              kerns=qkernel((r(i1)-dij)/bw,bw)
              END IF

              IF (kerns.ne.0) THEN
C     none   
                IF (edge(1).eq.1) THEN
             wij1=(kerns*TestF(i,j))
             wij2=lamd(i)*lamd(j)*r(i1)*ExpF*PCF(i1)
             wij=wij1/wij2     
                 FMCF(i1)=FMCF(i1)+wij 
                END IF                          
C    isotropic
                IF (edge(2).eq.1) THEN                  
            wij1=kerns*TestF(i,j)*cor(i,j)
            wij2=lamd(i)*lamd(j)*r(i1)*ExpF*PCF(i1)
            wij=wij1/wij2
                 FMCF(i1)=FMCF(i1)+wij
                END IF

              END IF
           END IF
         END DO
        END DO
        FMCF(i1)=FMCF(i1)/(2*3.141592654*Ar)
       END DO

        RETURN
        END
       


       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
C
C     functions called by :
C     -----------------------------------------
C
C     * boxkernel, ekernel, qkernel
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C--------------------------------------------------------------------
C
C     boxkernel
C
C--------------------------------------------------------------------

       FUNCTION boxkernel(x,h)

       IMPLICIT REAL*8 (a-h,o-z)

       REAL*8 x, h

       IF (abs(x).le.1d0) THEN
         boxkernel=1d0/2d0
       ELSE
         boxkernel=0d0
       END IF
       boxkernel=boxkernel/h
       RETURN
       END

C--------------------------------------------------------------------
C
C     Epanechnikov kernel
C
C--------------------------------------------------------------------

       FUNCTION ekernel(x,h)

     
       IMPLICIT REAL*8 (a-h,o-z)
       REAL*8 x, h

       IF (abs(x).le.1) THEN
           ekernel=(3d0/(4d0*h))*(1-x**2)
       ELSE
          ekernel=0.0
       ENDIF

C       IF (abs(x).le.1d0) THEN
C           ekernel=(3d0/4d0)*(1-x**2)
C       ELSE
C           ekernel=0d0
C       END IF
C       ekernel=ekernel/h

        RETURN
        END

 

C--------------------------------------------------------------------
C
C     quartic (biweight) kernel
C
C--------------------------------------------------------------------

       FUNCTION qkernel(x,h)

       IMPLICIT REAL*8 (a-h,o-z)
       REAL*8 x, h

       IF (abs(x).le.1d0) THEN
           qkernel=(15d0/16d0)*(1-x**2)**2
       ELSE
           qkernel=0d0
       END IF
       qkernel=qkernel/h

       RETURN
       END

C--------------------------------------------------------------------
C
C     Expected Test Function
C
C--------------------------------------------------------------------

      FUNCTION ETestF(n,TestF)

      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER n
      DOUBLE PRECISION TestF(n,n)

       ETestF=0.0
       DO i=1,n
          DO j=1,n
           IF (j.ne.i) THEN
           ETestF=ETestF+TestF(i,j)
           ENDIF
         ENDDO
       ENDDO
       ETestF=ETestF/(n**2-n)
       RETURN
       END

C--------------------------------------------------------------------
C
C     Integral Function using Simpson
C
C--------------------------------------------------------------------

        FUNCTION SIMP(fxs,a,b,nrf)
        IMPLICIT REAL*8 (a-h,o-z)
        INTEGER nrf,c1,c2
        DOUBLE PRECISION fxs(nrf)
        REAL*8 a, b, h
        
         h=(b-a)/(nrf-1.0)	
         SIMP=3.0*(fxs(1)+fxs(nrf))/8.0+7.0*(fxs(2)+
     &    fxs(nrf-1))/6.0+23.0*(fxs(3)+fxs(nrf-2))/24.0
         c1=4
         c2=nrf-3
         DO i=c1, c2
          SIMP=SIMP+fxs(i)
         ENDDO
         SIMP=SIMP*h
         RETURN
         END




