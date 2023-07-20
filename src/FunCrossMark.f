C
C     Author: Carles Comas
C     Date: 3th FEBRUARY 2023 
C
C     Function  Functional.mark.cor.f to obtain the Functional Mark Correlation Function (FMCF)
C     
C

       SUBROUTINE funcrossmarksub(f1,f2,d,n,r,nr,nrf,minrf,maxrf,
     &                 kr,corr,bw,h,lamd,cor,edge,Ar,FMCF)


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

       INTEGER  i, j, i1, n, nr, kr, corr, edge, nrf, nn  
       DOUBLE PRECISION d(n,n), r(nr), lamd(n), cor(n,n)
       DOUBLE PRECISION TestF(n,n), FMCF(nr), PCF(nr)
       DOUBLE PRECISION fxs(nrf), f1(nrf,n)
       DOUBLE PRECISION maxrf(n,n), f2(nrf,n), minrf(n,n)
       REAL*8 Ar, bw, dij, kerns, ExpF, a, b, h
       DIMENSION kr(3), edge(2),corr(2)
     
       
       
     


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


       DO i=1, n
         DO j=1, n
          TestF(i,j)=0.0
         ENDDO
       ENDDO

C   TEST FUNCTION PAPER TEST 2009

C        DO i=1, n
C         DO j=1, n
C           IF (j.ne.i) THEN 
C             DO i1=1, nrf
C             fxs(i1)=(f1(i1,i)-f1(i1,j))**2
C             ENDDO
C             TestF(i,j)=sqrt(SIMP(fxs, minrf(i,j), maxrf(i,j),nrf))
C           ENDIF
C         ENDDO
C        ENDDO


c    TEST FUNCTIONS

C   FUNCTIONAL MARK VARIOGRAM FUNCTION

      IF (corr(1).eq.1) THEN
       DO i=1, n
           DO j=1, n
           IF (j.ne.i) THEN 
            DO i1=1, nrf
              fxs(i1)=0.5*(f1(i1,i)-f2(i1,j))**2
            ENDDO
            TestF(i,j)=SIMP(fxs, minrf(i,j), maxrf(i,j),nrf,h)
           ENDIF
           ENDDO
          ENDDO
      ENDIF

C   FUNCTIONAL MARK CORRELATION FUNCTION
      IF (corr(2).eq.1) THEN
        DO i=1, n
           DO j=1, n 
           IF (j.ne.i) THEN 
            DO i1=1, nrf
              fxs(i1)=f1(i1,i)*f2(i1,j)
            ENDDO
             TestF(i,j)=SIMP(fxs, minrf(i,j), maxrf(i,j),nrf,h) 
            ENDIF
           ENDDO
         ENDDO
      ENDIF

   

C   MEAN FUNCTIONAL MARK CORRELATION FUNCTION
C      IF(corr(2).eq.1) THEN
      ExpF=ETestF(n,TestF)
C      ENDIF

C      IF(corr(1).eq.1) THEN
C      ExpF=1.0
C      ENDIF

C   GENERAL FUNCTIONAL MARK CORRELATION FUNCTIONS

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
c           IF (j.ne.i) THEN
           ETestF=ETestF+TestF(i,j)
c           ENDIF
         ENDDO
       ENDDO
       ETestF=ETestF/(n**2)
       RETURN
       END


C--------------------------------------------------------------------
C
C     Integral Function using Simpson 1/3 rule
C
C--------------------------------------------------------------------
        FUNCTION SIMP(fxs,a,b,nrf,h)
        IMPLICIT REAL*8 (a-h,o-z)
        INTEGER nrf, i2, i3
        DOUBLE PRECISION fxs(nrf),coef(nrf)
        REAL*8 h, ni, sum1, a, b

       
        DO i=1,nrf
           coef(i)=1
        ENDDO    
  
         SIMP=0.0
        IF(a.lt.(b-3)) THEN
        ni=b-a+1
c        h=(b-a)/(ni-1)
       
        sum1=0


        IF (MOD(ni,2.0) == 0) THEN

          i2=a+1
          i3=b-2
           DO i=i2,i3
             IF (MOD(i,2)== 0) THEN
                coef(i)=4
             ELSE
                coef(i)=2
             ENDIF
            sum1=coef(i)*fxs(i)+sum1
           ENDDO 
           i2=b-1
           i3=b
           i4=a
            SIMP=h/3*(sum1+fxs(i2)+fxs(i4))+h/2*(fxs(i2)+fxs(i3))
         ELSE
      
          i2=a+1
          i3=b-1
          i4=a
          i5=b
          DO i=i2,i3
             IF (MOD(i,2)== 0) THEN
                coef(i)=4
             ELSE
                coef(i)=2
             ENDIF
            sum1=coef(i)*fxs(i)+sum1
          ENDDO 
           SIMP=h/3*(sum1+fxs(i5)+fxs(i4))
         ENDIF
        ENDIF

         SIMP=SIMP*1

c          IF(SIMP == SIMP) THEN
c             SIMP=SIMP*1
c          ELSE
c             SIMP=0.0
c          ENDIF

         RETURN
         END
       


c        FUNCTION SIMP(fxs,a,b,nrf)
c        IMPLICIT REAL*8 (a-h,o-z)
c        INTEGER nrf,c1,c2
c        DOUBLE PRECISION fxs(nrf)
c        REAL*8 a, b, h
        
c         h=(b-a)/(nrf-1.0)	
c         SIMP=3.0*(fxs(1)+fxs(nrf))/8.0+7.0*(fxs(2)+
c     &    fxs(nrf-1))/6.0+23.0*(fxs(3)+fxs(nrf-2))/24.0
c         c1=4
c         c2=nrf-3
c         DO i=c1, c2
c          SIMP=SIMP+fxs(i)
c         ENDDO
c         SIMP=SIMP*h
c         RETURN
c         END




