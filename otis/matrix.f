************************************************************************
*     
*     Subroutine MATRIX               Called by: SSDIFF
*
*     compute the simultaneous solution to a system of linear algebraic
*     equations, where the coefficient matrix is tridiagonal. 
*
*     The algorithm used to solve the tridiagonal system is known
*     as the Thomas algorithm.  The algorithm has been modified to
*     include work vector AWORK.  This vector is included so that the
*     contents of A are preserved.
*
*     (for more information, see 'Numerical Methods for Engineers', 
*      Chapra and Canale, 1988, pages 286-288)
*
*************************************************************************
      SUBROUTINE MATRIX(IMAX,CONC,A,B,C,D,J)

      INCLUDE 'fmodules.inc'
*
*     subroutine arguments
*
      INTEGER*4 IMAX,J
      DOUBLE PRECISION A(*),B(*),C(*),D(*),CONC(MAXSEG,*)
*
*     Local variables
*
      INTEGER*4 I
      DOUBLE PRECISION AWORK(MAXSEG)
*
*     Thomas algorithm
*
      DO 10 I = 2,IMAX
         AWORK(I) = A(I)/B(I-1)
         B(I) = B(I) - AWORK(I)*C(I-1)
         D(I) = D(I) - AWORK(I)*D(I-1)
 10   CONTINUE

      CONC(IMAX,J) = D(IMAX) / B(IMAX)
      DO 20 I = IMAX-1,1,-1
         CONC(I,J) = (D(I) - C(I)*CONC(I+1,J))/B(I)
 20   CONTINUE
   
      RETURN
      END
