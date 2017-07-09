************************************************************************
*
*     Subroutine SUBSTIT               Called by: CONSER, REACT
*
*     simultaneous solution of a system of linear algebraic equations,
*     where the coefficient matrix is tridiagonal.
*
*     preforms the forward and backward substitution phases of the
*     Thomas algorithm.  This substitution step is preceded by matrix
*     decomposition (subroutine DECOMP).
*
*     The algorithm has been modified to include work vectors AWORK and
*     BWORK.  These vectors are included to preserve the contents of A
*     and B.
*
*     (for more information, see 'Numerical Methods for Engineers',
*      Chapra and Canale, 1988, pages 286-288)
*
*************************************************************************
      SUBROUTINE SUBSTIT(IMAX,CONC,AWORK,BWORK,C,D,J)

      INCLUDE 'fmodules.inc'
*
*     subroutine arguments
*
      INTEGER*4 IMAX,J
      DOUBLE PRECISION C(*),D(*)
      DOUBLE PRECISION AWORK(MAXSEG,*),BWORK(MAXSEG,*),CONC(MAXSEG,*)
*
*     Local variables
*
      INTEGER*4 I
*
*     forward substitution
*
      DO 10 I = 2,IMAX
         D(I) = D(I) - AWORK(I,J)*D(I-1)
 10   CONTINUE
*
*     backward substitution
*
      CONC(IMAX,J) = D(IMAX) / BWORK(IMAX,J)
      DO 20 I = IMAX-1,1,-1
         CONC(I,J) = (D(I) - C(I)*CONC(I+1,J))/BWORK(I,J)
 20   CONTINUE
   
      RETURN
      END
