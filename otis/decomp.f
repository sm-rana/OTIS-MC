************************************************************************
*
*     Subroutine DECOMP               Called by: PREPROC3
*
*     prepare to compute the simultaneous solution to a system of linear
*     algebraic equations, where the coefficient matrix is tridiagonal.
*
*     perform the decomposition phase of the Thomas algorithm.  The
*     decomposition phase is followed by forward and backward substition
*     (subroutine SUBSTIT).  Note that the solution to the system of
*     equations is provided by SUBSTIT.
*
*     The algorithm has been modified to include work vectors AWORK
*     and BWORK.  These vectors are included so that the contents of A
*     and B are preserved.
*
*     (for more information, see 'Numerical Methods for Engineers',
*      Chapra and Canale, 1988, pages 286-288)
*
*************************************************************************
      SUBROUTINE DECOMP(IMAX,A,B,C,AWORK,BWORK,J)
      INCLUDE 'fmodules.inc'
*
*     subroutine arguments
*
      INTEGER*4 IMAX,J
      DOUBLE PRECISION A(*),C(*)
      DOUBLE PRECISION B(MAXSEG,*),BWORK(MAXSEG,*),AWORK(MAXSEG,*)
*
*     Local variables
*
      INTEGER*4 I
*
*     L.U. decompositon of the coefficient matrix
*
        ! Masud: J = 1,NSOLUTE --> J = 1
        ! Masud: Both AWORK and BWORK is filled up in this subroutine using A and B
      BWORK(1,J) = B(1,J)

      DO 10 I = 2,IMAX
         AWORK(I,J) = A(I)/BWORK(I-1,J)
         BWORK(I,J) = B(I,J) - AWORK(I,J)*C(I-1)
 10   CONTINUE

      RETURN
      END
