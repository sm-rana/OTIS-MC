************************************************************************
*
*     Subroutine PREPROC2            Called by: PREPROC, SETUP3
*
*     compute the time-invariant parameter groups used within the
*     finite difference subroutines (PREPROC3, CONSER/REACT).  As these
*     groups are constant in time, they may be computed here on a one-
*     time basis.
*
************************************************************************
      SUBROUTINE PREPROC2(IMAX,NSOLUTE,TIMEB,LAMBDA2,LAM2DT,LAMHAT2,
     &                    LHAT2DT,LAMHAT,LHATDT,SGROUP,SGROUP2,RHOLAM,
     &                    CSBACK)
      INCLUDE 'fmodules.inc'
*
*     subroutine arguments
*
      INTEGER*4 IMAX,NSOLUTE
      DOUBLE PRECISION TIMEB
      DOUBLE PRECISION LAMBDA2(MAXSEG,*),LAM2DT(MAXSEG,*),
     &                 LAMHAT2(MAXSEG,*),LHAT2DT(MAXSEG,*),
     &                 LAMHAT(MAXSEG,*),LHATDT(MAXSEG,*),
     &                 SGROUP(MAXSEG,*),SGROUP2(MAXSEG,*),
     &                 RHOLAM(MAXSEG,*),CSBACK(MAXSEG,*)
*
*     local variables
*
      INTEGER*4 I,J
*
*     preprocess
*
      DO 20 I = 1,IMAX
         DO 10 J = 1,NSOLUTE
            LAM2DT(I,J) = LAMBDA2(I,J) / TIMEB
            LHATDT(I,J) = LAMHAT(I,J) / TIMEB
            LHAT2DT(I,J) = LAMHAT2(I,J) / TIMEB
            SGROUP(I,J) = RHOLAM(I,J) / (2.D0 + LHATDT(I,J))
            SGROUP2(I,J) = LHAT2DT(I,J)*CSBACK(I,J)
 10      CONTINUE
 20   CONTINUE

      RETURN
      END
