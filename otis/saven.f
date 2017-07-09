************************************************************************
*
*     Subroutine SAVEN            Called by: MAININIT, DYNAMIC, OTISRUN,
*                                            MAINRUN, DYNAMIC2
*
*     save flow variables from previous time step (time step 'N') and
*     preprocess parameter groups.
*
************************************************************************
      SUBROUTINE SAVEN(IMAX,NSOLUTE,AFACE,AFACEN,AREA,AREAN,GAMMA,
     &                 GAMMAN,Q,QLATIN,QLATINN,QN,CLATIN,CLATINN,AN,CN,
     &                 A,C,BTERMS,BTERMSN,TIMEB,TGROUP,BN,TWOPLUS,ALPHA,
     &                 IGROUP)
      INCLUDE 'fmodules.inc'
*
*     subroutine arguments
* SM: IMAX is the number of segments
      INTEGER*4 IMAX,NSOLUTE
      DOUBLE PRECISION TIMEB
      DOUBLE PRECISION AFACE(*),AFACEN(*),AREA(*),AREAN(*),GAMMA(*),
     &                 GAMMAN(*),Q(*),QLATIN(*),QLATINN(*),QN(*),AN(*),
     &                 CN(*),A(*),C(*),ALPHA(*)
      DOUBLE PRECISION CLATIN(MAXSEG,*),CLATINN(MAXSEG,*),
     &                 TGROUP(MAXSEG,*),BN(MAXSEG,*),TWOPLUS(MAXSEG,*),
     &                 BTERMS(MAXSEG,*),BTERMSN(MAXSEG,*),
     &                 IGROUP(MAXSEG,*)
*
*     local variables
*
      INTEGER*4 I,J
*
*     save variables from previous time step and preprocess
*
      DO 20 I = 1,IMAX
         AN(I) = A(I)
         AREAN(I) = AREA(I)
         CN(I) = C(I)
         GAMMAN(I) = GAMMA(I)
         QLATINN(I) = QLATIN(I)
         DO 10 J=1,NSOLUTE
            BTERMSN(I,J) = BTERMS(I,J)
            CLATINN(I,J) = CLATIN(I,J)
            BN(I,J) = TIMEB - BTERMSN(I,J)
     &                  - 0.5D0 * ALPHA(I)
     &                * (1.D0 - GAMMAN(I) / TWOPLUS(I,J))
            TGROUP(I,J) = 0.5D0 * ALPHA(I)*(4.D0+GAMMA(I)-GAMMAN(I))
     &                    / TWOPLUS(I,J)
            IGROUP(I,J) = .5D0*(QLATINN(I)*CLATINN(I,J)/AREAN(I) +
     &                     QLATIN(I)*CLATIN(I,J)/AREA(I))
 10      CONTINUE
 20   CONTINUE
*
*     these variables at 'N' only used for the first & last segments
*
      QN(1) = Q(1)
      QN(2) = Q(IMAX)
      AFACEN(1) = AFACE(1)
      AFACEN(2) = AFACE(IMAX)

      RETURN
      END
