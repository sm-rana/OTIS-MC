************************************************************************
*
*     Subroutine FLOWINIT            Called by: MAININIT
*
*     initialize the flow rates for the case of steady flow
*
************************************************************************
      SUBROUTINE FLOWINIT(IMAX,DELTAX,QSTART,Q,QLATIN,QLATOUT)
*
*     subroutine arguments
*
      INTEGER*4 IMAX
      DOUBLE PRECISION QSTART,DELTAX(*),Q(*),QLATIN(*),QLATOUT(*)
*
*     local variables
*
      INTEGER*4 I
*
*     initialize flow rates
*
      Q(1) = QSTART + 0.5D0 * DELTAX(1) * (QLATIN(1) - QLATOUT(1))

      DO 10 I = 2,IMAX
         Q(I) = Q(I-1) + .5D0 *
     &                   ( (QLATIN(I-1) - QLATOUT(I-1))*DELTAX(I-1) +
     &                     (QLATIN(I) - QLATOUT(I))*DELTAX(I) )
 10   CONTINUE

      RETURN
      END
