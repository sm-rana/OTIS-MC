************************************************************************
*
*     Subroutine DISTANCE            Called by: WEIGHTS
*
*     compute the distance at the center of each segment
*
************************************************************************
      SUBROUTINE DISTANCE(IMAX,DELTAX,XSTART,DIST)
*
*     local variables
*
      INTEGER*4 I
*
*     subroutine arguments
*
      INTEGER*4 IMAX
      DOUBLE PRECISION XSTART,DELTAX(*),DIST(*)
*
*     compute the distances at the center of each segment
*
      DIST(1) = XSTART + 0.5D0 * DELTAX(1)
      DO 10 I = 2,IMAX
         DIST(I) = DIST(I-1) + .5D0 * (DELTAX(I-1) + DELTAX(I))
 10   CONTINUE

      RETURN
      END
