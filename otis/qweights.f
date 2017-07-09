************************************************************************
*
*     Subroutine QWEIGHTS            Called by: MAININIT
*
*     determine the weights that are used to interpolate flow and area
*     values.  The linear interpolation formula used is of the form:
*
*     Q(I) = QVALUE(J-1) + WT * [ QVALUE(J) - QVALUE(J-1) ]
*
*          where
*                        [ DIST(I) - FLOWLOC(J-1) ] 
*                QWT = ------------------------------
*                       [ FLOWLOC(I) - FLOWLOC(I-1) ]
*
*     Also compute the distance from the center of each segment to the
*     nearest downstream flow location, DSDIST, and to the nearest
*     upstream flow location, USDIST.  This distances are normalized
*     using the segment length.
*
*     Variables
*     ---------
*     I         segment
*     J         flow location
*     DIST(I)   distance at the center of segment I
*
************************************************************************
      SUBROUTINE QWEIGHTS(IMAX,DELTAX,XSTART,QINDEX,QWT,FLOWLOC,NFLOW,
     &                    DSDIST,USDIST,DIST)
      INCLUDE 'fmodules.inc'
*
*     local variables
*
      INTEGER*4 I,J
*
*     subroutine arguments
*
      INTEGER*4 IMAX,NFLOW,QINDEX(*)
      DOUBLE PRECISION XSTART,DELTAX(*),QWT(*),FLOWLOC(*),DSDIST(*),
     &                 USDIST(*),DIST(*)
*
*     check to make sure that the first flow location is placed at the
*     upstream boundary.  Also check to make sure that the last flow
*     location is at (or downstream from) the last segment.
*
      IF (FLOWLOC(1) .NE. XSTART) CALL ERROR2(4,FLOWLOC(1),XSTART)
      IF (FLOWLOC(NFLOW) .LT. DIST(IMAX))
     &   CALL ERROR2(5,FLOWLOC(NFLOW),DIST(IMAX))
*
*     determine the interpolation weights.  If the distance at the
*     center of the current segment is greater than the current flow
*     location, go to the next flow location.  For each segment, set the 
*     interpolation weight (QWT) and the index corresponding to the 
*     appropriate flow location (QINDEX).  Also compute DSDIST and 
*     USDIST, the normalized distances to the nearest flow locations.
*
      J = 2

      DO 20 I = 1,IMAX
 10      IF (DIST(I) .GE. FLOWLOC(J)) THEN
            J = J + 1
            GOTO 10
         ENDIF
         QWT(I) = (DIST(I) - FLOWLOC(J-1))/(FLOWLOC(J)-FLOWLOC(J-1))
         QINDEX(I) = J
         DSDIST(I) = (FLOWLOC(J) - DIST(I)) / DELTAX(I)
         USDIST(I) = (DIST(I) - FLOWLOC(J-1)) / DELTAX(I)
 20   CONTINUE

      RETURN
      END
