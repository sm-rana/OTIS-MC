************************************************************************
*     
*     Subroutine QAINIT            Called by: QCHANGE, MAININIT, SETUP2
*
*     initialize the flow rates, lateral flows and areas for the case of
*     unsteady flow.
* 
*     Flow and area values for each segment are interpolated based on
*     the data that is available at the user-specified flow locations.
*
*     Lateral flow rates and inflow concentrations are set equal to
*     the values at the nearest downstream flow location.  Note that
*     QINDEX(I) denotes the nearest downstream flow location, DSDIST(I)
*     is the distance from the center of the segment to the nearest
*     downstream flow location (FLOWLOC(QINDEX(I))), and USDIST(I) is
*     the distance to the nearest upstream flow location
*     (FLOWLOC(QINDEX(I)-1)).  USDIST and DSDIST are normalized to the
*     segment length.
*
*     Note that the interpolation weights and indices are precomputed
*     in the subroutine QWEIGHTS.
*
************************************************************************
      SUBROUTINE QAINIT(IMAX,Q,AREA,QVALUE,AVALUE,QWT,QINDEX,QINVAL,
     &                  CLVAL,QLATIN,CLATIN,DSDIST,USDIST,NSOLUTE)
      INCLUDE 'fmodules.inc'
*
*     local variables
*
      INTEGER*4 I,J
*
*     subroutine arguments
*
      INTEGER*4 IMAX,NSOLUTE,QINDEX(*)
      DOUBLE PRECISION Q(*),AREA(*),QVALUE(*),AVALUE(*),QWT(*),
     &                 QINVAL(*),QLATIN(*),DSDIST(*),USDIST(*),
     &                 CLATIN(MAXSEG,*),CLVAL(MAXFLOWLOC+1,*)
*
*     function declaration
*
      DOUBLE PRECISION WTAVG
*
*     initialize flow rates and areas
*
       DO 10 I = 1,IMAX
         J = QINDEX(I)
         Q(I) = QVALUE(J-1) + QWT(I) * (QVALUE(J) - QVALUE(J-1))
         AREA(I) = AVALUE(J-1) + QWT(I) * (AVALUE(J) - AVALUE(J-1))
 10   CONTINUE
*
*     set lateral flows and inflow concentrations equal to the values at
*     the nearest downstream flow location.  If the flow location is w/i
*     the current segment, set the flows and concentrations equal to a
*     weighted average of the closest flow location and the 2nd closest.
*     In this way QLATIN and CLATIN are spatially constant between flow
*     locations.
*
      DO 50 I = 1,IMAX
         IF (USDIST(I) .LT. 0.5D0) THEN
*
*          a flow location is w/i the upstream half of the segment
* 
           QLATIN(I) = WTAVG(USDIST(I),QINVAL(QINDEX(I)-1),
     &                                  QINVAL(QINDEX(I)))
            DO 20 J = 1,NSOLUTE
               CLATIN(I,J) = WTAVG(USDIST(I),CLVAL(QINDEX(I)-1,J),
     &                                       CLVAL(QINDEX(I),J))
 20         CONTINUE
         ELSEIF (DSDIST(I) .LT. 0.5D0) THEN
*
*           a flow location is w/i the downstream half of the segment
*
            QLATIN(I) = WTAVG(DSDIST(I),QINVAL(QINDEX(I)+1),
     &                                  QINVAL(QINDEX(I)))
            DO 30 J = 1,NSOLUTE
               CLATIN(I,J) = WTAVG(DSDIST(I),CLVAL(QINDEX(I)+1,J),
     &                                       CLVAL(QINDEX(I),J))
 30         CONTINUE
         ELSE
            QLATIN(I) = QINVAL(QINDEX(I))
            DO 40 J = 1,NSOLUTE
               CLATIN(I,J) = CLVAL(QINDEX(I),J)
 40         CONTINUE

         ENDIF
 50   CONTINUE

      RETURN
      END
