************************************************************************
*
*     Subroutine WEIGHTS            Called by: INPUT6
*
*     Determine the weights that are used to interpolate linearly
*     between segment centroids.  The linear interpolation formula used
*     is of the form:
*
*     Conc@PRTLOC = Conc@I + WT * [ Conc@I+1 - Conc@I ]
*
*          where
*                 WT = [ PRTLOC - DIST(I) ] / [ DIST(I+1) - DIST(I) ]
*
*     If interpolation is not requested (IOPT=0), just find the nearest
*     upstream segment and set the print location accordingly.
*
*     Variables
*     ---------
*     RPRTLOC(J) requested print location - distance at which results
*                are desired
*     PRTLOC(J)  print location - distance for which results are printed
*     WT(J)      interpolation weight corresponding to RPRTLOC(J)
*     PINDEX     index corresponding to RPRTLOC(J) denoting segment I
*     DIST(I)    distance at the center of segment I
*
************************************************************************
      SUBROUTINE WEIGHTS(IMAX,DELTAX,XSTART,PINDEX,WT,PRTLOC,NPRINT,
     &                   RPRTLOC,IOPT,DIST)

      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
*     subroutine arguments
*
      INTEGER*4 IMAX,IOPT,NPRINT,PINDEX(*)
      DOUBLE PRECISION XSTART,DELTAX(*),WT(*),PRTLOC(*),RPRTLOC(*),
     &                 DIST(*)
*
*     local variables
*
      INTEGER*4 I,J
*
*     format statement
*
 100  FORMAT(//,10X,'WARNING: Print Location (',F14.5,
     &       ')',/,10X,'is less than distance of the first segment.',//,
     &       10X,'Print location set to ',F14.5,///)
*
*     compute the distances at the center of each segment
*
      CALL DISTANCE(IMAX,DELTAX,XSTART,DIST)
*
*     determine the interpolation weights & conduct error tests
*
      J = 1

      DO 20 I = 1,NPRINT
*
*        Error Checks
*
*        1) Test to make sure the print location does not exceed the
*           distance at the end of the stream.
*        2) Make sure the print location is not in the upstream half
*           of the first segment. If it is, set it equal to the distance
*           at the center of the first segment.
*
         IF (RPRTLOC(I) .GT. DIST(IMAX))
     &      CALL ERROR2(2,RPRTLOC(I),DIST(IMAX))
         IF (RPRTLOC(I) .LT. DIST(1)) THEN
            WRITE(LDECHO,100) RPRTLOC(I),DIST(1)
            RPRTLOC(I) = DIST(1)
         ENDIF
*
*        Sequentially search through the vector of distances until
*        the print location falls between two of the vector elements
*        (i.e. the print location is between the centers of two adjacent
*        segments).  Note that the test interval is inclusive
*        [.GE. DIST(J) and .LE. DIST(J+1)] to avoid searching past the
*        last segment. Set PINDEX equal to the number of the upstream
*        segment.  
*
*        If interpolation is requested (IOPT=1): compute the
*        interpolation weight and set the print location equal to the
*        requested print location.
*
*        If interpolation is not requested (IOPT=0):
*          1) If RPRTLOC=DIST(J+1), set the weight equal to 1 and set
*             the print location equal to the distance at the center of
*             segment J+1
*          2) else, set the weight equal to zero and set the print
*             location equal to the distance at the center of the
*             upstream segment.
*
 10      IF ((RPRTLOC(I).GE.DIST(J)).AND.(RPRTLOC(I).LE.DIST(J+1))) THEN
            PINDEX(I) = J
            IF (IOPT .EQ. 1) THEN
               WT(I) = (RPRTLOC(I) - DIST(J)) / (DIST(J+1) - DIST(J))
               PRTLOC(I) = RPRTLOC(I)
            ELSE
               IF (RPRTLOC(I) .EQ. DIST(J+1)) THEN
                  WT(I) = 1.D0
                  PRTLOC(I) = RPRTLOC(I)
               ELSE
                  WT(I) = 0.D0
                  PRTLOC(I) = DIST(J)
               ENDIF
            ENDIF
            J = 1
         ELSE
            J = J + 1
            GOTO 10
         ENDIF

 20   CONTINUE

      RETURN
      END
