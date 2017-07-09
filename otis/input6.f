************************************************************************
*
*     Subroutine INPUT6               Called by: INPUT
*
*     read the print parameters
*
************************************************************************
      SUBROUTINE INPUT6(NPRINT,PRTLOC,IMAX,DELTAX,XSTART,PINDEX,WT,DIST)
*
*     dimensional parameters and logical devices
*
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
*     subroutine arguments
*
      INTEGER*4 IMAX,NPRINT,PINDEX(*)
      DOUBLE PRECISION XSTART
      DOUBLE PRECISION PRTLOC(*),DELTAX(*),WT(*),DIST(*)
*
*     local variables
*
      INTEGER*4 I,IOPT
      DOUBLE PRECISION RPRTLOC(MAXPRINT)
      CHARACTER*500 BUFFER
*
*     Input Format Statements
*
 1000 FORMAT(2I5)
 1040 FORMAT(F13.2)
*
*     Output Format Statements
*
 2000 FORMAT(///,10X,'Number of Print Locations : ',I3,/,10X,
     &       'Interpolation Option      : ',I3,/)
 2020 FORMAT(10X,'Printing Information',/,10X,20('-'),/,10X,'Requested',
     &       9X,'Effective'/,10X,'Print',13X,'Print',/,10X,
     &       'Location [L]',6X,'Location [L]',/,10X,30('-'))
 2040 FORMAT(10X,F12.3,6X,F12.3)
*
*     Read the number of print locations and the interpolation option.
*     Test to see if the number of locations exceeds the maximum and for
*     a valid interpolation option.
*
      CALL GETLINE(LDPARAM,BUFFER)
      READ (BUFFER,1000) NPRINT,IOPT
      WRITE (LDECHO,2000) NPRINT,IOPT
      !SM:  NPRINT = 1
      IF (NPRINT .GT. MAXPRINT) CALL ERROR(3,NPRINT,MAXPRINT)
      IF (IOPT .NE. 0 .AND. IOPT .NE. 1) CALL ERROR3(2,IOPT)
*
*     Read the requested print locations, RPRTLOC [L].  Determine
*     the interpolation weights for printing (S/R Weights).  If
*     interpolation is not requested (IOPT=0), the weights are set to
*     zero and the print location used (PRTLOC) will be equal to the
*     distance of the nearest segment centroid that is upstream from
*     RPRTLOC.
*
      DO 10 I = 1,NPRINT
         CALL GETLINE(LDPARAM,BUFFER)
         READ (BUFFER,1040) RPRTLOC(I)
 10   CONTINUE

      CALL WEIGHTS(IMAX,DELTAX,XSTART,PINDEX,WT,PRTLOC,NPRINT,RPRTLOC,
     &             IOPT,DIST)

      WRITE(LDECHO,2020)

      DO 20 I=1,NPRINT
         WRITE(LDECHO,2040) RPRTLOC(I),PRTLOC(I)
 20   CONTINUE

      RETURN
      END
