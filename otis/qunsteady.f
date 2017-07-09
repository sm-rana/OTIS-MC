************************************************************************
*
*     Subroutine QUNSTEADY              Called by: INPUTQ
*
*     input the following parameters:
*
*     QINVAL, CLVAL, NFLOW, FLOWLOC, QVALUE, AVALUE
*
************************************************************************
      SUBROUTINE QUNSTEADY(QINVAL,CLVAL,NSOLUTE,NFLOW,FLOWLOC,QVALUE,
     &                     AVALUE)
*
*     dimensional parameters and logical devices
*
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
*     subroutine arguments
*
      INTEGER*4 NSOLUTE,NFLOW
      DOUBLE PRECISION QINVAL(*),FLOWLOC(*),QVALUE(*),AVALUE(*),
     &                 CLVAL(MAXFLOWLOC+1,*)
*
*     local variables
*
      INTEGER*4 I,J
      CHARACTER*500 BUFFER
*
*     Input Format Statements
*
 1000 FORMAT(I5)
 1020 FORMAT(F13.2)
*
*     Output Format Statements
*
 2000 FORMAT(10X,'Number of Flow Locations : ',I3,/)
 2020 FORMAT(10X,'Locations at which Q is Specified',/,10X,33('-'),
     &       /,10X,
     &       'Flow           Channel          Lateral          Area',/,
     &       10X,'Location       Flow             Inflow',/,10X,
     &       '   [L]         [L^3/sec]        [L^3/s-L]        [L^2]',
     &       /,10X,63('-'))
 2040 FORMAT(8X,F12.3,1P,3(5X,E12.6))
 2060 FORMAT(//,10X,'Initial Lateral Inflow Concentrations, Solute #',
     &       I3,/,10X,51('-'))
 2080 FORMAT(10X,'Flow             Inflow',/,10X,
     &       'Location         Concent.',/,10X,
     &       '  [L]            [mass/L^3]',/,10X,30('-'))
 2100 FORMAT(8X,F12.3,7X,1PE11.5)
*
*     Read the number of flow locations.  Test to see if
*     number of locations exceeds the maximum.
*
      CALL GETLINE(LDFLOW,BUFFER)
      READ (BUFFER,1000) NFLOW
      WRITE (LDECHO,2000) NFLOW
      IF (NFLOW .GT. MAXFLOWLOC) CALL ERROR(6,NFLOW,MAXFLOWLOC)
*
*     Read the flow locations [L].  Note that the flow locations must be
*     specified in ascending order.
*
      DO 10 I = 1,NFLOW
         CALL GETLINE(LDFLOW,BUFFER)
         READ (BUFFER,1020) FLOWLOC(I)
 10   CONTINUE
      DO 20 I = 2,NFLOW
         IF (FLOWLOC(I) .LT. FLOWLOC(I-1))
     &      CALL ERROR2(3,FLOWLOC(I-1),FLOWLOC(I))
 20   CONTINUE
*
*     Read lateral inflow, lateral outflow, lateral inflow
*     concentrations, flow and area values for each flow location.
*
      CALL READQ(QINVAL,CLVAL,QVALUE,AVALUE,NSOLUTE,NFLOW)
*
*     echo flow, lateral inflow and area values at the flow locations
*
      WRITE(LDECHO,2020)
      DO 30 I = 1,NFLOW
         WRITE(LDECHO,2040) FLOWLOC(I),QVALUE(I),QINVAL(I),AVALUE(I)
 30   CONTINUE
*
*     echo lateral inflow concentrations
*
      DO 50 J = 1,NSOLUTE
         WRITE(LDECHO,2060) J
         WRITE(LDECHO,2080)
         DO 40 I = 1,NFLOW
            WRITE(LDECHO,2100) FLOWLOC(I),CLVAL(I,J)
 40      CONTINUE
 50   CONTINUE

      RETURN
      END
