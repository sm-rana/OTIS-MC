************************************************************************
*
*     Subroutine INPUTQ              Called by: INPUT
*
*     input the following parameters:
*
*     QSTEP, QSTART, QLATIN, QLATOUT, AREA, CLATIN
*
************************************************************************
      SUBROUTINE INPUTQ(QSTEP,QSTART,QLATIN,QLATOUT,CLATIN,AREA,NREACH,
     &                  LASTSEG,NSOLUTE,FLOWLOC,QVALUE,AVALUE,NFLOW,
     &                  QINVAL,CLVAL)
*
*     dimensional parameters and logical devices
*
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
*     subroutine arguments
*
      INTEGER*4 NREACH,NSOLUTE,NFLOW,LASTSEG(0:MAXREACH)
      DOUBLE PRECISION QSTEP,QSTART,QLATIN(*),QLATOUT(*),AREA(*),
     &                 FLOWLOC(*),QVALUE(*),AVALUE(*),QINVAL(*),
     &                 CLATIN(MAXSEG,*),CLVAL(MAXFLOWLOC+1,*)
*
*     local variable
*
      CHARACTER*500 BUFFER
*
*     Input Format Statements
*
 1010 FORMAT(D13.5)
*
*     Output Format Statements
*
 2000 FORMAT(///,10X,'Flow and Area Data',/,10X,18('-'))
 2020 FORMAT(10X,'Flow Option                              : ',2X,A8)
 2040 FORMAT(10X,'Flow interval                       [hr] : ',1PE13.5)
 2060 FORMAT(//,10X,'Initial Flow and Area Values',/,10X,28('-'),/)
*
*     Read interval between changes in flow [hr].
*
      CALL GETLINE(LDFLOW,BUFFER)
      READ (BUFFER,1010) QSTEP
*
*     If the interval between changes in flow is zero, the steady flow
*     option has been selected.  Otherwise the flow is unsteady.  Call
*     the appropriate routine to read the flowrate at upstream boundary,
*     lateral inflows, lateral outflows, cross sectional areas, and 
*     inflow concentrations.  Also set QSTART equal to QVALUE(1) for use
*     in INPUT5.
*
      WRITE(LDECHO,2000)
      IF (QSTEP .EQ. 0.D0) THEN
         WRITE(LDECHO,2020) 'Steady  '
         WRITE(LDECHO,2040) QSTEP 
         WRITE(LDECHO,2060)
         CALL QSTEADY(QSTART,QLATIN,QLATOUT,CLATIN,AREA,NREACH,LASTSEG,
     &                NSOLUTE)
      ELSE
         WRITE(LDECHO,2020) 'Unsteady'
         WRITE(LDECHO,2040) QSTEP
         WRITE(LDECHO,2060)
         CALL QUNSTEADY(QINVAL,CLVAL,NSOLUTE,NFLOW,FLOWLOC,QVALUE,
     &                  AVALUE)
         QSTART = QVALUE(1)
      ENDIF

      RETURN
      END

