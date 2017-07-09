************************************************************************
*
*     Subroutine READQ           Called by: QCHANGE, QUNSTEADY, SETUP2
*
*     input the following time-variable parameters:
*
*     QINVAL, CLVAL, QVALUE, AVALUE
*
************************************************************************
      SUBROUTINE READQ(QINVAL,CLVAL,QVALUE,AVALUE,NSOLUTE,NFLOW)
*
*     dimensional parameters and logical devices
*
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
*     subroutine arguments
*
      INTEGER*4 NSOLUTE,NFLOW
      DOUBLE PRECISION QINVAL(*),QVALUE(*),AVALUE(*),
     &                 CLVAL(MAXFLOWLOC+1,*)
*
*     local variables
*
      INTEGER*4 I,J
      CHARACTER*500 BUFFER
*
*     input format statements
*
 1000 FORMAT(70D13.5)
*
*     read the lateral flows, main-channel flows and areas at the
*     various flow locations.
*
      CALL GETLINE(LDFLOW,BUFFER)
      READ(BUFFER,1000) (QINVAL(I),I=1,NFLOW)
      CALL GETLINE(LDFLOW,BUFFER)
      READ(BUFFER,1000) (QVALUE(I),I=1,NFLOW)
      CALL GETLINE(LDFLOW,BUFFER)
      READ(BUFFER,1000) (AVALUE(I),I=1,NFLOW)
      QINVAL(NFLOW+1) = QINVAL(NFLOW)
*
*     Read lateral inflow concentrations at each flow location
*
      DO 10 J=1,NSOLUTE
         CALL GETLINE(LDFLOW,BUFFER)
         READ(BUFFER,1000) (CLVAL(I,J),I=1,NFLOW)
         CLVAL(NFLOW+1,J) = CLVAL(NFLOW,J)
 10   CONTINUE

      RETURN
      END
