************************************************************************
*
*     Subroutine QCHANGE             Called by: UPDATE
*
*     Read in the new flow variables, initialize the vectors, determine
*     the next time the flow changes, and set the concentration at the
*     upstream boundary.
*
************************************************************************
      SUBROUTINE QCHANGE(IMAX,Q,AREA,QLATIN,CLATIN,NSOLUTE,QVALUE,
     &                   AVALUE,QWT,QINDEX,NFLOW,JBOUND,QSTOP,QSTEP,
     &                   TSTOP,USCONC,QINVAL,CLVAL,DSDIST,USDIST,IBOUND,
     &                   USBC,NBOUND)
      INCLUDE 'fmodules.inc'
*
*     subroutine arguments
*
      INTEGER*4 IMAX,NSOLUTE,NFLOW,NBOUND,JBOUND,IBOUND,QINDEX(*)
      DOUBLE PRECISION QSTOP,QSTEP,TSTOP
      DOUBLE PRECISION Q(*),AREA(*),QLATIN(*),QVALUE(*),AVALUE(*),
     &                 QWT(*),QINVAL(*),DSDIST(*),USDIST(*)
      DOUBLE PRECISION CLATIN(MAXSEG,*),USCONC(MAXBOUND,*),           
     &                 CLVAL(MAXFLOWLOC+1,*),
     &                 USBC(MAXBOUND,*)
*
*     local variables
*
      INTEGER*4 I,J
*
*     Read in the new flow variables and initialize the vectors.
*
      CALL READQ(QINVAL,CLVAL,QVALUE,AVALUE,NSOLUTE,NFLOW)
      CALL QAINIT(IMAX,Q,AREA,QVALUE,AVALUE,QWT,QINDEX,QINVAL,CLVAL,
     &            QLATIN,CLATIN,DSDIST,USDIST,NSOLUTE)
*
*     determine the next time the flow changes
*
      QSTOP = TSTOP + QSTEP
*
*     if a flux boundary condition is in effect, update the
*     concentration at the upstream boundary based on the new flowrate
*
      IF (IBOUND .EQ. 2) THEN
         DO 10 J = 1,NSOLUTE
            DO 20 I = JBOUND,NBOUND
               USCONC(I,J) = USBC(I,J) / QVALUE(1)
 20         CONTINUE
 10      CONTINUE 
      ENDIF

      RETURN
      END
