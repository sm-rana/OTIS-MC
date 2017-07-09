************************************************************************
*
*     Subroutine SETSTOP            Called by: PREPROC, UPDATE, SETUP3
*
*     determine the endpoint of the integration, TSTOP.  TSTOP is
*     generally equal to the time at which the boundary condition or
*     flow variables change.
*
************************************************************************
      SUBROUTINE SETSTOP(QSTOP,TSTOP,STOPTYPE,BCSTOP,IBOUND)
*
*     subroutine arguments
*
      INTEGER*4 IBOUND
      DOUBLE PRECISION QSTOP,TSTOP,BCSTOP
      CHARACTER*(*) STOPTYPE
*
*     For a continuous boundary condition (IBOUND=3), force a boundary
*     condition change for every change in the flow.
*       ! masud: BCSTOP == TSTART, QSTOP = 999e99
      IF (IBOUND .EQ. 3 .AND. QSTOP .LT. BCSTOP) BCSTOP = QSTOP
*
*     If the boundary condition and the flow change at the same time,
*     set the stopping point and type accordingly.  If the flow changes
*     before the boundary condition, the stopping point is set to the
*     time at which the flow changes.  If the boundary condition changes
*     first, set the stopping point equal to the time at which the
*     boundary condition changes.
*
      IF (ABS(QSTOP-BCSTOP).LT. 1.D-7) THEN
         TSTOP = QSTOP
         STOPTYPE = 'Both'
      ELSEIF (QSTOP.LT.BCSTOP) THEN
         TSTOP = QSTOP
         STOPTYPE = 'QChange'
      ELSEIF (BCSTOP.LT.QSTOP) THEN
         TSTOP = BCSTOP
         STOPTYPE = 'BCondition'
      ENDIF

      RETURN
      END
