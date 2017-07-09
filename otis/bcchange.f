************************************************************************
*
*     Subroutine BCCHANGE             Called by: UPDATE
*
*     Update the upstream boundary condition and set the time of the
*     next boundary condition change.
*
************************************************************************
      SUBROUTINE BCCHANGE(NSOLUTE,JBOUND,IBOUND,TSTOP,BCSTOP,TSTEP,TIME,
     &                    USTIME,USCONC,USBC)
*
*     dimensional parameters
*
      INCLUDE 'fmodules.inc'
*
*     subroutine arguments
*
      INTEGER*4 NSOLUTE,JBOUND,IBOUND
      DOUBLE PRECISION TSTOP,BCSTOP,TSTEP,TIME,USTIME(*)
      DOUBLE PRECISION USCONC(MAXBOUND,*),USBC(MAXBOUND,*)
*
*     function definition
*
      DOUBLE PRECISION YINTERP
*
*     local variables
*
      INTEGER*4 J
*
*     initialize the upstream boundary condition variables.  Note that
*     for a step boundary condition (IBOUND = 1 or 2), the boundary
*     condition changes at user-specified times indicated by USTIME.
*     For a continous boundary condition (IBOUND =3), the boundary
*     condition changes every timestep.
*
*     For a step boundary condition, simply increment the boundary
*     condition counter (JBOUND) and set the time of the next boundary
*     condition change (BCSTOP). 
*
*     For a continuous boundary condition, increment JBOUND such that
*     the endpoints for interpolation are set correctly, set the
*     boundary condition for the current time step by interpolating
*     between user-specified data, and set BCSTOP.
*
      IF (IBOUND .EQ. 1 .OR. IBOUND .EQ. 2) THEN
         JBOUND = JBOUND + 1
         BCSTOP = USTIME(JBOUND+1)
      ELSE
 10      IF (TIME .GT. USTIME(JBOUND+1)) THEN
            JBOUND = JBOUND + 1
            GOTO 10
         ENDIF
         DO 20 J = 1,NSOLUTE
            USCONC(JBOUND,J) = YINTERP(TIME,USTIME(JBOUND),
     &                                 USTIME(JBOUND+1),USBC(JBOUND,J),
     &                                 USBC(JBOUND+1,J))
 20      CONTINUE

         BCSTOP = TSTOP + TSTEP
      ENDIF
  
      RETURN
      END

