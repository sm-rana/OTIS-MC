************************************************************************
*
*     Subroutine SETTYPE               Called by: MAININIT
*
*     determine the type of simulation to conduct so that the proper
*     preprocessing and finite difference routines are called.  The
*     simulation will be either:
*
*     1) 'Steady-State' with conservative or reactive solutes
*     2) dynamic, with 'Conservative' solutes
*     3) dynamic, with 'Reactive' solutes
*
************************************************************************
      SUBROUTINE SETTYPE(IDECAY,ISORB,TSTEP,CHEM)
*
*     dimensional parameters and logical devices
*
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
*     subroutine arguments
*
      INTEGER*4 IDECAY,ISORB
      DOUBLE PRECISION TSTEP
      CHARACTER*(*) CHEM
*
*     set chemistry indicator, CHEM, so that the proper preprocessing
*     and finite difference routines are called.
*     SM: my work is not Steady-State simulation --> for my work it is 'Conservative'
*
      IF (TSTEP .LT. 1.D-10) THEN
         CHEM = 'Steady-State'
      ELSEIF (IDECAY .EQ. 0 .AND. ISORB .EQ. 0) THEN
         CHEM = 'Conservative'
      ELSE
         CHEM = 'Reactive'
      ENDIF

      RETURN
      END
