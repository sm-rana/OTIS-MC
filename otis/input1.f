************************************************************************
*
*     Subroutine INPUT1             Called by: INPUT
*
*     input the following parameters:
*
*     PRTOPT, TSTEP, PSTEP, TSTART, TFINAL, XSTART, DSBOUND
*
************************************************************************
      SUBROUTINE INPUT1(TSTEP,TSTART,TFINAL,IPRINT,XSTART,DSBOUND,
     &                  PRTOPT)
*
*     dimensional parameters and logical devices
*
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
*     local variables
*
      CHARACTER*80 TITLE
      DOUBLE PRECISION PSTEP
      CHARACTER*500 BUFFER
*
*     subroutine arguments
*
      INTEGER*4 IPRINT,PRTOPT
      DOUBLE PRECISION TSTEP,TSTART,TFINAL,XSTART,DSBOUND
*
*     Input Format Statements
*
 1010 FORMAT(A80)
 1020 FORMAT(I5)
 1030 FORMAT(D13.5)
*
*     Output Format Statements
*
 2000 FORMAT(///,10X,'Input Data',/,10X,10('-'),//,10X,A80,//,10X)
 2020 FORMAT(10X,'Print Option                              : ',I5)
 2040 FORMAT(10X,'Requested Print Interval              [hr]: ',1PE13.5)
 2060 FORMAT(10X,'Integration Time Step                 [hr]: ',1PE13.5)
 2080 FORMAT(10X,'Starting Time                         [hr]: ',1PE13.5)
 2100 FORMAT(10X,'Ending Time                           [hr]: ',1PE13.5)
 2120 FORMAT(10X,'Starting Distance                      [L]: ',1PE13.5)
 2140 FORMAT(10X,'Downstream Boundary Condition [mass/L^2-s]: ',1PE13.5)
 2160 FORMAT(10X,'Simulation Type  [dynamic or steady-state]:  ',A12)
*
*     read the simulation title
*
      CALL GETLINE(LDPARAM,BUFFER)
      READ (BUFFER,1010) TITLE
      WRITE (LDECHO,2000) TITLE
*
*     Read the print option
*
      CALL GETLINE(LDPARAM,BUFFER)
      READ (BUFFER,1020) PRTOPT
      WRITE (LDECHO,2020) PRTOPT
      IF ((PRTOPT .LT. 1).OR.(PRTOPT .GT. 2)) CALL ERROR3(1,PRTOPT)
*
*     read the time interval for the printing of results, PSTEP, and the
*     time step, TSTEP.
*
      CALL GETLINE(LDPARAM,BUFFER)
      READ (BUFFER,1030) PSTEP
      WRITE (LDECHO,2040) PSTEP
      CALL GETLINE(LDPARAM,BUFFER)
      READ (BUFFER,1030) TSTEP
      WRITE (LDECHO,2060) TSTEP
*
*     Compute IPRINT, the number of time steps corresponding to PSTEP.
*     Note that for steady-state runs (TSTEP=0), IPRINT is set to zero
*     as it is not used.
*
      IF (TSTEP .NE. 0.D0) THEN
         IPRINT = NINT(PSTEP/TSTEP)
      ELSE
         IPRINT = 0
      ENDIF
      IF (IPRINT .LT. 1) IPRINT = 1
*
*     read the starting time and the final time [hr]
*
      CALL GETLINE(LDPARAM,BUFFER)
      READ (BUFFER,1030) TSTART
      WRITE (LDECHO,2080) TSTART     
      CALL GETLINE(LDPARAM,BUFFER)
      READ (BUFFER,1030) TFINAL
      WRITE (LDECHO,2100) TFINAL
*
*     Read the distance at the upstream boundary [L] and the 
*     downstream boundary condition [mass/L^2-sec].
*
      CALL GETLINE(LDPARAM,BUFFER)
      READ (BUFFER,1030) XSTART
      CALL GETLINE(LDPARAM,BUFFER)
      READ (BUFFER,1030) DSBOUND
      WRITE (LDECHO,2120) XSTART
      WRITE (LDECHO,2140) DSBOUND
*
*     Determine the simulation type and echo.  Note that the user
*     selects the steady-state option by setting TSTEP to zero.
*
      IF (TSTEP .LT. 1.D-10) THEN
         WRITE(LDECHO,2160) 'Steady-State'
      ELSE
         WRITE(LDECHO,2160) 'Dynamic'
      ENDIF
      
      RETURN
      END
