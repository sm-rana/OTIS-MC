************************************************************************
*
*     Subroutine ERROR2               Called by: INPUT2, QUNSTEADY
*                                                QWEIGHTS, WEIGHTS
*     print error messages
*
************************************************************************
      SUBROUTINE ERROR2(NUMBER,PARAM1,PARAM2)
*
*     dimensional parameters and logical devices
*
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
      INTEGER*4 NUMBER
      DOUBLE PRECISION PARAM1,PARAM2

 1000 FORMAT(//,2X,
     &'Error: The storage area specified must be greater than zero.',
     &//,2X,'Storage Area Specified: ',1PE16.8,//)
 1010 FORMAT(//,2X,'Error:',//,2X,
     &'The print location specified exceeds the length of the stream.',
     &//,2X,'Print Location Specified  : ',1PE16.8,
     & /,2X,'End of Stream             : ',1PE16.8,//)
 1020 FORMAT(//,2X,'Error:',//,2X,
     &'The flow locations must be entered in ascending order.',//,2X,
     &'Flow Location J-1 {FLOWLOC(J-1)}  : ',1PE16.8,/,2X,
     &'Flow Location J   {FLOWLOC(J)}    : ',1PE16.8,//)
 1030 FORMAT(//,2X,'Error:',//,2X,
     &'The first flow location must be equal to the distance at',/,2X,
     &'the upstream boundary.',//,2X,
     &'First Flow Location {FLOWLOC(1)}            : ',1PE16.8,/,2X,
     &'Distance at the Upstream Boundary {XSTART}  : ',1PE16.8,//)
 1040 FORMAT(//,2X,'Error:',//,2X,
     &'The last flow location must be at (or below) the ',
     &'downstream boundary.',//,2X,
     &'Last Flow Location {FLOWLOC()}       : ',1PE16.8,/,2X,
     &'Distance at the Downstream Boundary  : ',1PE16.8,//)
 1050 FORMAT(//,2X,'Error:',//,2X,
     &'For IBOUND=3, the time of the last boundary condition ',/,2X,
     &'must be greater than or equal to the simulation end time.',//,2X,
     &'Time of Last Boundary Condition {USTIME(NBOUND)} : ',1PE16.8,/,
     &2X,'Simulation End Time {TFINAL}                     : ',1PE16.8,
     &//)
 2000 FORMAT(2X,'**** Fatal Input Error, See file echo.out ****',//)
      IF (NUMBER .EQ. 1) THEN
         WRITE(*,1000) PARAM1
         WRITE(LDECHO,1000) PARAM1
      ELSEIF (NUMBER .EQ. 2) THEN
         WRITE(*,1010) PARAM1,PARAM2
         WRITE(LDECHO,1010) PARAM1,PARAM2
      ELSEIF (NUMBER .EQ. 3) THEN
         WRITE(*,1020) PARAM1,PARAM2
         WRITE(LDECHO,1020) PARAM1,PARAM2
      ELSEIF (NUMBER .EQ. 4) THEN
         WRITE(*,1030) PARAM1,PARAM2
         WRITE(LDECHO,1030) PARAM1,PARAM2
      ELSEIF (NUMBER .EQ. 5) THEN
         WRITE(*,1040) PARAM1,PARAM2
         WRITE(LDECHO,1040) PARAM1,PARAM2
      ELSEIF (NUMBER .EQ. 6) THEN
         WRITE(*,1050) PARAM1,PARAM2
         WRITE(LDECHO,1050) PARAM1,PARAM2
      ENDIF

      WRITE(*,2000)
      STOP

      END
