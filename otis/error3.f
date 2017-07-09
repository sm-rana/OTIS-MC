************************************************************************
*
*     Subroutine ERROR3           Called by: INPUT1,3,6,7
*
*     print input error messages
*
************************************************************************
      SUBROUTINE ERROR3(NUMBER,PARAM)
*
*     dimensional parameters and logical devices
*
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'

      INTEGER*4 NUMBER,PARAM

 1000 FORMAT(//,2X,'Error: The Print Option must equal 1 or 2.',
     &//,2X,'Print Option specified (PRTOPT) : ',I5,//)
 2000 FORMAT(//,2X,'Error: The Interpolation Option must equal 0 or 1.',
     &//,2X,'Interpolation Option specified (IOPT) : ',I5,//)
 3000 FORMAT(//,2X,'Error: The Boundary Condition Option must equal 1,',
     &       ' 2 or 3.',//,2X,'Boundary Condition Option specified ',
     &       '(IBOUND) : ',I5,//)
 4000 FORMAT(//,2X,'Error: The Decay Option must equal 0 or 1.',
     &//,2X,'Decay Option specified (IDECAY) : ',I5,//)
 5000 FORMAT(//,2X,'Error: The Sorption Option must equal 0 or 1.',
     &//,2X,'Sorption Option specified (ISORB) : ',I5,//)

      IF (NUMBER .EQ. 1) THEN
         WRITE(*,1000) PARAM
         WRITE(LDECHO,1000) PARAM
      ELSEIF (NUMBER .EQ. 2) THEN
         WRITE(*,2000) PARAM
         WRITE(LDECHO,2000) PARAM
      ELSEIF (NUMBER .EQ. 3) THEN
         WRITE(*,3000) PARAM
         WRITE(LDECHO,3000) PARAM
      ELSEIF (NUMBER .EQ. 4) THEN
         WRITE(*,4000) PARAM
         WRITE(LDECHO,4000) PARAM
      ELSEIF (NUMBER .EQ. 5) THEN
         WRITE(*,5000) PARAM
         WRITE(LDECHO,5000) PARAM
      ENDIF

      WRITE(*,*) '  **** Fatal Input Error, See file echo.out ****'
      STOP

      END


