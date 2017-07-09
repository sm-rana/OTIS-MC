************************************************************************
*
*     Subroutine ERROR               Called by: INPUT2,3,6,7, QUNSTEADY
*
*     print input error messages
*
************************************************************************
      SUBROUTINE ERROR(NUMBER,PARAM1,PARAM2)

      INCLUDE 'fmodules.inc'
*
*     logical devices
*
      INCLUDE 'lda.inc'
*
      INTEGER*4 NUMBER,PARAM1,PARAM2

 1000 FORMAT(//,2X,
     &'Error: The number of reaches specified exceeds the maximum.',
     &//,2X,
     &'Decrease the number of reaches or increase the maximum.',/,2X
     &,'To increase the allowable number of reaches, change MAXREACH'
     &,/,2X,'and recompile.',//,10X,
     &'Number of Reaches Specified       : ',I5,/,10X,
     &'Maximum number allowed (MAXREACH) : ',I5,//)
 1010 FORMAT(//,2X,
     &'Error: The total number of segments exceeds the maximum.',//,2X,
     &'Decrease the number of segments or increase the maximum.  To',/,
     &2X,'increase the allowable number of segments, change MAXSEG ',
     &'and recompile.',//,10X,'Number of Segments Specified    : ',I5,/,
     &10X,'Maximum number allowed (MAXSEG) : ',I5,//)
 1020 FORMAT(//,2X,
     &'Error: The number of print locations specified exceeds the ',
     &'maximum.',//,2X,
     &'Decrease the number of print locations or increase the maximum',/
     &,2X,'To increase the allowable number of print locations, change'
     &,/,2X,'MAXPRINT and recompile.',//,10X,
     &'Number of Print Locations Specified : ',I5,/,10X,
     &'Maximum number allowed (MAXPRINT)   : ',I5,//)
 1030 FORMAT(//,2X,
     &'Error: The number of solutes specified exceeds the ',
     &'maximum.',//,2X,
     &'Decrease the number of solutes or increase the maximum.',/
     &,2X,'To increase the allowable number of solutes, change'
     &,/,2X,'MAXSOLUTE and recompile.',//,10X,
     &'Number of Solutes Specified        : ',I5,/,10X,
     &'Maximum number allowed (MAXSOLUTE) : ',I5,//)
 1040 FORMAT(//,2X,
     &'Error: The number of upstream boundary conditions specified',/,
     &9X,'exceeds the maximum.',//,2X,
     &'Decrease the number of boundary conditions or increase the',/,2X,
     &'maximum.  To increase the allowable number of conditions,',/,2X,
     &'change MAXBOUND and recompile.',//,10X,
     &'Number of Conditions Specified  : ',I5,/,10X,
     &'Maximum number allowed (MAXBOUND) : ',I5,//)
 1050 FORMAT(//,2X,
     &'Error: The number of flow locations specified exceeds the ',
     &'maximum.',//,2X,
     &'Decrease the number of flow locations or increase the maximum.',/
     &,2X,'To increase the allowable number of flow locations, change'
     &,/,2X,'MAXFLOWLOC and recompile.',//,10X,
     &'Number of Flow Locations Specified  : ',I5,/,10X,
     &'Maximum number allowed (MAXFLOWLOC) : ',I5,//)

      IF (NUMBER .EQ. 1) THEN
         WRITE(*,1000) PARAM1,PARAM2
         WRITE(LDECHO,1000) PARAM1,PARAM2
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

      WRITE(*,*) '  **** Fatal Input Error, See file echo.out ****'
      STOP

      END

