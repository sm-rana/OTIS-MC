************************************************************************
*
*     Subroutine OPENFILE2            Called by: MAININIT,MAINRUN
*
*     Open a simulation output file for each solute.  If sorption is
*     modeled, also open a sorption output file for each solute.
*
************************************************************************
      SUBROUTINE OPENFILE2(NSOLUTE,ISORB)
*
*     dimensional parameters and logical devices
*
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
*     subroutine arguments
*
      INTEGER*4 NSOLUTE,ISORB
*
*     local variables
*
      INTEGER*4 I
      CHARACTER*40 OUTFILE
      CHARACTER*500 BUFFER

 100  FORMAT(//,10X,'Output Files for Solutes',/,10X,24('-'),/,
     &       10X,'Solute #',12X,'File Name',/,10X,33('-'))
 200  FORMAT(//,10X,'Output Files for Sorbed Concentrations',/,10X,
     &       24('-'),/,10X,'Solute #',12X,'File Name',/,10X,47('-'))
 300  FORMAT(A40)
 400  FORMAT(12X,I5,13X,A40)
*
*     open solute output files
*
      WRITE(LDECHO,100)

      DO 10 I = 1,NSOLUTE
         CALL GETLINE(LDCTRL,BUFFER)
         READ (BUFFER,300) OUTFILE
         OPEN (UNIT=LDFILES(I),FILE=OUTFILE)
         WRITE (LDECHO,400) I,OUTFILE
 10   CONTINUE
*
*     open sorption output files
*
       IF (ISORB .EQ. 1) THEN
          WRITE(LDECHO,200)
          DO 20 I = 1,NSOLUTE
             CALL GETLINE(LDCTRL,BUFFER)
             READ (BUFFER,300) OUTFILE
             OPEN (UNIT=LDSORB(I),FILE=OUTFILE)
             WRITE (LDECHO,400) I,OUTFILE
 20      CONTINUE
      ENDIF

      RETURN
      END

