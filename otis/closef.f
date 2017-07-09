************************************************************************
*
*     Subroutine CLOSEF               Called by: MAIN
*
*     close the flow input file, the echo output file, the solute output
*     files and the sorption output files.
*
************************************************************************
      SUBROUTINE CLOSEF(NSOLUTE,ISORB)
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
*
*     close the flow input file
*
      CLOSE(UNIT=LDFLOW)
*
*     close the echo output file
*
      CLOSE (UNIT=LDECHO)
*
*     close the solute output files
*
      DO 10 I = 1,NSOLUTE
         CLOSE(UNIT=LDFILES(I))
 10   CONTINUE
*
*     close the sorption output files, if they exist.
*
      IF (ISORB .EQ. 1) THEN
         DO 20 I = 1,NSOLUTE
            CLOSE (UNIT=LDSORB(I))
 20      CONTINUE
      ENDIF

      RETURN
      END




