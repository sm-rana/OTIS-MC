************************************************************************
*
*     Subroutine OPENIN               Called by: OPENFILE
*
*     open specified input file.  Print error message and terminate
*     execution if the specified file does not exist.
*
*     Subroutine Arguments
*     --------------------
*     LDNUM    - logical device unit number assigned to file
*     FILENAME - name of input file to OPEN
*
************************************************************************
      SUBROUTINE OPENIN(LDNUM,FILENAME)
*
*     dimensional parameters and logical devices
*
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
*     subroutine arguments
*
      CHARACTER FILENAME*(*)
      INTEGER*4 LDNUM
*
*     local variables
*
      LOGICAL FILEXIST
*
*     format statement
*
1000  FORMAT(//,2X,'Error: Input File Not Found.',//,2X,
     &'Missing Input File: ',A40,//)

*
*     Open the specified input file.  If the file does not exist, print
*     an error message and terminate execution.
*
      INQUIRE (FILE=FILENAME,EXIST=FILEXIST)

      IF (FILEXIST) THEN
         OPEN (LDNUM,FILE=FILENAME)
      ELSE
         WRITE(*,1000) FILENAME
         WRITE(LDECHO,1000) FILENAME
         STOP
      ENDIF

      RETURN
      END
