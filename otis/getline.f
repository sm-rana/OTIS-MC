************************************************************************
*
*     Subroutine GETLINE               Called by: INPUT1-7 etc.
*
*     read from the input file designated by LDUNIT, discarding comment
*     lines and returning the first line of valid input data.
*
************************************************************************
      SUBROUTINE GETLINE(LDUNIT,BUFFER)
*
*     subroutine arguments
*
      INTEGER*4 LDUNIT
      CHARACTER*500 BUFFER
*
*     format statement
*
 1000 FORMAT(A500)
*
*     read line from input file.  If a # appears in the first column the
*     line is a comment and the next line is read.  Otherwise return to
*     calling routine with BUFFER containing input data
*
 10   READ (LDUNIT,1000) BUFFER
      IF (INDEX(BUFFER,'#') .EQ. 1) GOTO 10

      RETURN
      END
