************************************************************************
*
*     Subroutine HEADER              Called by: MAININIT
*
*     Open file for input echoing (echo.out) & print model info
*
************************************************************************
      SUBROUTINE HEADER(MODDESC)
*
*     dimensional parameters and logical devices
*
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
*     model descriptive string
*
      CHARACTER*(*) MODDESC
*
*     local variables
*
      CHARACTER*30 VERSION

      DATA VERSION /'Version: MOD40 (May 1997)'/

 1000 FORMAT(//,10X,A51,/,22X,A30,/,22X,26('-'))
*
*     Open output file for run information and input echoing
*
      OPEN (UNIT=LDECHO,FILE='echo.out',STATUS='UNKNOWN')
*
*     Write Info to file
*
      WRITE(LDECHO,1000) MODDESC,VERSION

      RETURN
      END
