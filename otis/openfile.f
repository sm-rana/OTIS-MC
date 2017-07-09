************************************************************************
*
*     Subroutine OPENFILE            Called by: MAININIT
*
*     Open input files
*
************************************************************************
      SUBROUTINE OPENFILE
*
*     dimensional paramters and logical devices
*
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
*     local variables
*
      CHARACTER*40 FILE
      CHARACTER*500 BUFFER
*
*     format statements
*
 1000 FORMAT(A40)
 2000 FORMAT(///,10X,'Input Files',/,10X,11('-'))
 2010 FORMAT(10X,'PARAMETER FILE: ',A40)
 2020 FORMAT(10X,'FLOW FILE     : ',A40)
*
*     Open the control file
*
      CALL OPENIN(LDCTRL,'control.inp')
*
*     write header, determine the parameter file and open
*     sm : reading from the control.inp file -- Flow file and Params File
      WRITE (LDECHO,2000)
      CALL GETLINE(LDCTRL,BUFFER)
      READ (BUFFER,1000) FILE
      WRITE (LDECHO,2010) FILE
      CALL OPENIN(LDPARAM,FILE)
*
*     Determine the flow file and open
*
      CALL GETLINE(LDCTRL,BUFFER)
      READ (BUFFER,1000) FILE
      WRITE (LDECHO,2020) FILE
      CALL OPENIN(LDFLOW,FILE)

      RETURN
      END
