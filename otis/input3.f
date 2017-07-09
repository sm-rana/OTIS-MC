************************************************************************
*
*     Subroutine INPUT3               Called by: INPUT
*
*     input the number of solutes and the flags indicating the type of
*     chemistry being considered.
*
************************************************************************
      SUBROUTINE INPUT3(NSOLUTE,IDECAY,ISORB)
*
*     dimensional parameters and logical devices
*
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
*     subroutine arguments
*
      INTEGER*4 NSOLUTE,IDECAY,ISORB
*
*     local variables
*
      CHARACTER*500 BUFFER
*
*     input Format Statements
*
 1000 FORMAT(3I5)
*
*     Output Format Statements
*
 2000 FORMAT(///,10X,'Number of Solutes: ',I3)
 2020 FORMAT(10X,'Decay Option     : ',I3)
 2040 FORMAT(10X,'Sorption Option  : ',I3)
*
*     Read the number of solutes, decay flag and sorption flag
*
      CALL GETLINE(LDPARAM,BUFFER)
      READ (BUFFER,1000) NSOLUTE,IDECAY,ISORB
      WRITE (LDECHO,2000) NSOLUTE
      WRITE (LDECHO,2020) IDECAY
      WRITE (LDECHO,2040) ISORB
*
*     error check
*
      IF (NSOLUTE .GT. MAXSOLUTE) CALL ERROR(4,NSOLUTE,MAXSOLUTE)
      IF (IDECAY .LT. 0 .OR. IDECAY .GT. 1) CALL ERROR3(4,IDECAY)
      IF (ISORB .LT. 0 .OR. ISORB .GT. 1) CALL ERROR3(5,ISORB)
 
      RETURN
      END
