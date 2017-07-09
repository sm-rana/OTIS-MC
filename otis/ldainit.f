************************************************************************
*
*     Subroutine LDAINIT             Called by: MAININIT
*
*     initialize Logical Device Assignments (LDAs)
*
************************************************************************
      SUBROUTINE LDAINIT
*
*     dimensional parameters and logical devices
*
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
      INTEGER*4 I,J
*
*     Assign Input Unit Numbers
*
      LDCTRL = 10
      LDPARAM = 11
      LDFLOW = 12
*
*     Assign Ouput Unit Numbers
*
      LDECHO = 30
      J = 40

      DO 10 I = 1, MAXSOLUTE
         LDFILES(I) = J
         LDSORB(I) = J + 20
         J = J + 1
 10   CONTINUE

      RETURN
      END
