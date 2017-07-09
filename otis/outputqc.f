************************************************************************
*
*     Subroutine OUTPUTQC              Called by: MAINRUN 
*
*     output solute flux, Q*C, such that mass conservation may be
*     checked
*
************************************************************************
      SUBROUTINE OUTPUTQC(NPRINT,PINDEX,CONC,TIME,WT,JSOLUTE,QCNC,Q)
*
*     dimensional parameters and logical devices
*
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
*     subroutine arguments
*
      INTEGER*4 NPRINT,JSOLUTE,PINDEX(*)
      DOUBLE PRECISION TIME,WT(*),QCNC(*),Q(*),CONC(MAXSEG,*)
*
*     local variables
*
      INTEGER*4 I,J
 100  FORMAT(70E14.6)
*
*     write out solute concentrations
*
         DO 10 J = 1,NPRINT
            I = PINDEX(J)
            QCNC(J) = CONC(I,JSOLUTE) * Q(I)
     &             + WT(J) * (CONC(I+1,JSOLUTE)*Q(I+1)
     &                         - CONC(I,JSOLUTE)*Q(I))
 10      CONTINUE
         WRITE (LDFILES(JSOLUTE),100) TIME,(QCNC(J),J=1,NPRINT)

      RETURN
      END
