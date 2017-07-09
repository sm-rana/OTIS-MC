************************************************************************
*
*     Subroutine OUTPUT              Called by: MAINRUN, OUTINIT
*
*     output solute concentrations
*
************************************************************************
      SUBROUTINE OUTPUT(NPRINT,PINDEX,CONC,CONC2,TIME,WT,JSOLUTE,PRTOPT,
     &                  CNC,CNC2,ISORB,SORB)
*
*     dimensional parameters and logical devices
*
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
*     subroutine arguments
*
      INTEGER*4 NPRINT,JSOLUTE,PRTOPT,ISORB,PINDEX(*)
      DOUBLE PRECISION TIME,WT(*),CNC(*),CNC2(*),CONC(MAXSEG,*),
     &                 CONC2(MAXSEG,*),SORB(MAXSEG,*)
*
*     local variables
*
      INTEGER*4 I,J
 100  FORMAT(1P70E14.6)
*
*     write out solute concentrations
*
      IF (PRTOPT .EQ. 1) THEN
         DO 10 J = 1,NPRINT
            I = PINDEX(J)
            CNC(J) = CONC(I,JSOLUTE)
     &             + WT(J) * (CONC(I+1,JSOLUTE) - CONC(I,JSOLUTE))
 10      CONTINUE
         WRITE (LDFILES(JSOLUTE),100) TIME,(CNC(J),J=1,NPRINT)
      ELSE
         DO 20 J = 1,NPRINT
            I = PINDEX(J)
            CNC(J) = CONC(I,JSOLUTE)
     &             + WT(J) * (CONC(I+1,JSOLUTE) - CONC(I,JSOLUTE))
            CNC2(J) = CONC2(I,JSOLUTE) + WT(J)
     &                     * (CONC2(I+1,JSOLUTE)-CONC2(I,JSOLUTE))
 20      CONTINUE

         WRITE (LDFILES(JSOLUTE),100) TIME,(CNC(J),J=1,NPRINT),
     &                                     (CNC2(J),J=1,NPRINT)
      ENDIF
*
*     write out sorbed concentrations
*
      IF (ISORB .EQ. 1) THEN
         DO 30 J = 1,NPRINT
            I = PINDEX(J)
            CNC(J) = SORB(I,JSOLUTE)
     &             + WT(J) * (SORB(I+1,JSOLUTE) - SORB(I,JSOLUTE))
 30      CONTINUE
         WRITE (LDSORB(JSOLUTE),100) TIME,(CNC(J),J=1,NPRINT)
      ENDIF

      RETURN
      END
