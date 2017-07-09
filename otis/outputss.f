************************************************************************
*
*     Subroutine OUTPUTSS              Called by: MAININIT, MAINRUN
*
*     output steady-state solute concentrations
*
************************************************************************
      SUBROUTINE OUTPUTSS(CONC,CONC2,PRTOPT,ISORB,SORB,NSOLUTE,IMAX,
     &                    DIST)
*
*     dimensional parameters and logical devices
*
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
*     subroutine arguments
*
      INTEGER*4 NSOLUTE,PRTOPT,ISORB,IMAX
      DOUBLE PRECISION DIST(*)
      DOUBLE PRECISION CONC(MAXSEG,*),CONC2(MAXSEG,*),SORB(MAXSEG,*)
*
*     local variables
*
      INTEGER*4 I,J
 100  FORMAT(1P70E14.6)
*
*     for each solute
*
      DO 40 J=1,NSOLUTE
*
*        write out solute concentrations
*
         IF (PRTOPT .EQ. 1) THEN
            DO 10 I = 1,IMAX
               WRITE (LDFILES(J),100) DIST(I),CONC(I,J)
 10         CONTINUE
         ELSE
            DO 20 I = 1,IMAX
               WRITE (LDFILES(J),100) DIST(I),CONC(I,J),CONC2(I,J)
 20         CONTINUE
         ENDIF
*
*        write out sorbed concentrations
*
         IF (ISORB .EQ. 1) THEN
            DO 30 I = 1,IMAX
               WRITE (LDSORB(J),100) DIST(I),SORB(I,J)
 30         CONTINUE
         ENDIF

 40   CONTINUE

      RETURN
      END
