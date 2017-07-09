************************************************************************
*
*     Subroutine OUTINIT             Called by: MAININIT, MAINRUN
*
*     write the initial conditions to the echo file and to the output
*     files.
*
************************************************************************
      SUBROUTINE OUTINIT(NPRINT,PINDEX,PRTLOC,Q,CONC,CONC2,TSTART,WT,
     &                   NSOLUTE,PRTOPT,CHEM,ISORB,SORB)
*
*     dimensional parameters and logical devices
*
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
*     subroutine arguments
*
      INTEGER*4 NPRINT,NSOLUTE,PRTOPT,ISORB,PINDEX(*)
      DOUBLE PRECISION TSTART,PRTLOC(*),Q(*),WT(*),CONC(MAXSEG,*),
     &                 CONC2(MAXSEG,*),SORB(MAXSEG,*)
      CHARACTER*(*) CHEM
*
*     local variables
*
      INTEGER*4 I,J,JSOLUTE
      DOUBLE PRECISION CNC(MAXPRINT),CNC2(MAXPRINT),QQ
*
*     output format statements
*
 100  FORMAT(//,10X,'Initial Conditions at the Print Locations',/,10X,
     &       41('-'),//,10X,'Location  Prev.    Inter.         Flow',
     &/,10X,'[L]       Seg.     Weight        [L^3/sec]',/,10X,46('-'))
 200  FORMAT(10X,F9.2,1X,I5,2(3X,1PE12.6))
 300  FORMAT(//,10X,'Initial Conditions, Solute #',I3,/,10X,32('-'),/,
     &30X,'Concentration',/,10X,'Location        Channel       Storage',
     &/,10X,'[L]            [mass/L^3]    [mass/L^3]',/,10X,41('-'))
 400  FORMAT(10X,F9.2,7X,1P,E9.3,5X,E9.3)
*
*     write out interpolated initial flows
*
      WRITE (LDECHO,100)
      DO 10 J = 1,NPRINT
         I = PINDEX(J)
         QQ = Q(I) + WT(J) * (Q(I+1) - Q(I))
         WRITE (LDECHO,200) PRTLOC(J),I,WT(J),QQ
 10   CONTINUE
*
*     If this is not a steady-state run, write out interpolated initial
*     concentrations.
*
      IF (CHEM .NE. 'Steady-State') THEN
         DO 30 JSOLUTE = 1,NSOLUTE
            WRITE(LDECHO,300) JSOLUTE
            CALL OUTPUT(NPRINT,PINDEX,CONC,CONC2,TSTART,WT,JSOLUTE,
     &                  PRTOPT,CNC,CNC2,ISORB,SORB)
            DO 20 J = 1,NPRINT
               IF (PRTOPT .EQ. 1) THEN
                  WRITE (LDECHO,400) PRTLOC(J),CNC(J)
               ELSE
                  WRITE (LDECHO,400) PRTLOC(J),CNC(J),CNC2(J)
               ENDIF
 20         CONTINUE
 30      CONTINUE
      ENDIF

      RETURN
      END
