************************************************************************
*
*     Subroutine QSTEADY              Called by: INPUTQ
*
*     A subroutine to input the following parameters:
*
*     QSTART, QLATIN, QLATOUT, AREA, CLATIN
*
************************************************************************
      SUBROUTINE QSTEADY(QSTART,QLATIN,QLATOUT,CLATIN,AREA,NREACH,
     &                   LASTSEG,NSOLUTE)
*
*     dimensional parameters and logical devices
*
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
*     subroutine arguments
*
      INTEGER*4 NREACH,NSOLUTE,LASTSEG(0:MAXREACH)
      DOUBLE PRECISION QSTART,AREA(*),QLATIN(*),QLATOUT(*),
     &                 CLATIN(MAXSEG,*)
*
*     local variables
*
      INTEGER*4 I,I2,I3,J
      CHARACTER*500 BUFFER
*
*     Input Format Statements
*
 1000 FORMAT(70D13.5)
*
*     Output Format Statements
*
 2000 FORMAT(10X,'Flowrate at Upstream Boundary    [L^3/sec]: ',1PE13.5)
 2020 FORMAT(//,20X,'Lateral         Lateral         X-sect.',
     & /,10X,
     &'Reach     Inflow          Outflow         Area',
     & /,10X,
     &'No.       [L^3/s-L]       [L^3/s-L]       [L^2]',
     & /,10X,55('-'))
 2040 FORMAT(10X,I5,1P,3(5X,E11.5))
 2060 FORMAT(//,10X,'Initial Lateral Inflow Concentrations, Solute #',
     &       I3,/,10X,51('-'))
 2080 FORMAT(21X,'Inflow',/,10X,'Reach      Concent.',/,10X,
     &'No.        [mass/L^3]',/,10X,23('-'))
 2100 FORMAT(10X,I4,7X,1PE11.5)
*
*     Read the flow rate at the upstream boundary [L^3/sec].
*
      CALL GETLINE(LDFLOW,BUFFER)
      READ (BUFFER,1000) QSTART
*
*     Read lateral inflow, lateral outflow, area, and lateral
*     inflow concentrations for each reach.
*
*     Fill vectors with reach values, and echo flow values
*
      WRITE(LDECHO,2000) QSTART
      WRITE(LDECHO,2020)
      I = 1

      DO 30 I2 = 1,NREACH
         CALL GETLINE(LDFLOW,BUFFER)
         READ(BUFFER,1000) QLATIN(I),QLATOUT(I),AREA(I),
     &                     (CLATIN(I,J),J=1,NSOLUTE)
        !SM:
*         print *, "In qsteady: Area(I) = " , AREA(I), "I" ,I
*         print *,"In qsateady: Nreach = ", NREACH
         DO 20 I3 = I+1,LASTSEG(I2)
            AREA(I3) = AREA(I)
            QLATIN(I3) = QLATIN(I)
            QLATOUT(I3) = QLATOUT(I)
            DO 10 J = 1,NSOLUTE
               CLATIN(I3,J) = CLATIN(I,J)
 10         CONTINUE
 20      CONTINUE
         WRITE(LDECHO,2040) I2,QLATIN(I),QLATOUT(I),AREA(I)
         I = LASTSEG(I2) + 1
 30   CONTINUE
*
*     echo lateral inflow concentrations
*
      DO 50 J = 1, NSOLUTE
         WRITE(LDECHO,2060) J
         WRITE(LDECHO,2080)
         I = 1
         DO 40 I2 = 1, NREACH
            WRITE(LDECHO,2100) I2,CLATIN(I,J)
            I = LASTSEG(I2) + 1
 40      CONTINUE
 50   CONTINUE

      RETURN
      END
