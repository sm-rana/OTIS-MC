************************************************************************
*
*     Subroutine change_area
*     It replicates the behavior of QSTEADY to change X-Sectional Area for
*     Monte Carlo Simulation
*
*     Called by: INPUTQ
*
*
************************************************************************
      SUBROUTINE change_area(NEW_A, AREA, NREACH, LASTSEG)
*     SUBROUTINE QSTEADY2(QSTART,QLATIN,QLATOUT,CLATIN,AREA,NREACH,
*     &                   LASTSEG,NSOLUTE)
*
*     dimensional parameters and logical devices
*
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
*     subroutine arguments
*
      INTEGER*4 NREACH,LASTSEG(0:MAXREACH)
      DOUBLE PRECISION AREA(*)

      double precision NEW_A
*
*     local variables
*
      INTEGER*4 I,I2,I3

      !write(*,*) "change_area()--> NEW_A", NEW_A
      !print *,"In change_area: Nreach = ", NREACH

      I = 1
      DO 870 I2 = 1,NREACH
         AREA(I) = NEW_A
         !print *, "In change_area: Area(I) = " , AREA(I), "I" ,I
         DO 871 I3 = I+1,LASTSEG(I2)
            AREA(I3) = AREA(I)
 871     CONTINUE
         I = LASTSEG(I2) + 1
 870  CONTINUE

      !print *, "End of change_area()"

      RETURN
      END
