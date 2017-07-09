************************************************************************
*
*     Subroutine INPUT2_MOD              Called by: MAININIT_MOD
*
*     input the parameters for each reach
*
*     Local Variables
*     ---------------
*     NSEG    number of segments in reach
*     END     distance at end of reach [L]
*
************************************************************************
      SUBROUTINE input2_mod(IMAX,DISP,AREA2,ALPHA,XSTART,DELTAX,NREACH,
     &                  LASTSEG, new_disp, new_area2, new_alpha,
     &                  nseg_tmp, rchlen_tmp)
*
*     dimensional parameters and logical devices
*
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'

      integer*4 nseg_tmp
      double precision new_disp, new_area2, new_alpha, rchlen_tmp
*
*     local variables
*
      INTEGER*4 I,I2,I3,NSEG(MAXREACH)


*
*     subroutine arguments
*
      INTEGER*4 IMAX,NREACH,LASTSEG(0:MAXREACH)
      DOUBLE PRECISION XSTART,AREA2(*),DELTAX(*),ALPHA(*),DISP(*)
      !SM: needs to be passed on from initial run of MAININIT
      DOUBLE PRECISION RCHLEN(MAXREACH),END(0:MAXREACH)


      ! SM: make sure NRECH is passed on to this function with value
      ! READ (BUFFER,1000) NREACH

      !write(*,*) "input2_mod(): nseg_tmp -->", nseg_tmp
      !write(*,*) "input2_mod(): NREACH -->", NREACH
      I = 1
      LASTSEG(0) = 0
      END(0) = XSTART

      DO 10 I2 = 1,NREACH
         !SM:
         DISP(I) = new_disp
         AREA2(I) = new_area2
         ALPHA(I) = new_alpha
         NSEG(I2) = nseg_tmp
         RCHLEN(I2) = rchlen_tmp
         !READ(BUFFER,1020) NSEG(I2),RCHLEN(I2),DISP(I),AREA2(I),ALPHA(I)
         !write(*,*) "input2_mod()-- DO LOOP:"
         !write(*,*) "I2, NSEG(I2) -->", I2, NSEG(I2)
         LASTSEG(I2) = LASTSEG(I2-1) + NSEG(I2)
         IF (AREA2(I).LE.0.D0) THEN
            CALL ERROR2(1,AREA2(I),AREA2(I))
         ELSEIF (LASTSEG(I2).GT.MAXSEG) THEN
            CALL ERROR(2,LASTSEG(I2),MAXSEG)
         ENDIF
         DELTAX(I) = RCHLEN(I2) / DBLE(NSEG(I2))
         END(I2) = END(I2-1) + RCHLEN(I2)
         !WRITE(LDECHO,2040) I2,I,LASTSEG(I2),END(I2-1),END(I2)
         I = LASTSEG(I2) + 1
 10   CONTINUE
*
*     Compute the number of segments
*
      IMAX = LASTSEG(NREACH)
      !WRITE (LDECHO, 2060) IMAX
*
*     Echo input values
*
      ! SM: Finding out what LASTSEG contains:
      ! print *, "Inside input2.f -- printing LASTSEG"
      ! print *, LASTSEG(0:MAXREACH)
      ! LASTSEG(0) = 0; LASTSEG(1) = 500

      I = 1
      !WRITE(LDECHO, 2080)

      DO 20 I2 = 1, NREACH
         !WRITE(LDECHO, 2100) I2,NSEG(I2),RCHLEN(I2),DELTAX(I),DISP(I)
         I = LASTSEG(I2) + 1
 20   CONTINUE
*
*      Fill vectors with reach values, echo storage parameters
*
      I = 1
      !WRITE(LDECHO, 2120)

      DO 40 I2 = 1,NREACH
         !WRITE(LDECHO,2140) I2,AREA2(I),ALPHA(I)
         DO 30 I3 = I+1,LASTSEG(I2)
            DELTAX(I3) = DELTAX(I)
            DISP(I3) = DISP(I)
            AREA2(I3) = AREA2(I)
            ALPHA(I3) = ALPHA(I)
 30      CONTINUE
         I = LASTSEG(I2) + 1
 40   CONTINUE
      !write(*,*) "Successful END of input2_mod()"
      RETURN
      END
