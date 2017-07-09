************************************************************************
*
*     Subroutine INPUT2              Called by: INPUT
*
*     input the parameters for each reach
*
*     Local Variables
*     ---------------
*     NSEG    number of segments in reach
*     END     distance at end of reach [L]
*
************************************************************************
      SUBROUTINE INPUT2(IMAX,DISP,AREA2,ALPHA,XSTART,DELTAX,NREACH,
     &                  LASTSEG, nseg, rchlen)
*
*     dimensional parameters and logical devices
*
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
*     local variables
*
      ! nseg, rchlen
      INTEGER*4 I,I2,I3,NSEG(MAXREACH)
      DOUBLE PRECISION RCHLEN(MAXREACH),END(0:MAXREACH)
      CHARACTER*500 BUFFER
*
*     subroutine arguments
*
      INTEGER*4 IMAX,NREACH,LASTSEG(0:MAXREACH)
      DOUBLE PRECISION XSTART,AREA2(*),DELTAX(*),ALPHA(*),DISP(*)
*
*     Input Format Statements
*
 1000 FORMAT(I5)
 1020 FORMAT(I5,4F13.0)
*
*     Output Format Statements
*
 2000 FORMAT(///,10X,'Number of Reaches  : ',I5,//)
 2020 FORMAT(10X,'Reach Information',/,10X,17('-'),//,10X,
     &'Reach   Start    End    Start           End',/,10X,
     &'No.     Seg.     Seg.   Distance [L]    Distance [L]'
     &,/,10X,52('-'))
 2040 FORMAT(10X,3(I5,3X),F9.3,7X,F9.3)
 2060 FORMAT(///,10X,'Number of Spatial Segments    : ',I4)
 2080 FORMAT(//,10X,'Reach Parameters',/,10X,16('-'),//,10X,
     &'Reach  No.     Reach         Segment      Dispersion',/,10X,
     &'No.    Segs.   Length [L]    Length [L]   Coeff.[L^2/s]',
     &/,10X,56('-'))
 2100 FORMAT (9X,I5,3X,I5,2(1X,F12.3),2X,1PE13.5)
 2120 FORMAT (//,10X,'Storage Parameters',/,10X,18('-'),/,17X,
     & 'Storage      Storage',/,10X,'Reach  Zone Area    Rate',/,10X,
     &'No.    [L^2]        [/sec]',/,10X,33('-'))
 2140 FORMAT (10X,I5,2(2X,1PE11.5))
*
*     Read the number of reaches
*
      CALL GETLINE(LDPARAM,BUFFER)
      READ (BUFFER,1000) NREACH
      WRITE (LDECHO,2000) NREACH
      IF (NREACH .GT. MAXREACH) CALL ERROR(1,NREACH,MAXREACH)
*
*     Read number of segments, reach lengths [L], dispersion coefficient
*     [L^2/sec], area of storage zone [L^2] and storage exchange
*     coefficient [/sec] for each reach.
*
*     Compute the index for the last segment of each reach, LASTSEG,
*     the segment length, DELTAX, and the distance at the end of each
*     reach, END.
*
*     Perform an error test - storage area must be nonzero.
*     Note:  to include a reach that does not have a storage
*            zone, Alpha should be set to zero.  The storage zone
*            area parameter must be nonzero, however, to avoid
*            division by zero.
*
*     Note that the DELTAX, DISP, AREA2, and ALPHA vectors are
*     filled in a subsequent loop.
*

      I = 1
      LASTSEG(0) = 0

      END(0) = XSTART
      WRITE(LDECHO,2020)

      DO 10 I2 = 1,NREACH
         CALL GETLINE(LDPARAM,BUFFER)
         !write(*,*) "inside input2()"
         READ(BUFFER,1020) NSEG(I2),RCHLEN(I2),DISP(I),AREA2(I),ALPHA(I)
         !write(*,*) "I2, NSEG(I2) -->", I2, NSEG(I2)
         !write(*,*) "inside input2()--2"

         LASTSEG(I2) = LASTSEG(I2-1) + NSEG(I2)
         IF (AREA2(I).LE.0.D0) THEN
            CALL ERROR2(1,AREA2(I),AREA2(I))
         ELSEIF (LASTSEG(I2).GT.MAXSEG) THEN
            CALL ERROR(2,LASTSEG(I2),MAXSEG)
         ENDIF
         DELTAX(I) = RCHLEN(I2) / DBLE(NSEG(I2))
         END(I2) = END(I2-1) + RCHLEN(I2)
         WRITE(LDECHO,2040) I2,I,LASTSEG(I2),END(I2-1),END(I2)
         I = LASTSEG(I2) + 1
         ! SM: I = 501
 10   CONTINUE


*
*     Compute the number of segments
*
      IMAX = LASTSEG(NREACH)
      WRITE (LDECHO, 2060) IMAX
*
*     Echo input values
*
      ! SM: Finding out what LASTSEG contains:
      ! print *, "Inside input2.f -- printing LASTSEG"
      ! print *, LASTSEG(0:MAXREACH)
      ! LASTSEG(0) = 0; LASTSEG(1) = 500

      I = 1
      WRITE(LDECHO, 2080)

      DO 20 I2 = 1, NREACH
         WRITE(LDECHO, 2100) I2,NSEG(I2),RCHLEN(I2),DELTAX(I),DISP(I)
         I = LASTSEG(I2) + 1
 20   CONTINUE
*
*      Fill vectors with reach values, echo storage parameters
*
      I = 1
      WRITE(LDECHO, 2120)

      DO 40 I2 = 1,NREACH
         WRITE(LDECHO,2140) I2,AREA2(I),ALPHA(I)
         DO 30 I3 = I+1,LASTSEG(I2)
            DELTAX(I3) = DELTAX(I)
            DISP(I3) = DISP(I)
            AREA2(I3) = AREA2(I)
            ALPHA(I3) = ALPHA(I)
 30      CONTINUE
         I = LASTSEG(I2) + 1
 40   CONTINUE

      RETURN
      END
