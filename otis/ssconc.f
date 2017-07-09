************************************************************************
*
*     Subroutine SSCONC               Called by: STEADY
*
*     compute steady-state storage zone and streambed sediment
*     concentrations
*
************************************************************************
      SUBROUTINE SSCONC(IMAX,CONC,CONC2,ALPHA,AREA,AREA2,LAMBDA2,
     &                  NSOLUTE,LAMHAT2,CSBACK,KD,SORB)
      INCLUDE 'fmodules.inc'
*
*     subroutine arguments
*
      INTEGER IMAX,NSOLUTE
      DOUBLE PRECISION ALPHA(*),AREA(*),AREA2(*),CONC(MAXSEG,*),
     &                 CONC2(MAXSEG,*),LAMBDA2(MAXSEG,*),
     &                 LAMHAT2(MAXSEG,*),CSBACK(MAXSEG,*),KD(MAXSEG,*),
     &                 SORB(MAXSEG,*)
*
*     local variables
*
      INTEGER*4 I,J
*
*     compute the steady-state storage zone and sorbate concentrations.
*     If the exchange coefficient is zero, the storage zone
*     concentration is zero.
*
      DO 20 J = 1,NSOLUTE
         DO 10 I = 1,IMAX
            IF (ALPHA(I) .NE. 0.D0 .OR. LAMHAT2(I,J) .NE. 0.D0) THEN
               CONC2(I,J) = (ALPHA(I)*AREA(I)*CONC(I,J)
     &                        + AREA2(I)*LAMHAT2(I,J)*CSBACK(I,J)) /
     &                       (ALPHA(I)*AREA(I) +
     &                        AREA2(I)*(LAMBDA2(I,J)+LAMHAT2(I,J)))
            ELSE
               CONC2(I,J) = 0.D0
            ENDIF
            SORB(I,J) = KD(I,J)*CONC(I,J)
 10      CONTINUE
 20   CONTINUE

      RETURN
      END
