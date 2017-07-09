************************************************************************
*
*     Subroutine INPUT5               Called by: INPUT
*
*     read the sorption parameters
*
************************************************************************
      SUBROUTINE INPUT5(NREACH,LASTSEG,NSOLUTE,RHOLAM,LAMHAT,LAMHAT2,KD,
     &                  CSBACK)
*
*     dimensional parameters and logical devices
*
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
*     subroutine arguments
*
      INTEGER*4 NREACH,NSOLUTE,LASTSEG(0:MAXREACH)
      DOUBLE PRECISION RHOLAM(MAXSEG,*),LAMHAT(MAXSEG,*),
     $                 LAMHAT2(MAXSEG,*),KD(MAXSEG,*),CSBACK(MAXSEG,*)
*
*     local variables
*
      INTEGER*4 I,I2,I3,J
      DOUBLE PRECISION RHO(MAXSEG,MAXSOLUTE)
      CHARACTER*500 BUFFER
*
*     input Format Statements
*
 1000 FORMAT(5D13.5)
*
*     Output Format Statements
*
 2000 FORMAT(//,10X,'Sorption Parameters, Solute #',I3,/,10X,70('-'),/,
     &19X,'First-order Rate Coeff.              ',
     &'Partition   Background',/,10X,'Reach      Channel     ',
     &'Storage      Rho         Coeff.    S.Z. Conc.',
     &/,10X,'No.       [/second]   [/second]    ',
     &   '[M/L^3]     [L^3/M]     [M/L^3]',/,10X,70('-'))
 2020 FORMAT(10X,I5,2X,5(3X,1PE9.3))
*
*     Read the sorption parameters for each reach and solute.  Note that
*     the parameter RHO is always multiplied by LAMHAT -- therefore do
*     the multiplication here and store in RHOLAM.  Fill vectors with
*     reach values.
*
      DO 30 J=1,NSOLUTE
         I = 1
         DO 20 I2 = 1,NREACH
            CALL GETLINE(LDPARAM,BUFFER)
            READ(BUFFER,1000) LAMHAT(I,J),LAMHAT2(I,J),RHO(I,J),KD(I,J),
     &                        CSBACK(I,J)
            RHOLAM(I,J) = LAMHAT(I,J) * RHO(I,J)
            DO 10 I3 = I+1,LASTSEG(I2)
               LAMHAT(I3,J) = LAMHAT(I,J)
               LAMHAT2(I3,J) = LAMHAT2(I,J)
               RHOLAM(I3,J) = RHOLAM(I,J)
               KD(I3,J) = KD(I,J)
               CSBACK(I3,J) = CSBACK(I,J)
 10         CONTINUE
            I = LASTSEG(I2) + 1
 20      CONTINUE
 30   CONTINUE
*
*     echo sorption parameters
*
      DO 50 J = 1,NSOLUTE
         WRITE(LDECHO,2000) J
         I = 1
         DO 40 I2 = 1, NREACH
            WRITE(LDECHO,2020) I2,LAMHAT(I,J),LAMHAT2(I,J),RHO(I,J),
     &                         KD(I,J),CSBACK(I,J)
            I = LASTSEG(I2) + 1
 40      CONTINUE
 50   CONTINUE
 
      RETURN
      END
