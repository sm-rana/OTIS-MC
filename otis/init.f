************************************************************************
*
*     Subroutine INIT              Called by: MAININIT
*
*     initialize chemical parameters
*
************************************************************************
      SUBROUTINE INIT(LAMBDA,LAMBDA2,LAMHAT,LAMHAT2,RHOLAM,KD,CSBACK)
*
*     dimensional parameters
*
      INCLUDE 'fmodules.inc'
*
*     subroutine arguments
* 
      DOUBLE PRECISION LAMBDA(MAXSEG,*),LAMBDA2(MAXSEG,*),
     &                 LAMHAT(MAXSEG,*),LAMHAT2(MAXSEG,*),
     $                 RHOLAM(MAXSEG,*),KD(MAXSEG,*),CSBACK(MAXSEG,*)
*
*     local variables
*
      INTEGER*4 I,J
*
*     initialize
*
      DO 20 I=1,MAXSEG
         DO 10 J=1,MAXSOLUTE
            LAMBDA(I,J) = 0.D0
            LAMBDA2(I,J) = 0.D0
            LAMHAT(I,J) = 0.D0
            LAMHAT2(I,J) = 0.D0
            RHOLAM(I,J) = 0.D0
            KD(I,J) = 0.D0
            CSBACK(I,J) = 0.D0
 10      CONTINUE
 20   CONTINUE

      RETURN
      END
