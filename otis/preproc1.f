************************************************************************
*     
*     Subroutine PREPROC1            Called by: PREPROC
*
*     compute the time invariant parameter groups used within the finite
*     difference subroutines (PREPROC3, CONSER/REACT).  As these groups 
*     are constant in time, they may be computed here on a one-time 
*     basis, rather than at each timestep.
*
************************************************************************
      SUBROUTINE PREPROC1(DISP,DFACE,DELTAX,HPLUSF,HPLUSB,HMULTF,HMULTB,
     &                    HDIV,HDIVF,HADV,IMAX)
      INCLUDE 'fmodules.inc'
*
*     subroutine arguments
*
      INTEGER*4 IMAX
      DOUBLE PRECISION DISP(*),DFACE(*),DELTAX(*),HPLUSF(*),HPLUSB(*),
     &                 HMULTF(*),HMULTB(*),HDIV(*),HDIVF(*),HADV(*)
*
*     local variables
*
      INTEGER*4 I
*
*     compute the interfacial dispersion coefficient
*     
      DO 10 I = 1,IMAX-1
         DFACE(I) =  DELTAX(I)/(DELTAX(I+1)+DELTAX(I)) * DISP(I+1) 
     &             + DELTAX(I+1)/(DELTAX(I+1)+DELTAX(I)) * DISP(I) 
 10   CONTINUE

      DFACE(IMAX) = DISP(IMAX)
*
*     compute the H-groups
*
      HPLUSF(1) = DELTAX(1) + DELTAX(2)
      HMULTF(1) = DELTAX(1) * HPLUSF(1)
      HDIV(1) = DELTAX(1) / HPLUSF(1)
      HDIVF(1) = DELTAX(2) / HPLUSF(1)

      DO 20 I = 2, IMAX-1
         HPLUSF(I) = DELTAX(I) + DELTAX(I+1)
         HPLUSB(I) = DELTAX(I) + DELTAX(I-1) 
         HMULTF(I) = DELTAX(I) * HPLUSF(I)
         HMULTB(I) = DELTAX(I) * HPLUSB(I)
         HDIV(I) = DELTAX(I) / HPLUSF(I)
         HDIVF(I) = DELTAX(I+1) / HPLUSF(I)         
         HADV(I) = 0.5D0 * (HDIVF(I)-DELTAX(I-1)/HPLUSB(I)) / DELTAX(I)
 20   CONTINUE

      I = IMAX
      HPLUSB(I) = DELTAX(I) + DELTAX(I-1) 
      HMULTB(I) = DELTAX(I) * HPLUSB(I)

      RETURN
      END



