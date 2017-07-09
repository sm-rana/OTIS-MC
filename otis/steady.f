************************************************************************
*
*     Subroutine STEADY           Called by: MAININIT, MAINRUN, OTISRUN
*
*     compute the initial conditions (for dynamic simulations) or the
*     steady-state concentrations (for steady-state simulations).
*
************************************************************************
      SUBROUTINE STEADY(IMAX,AREA2,DELTAX,Q,USCONC,AREA,ALPHA,DSBOUND,
     &                  QLATIN,CLATIN,CONC,CONC2,NSOLUTE,LAMBDA,LAMBDA2,
     &                  DFACE,HPLUSF,HPLUSB,HMULTF,HMULTB,HDIV,HDIVF,
     &                  HADV,LAMHAT2,CSBACK,KD,SORB)
      INCLUDE 'fmodules.inc'
*
*     subroutine arguments
*
      INTEGER*4 IMAX,NSOLUTE
      DOUBLE PRECISION DSBOUND,AREA2(*),DELTAX(*),Q(*),AREA(*),ALPHA(*),
     &                 QLATIN(*),DFACE(*),HPLUSF(*),HPLUSB(*),HMULTF(*),
     &                 HMULTB(*),HDIV(*),HDIVF(*),HADV(*),
     &                 USCONC(MAXBOUND,*),CLATIN(MAXSEG,*),
     &                 CONC(MAXSEG,*),CONC2(MAXSEG,*),LAMBDA(MAXSEG,*),
     &                 LAMBDA2(MAXSEG,*),LAMHAT2(MAXSEG,*),
     &                 CSBACK(MAXSEG,*),KD(MAXSEG,*),SORB(MAXSEG,*)
*
*     compute main channel concentrations using finite differences
*
      CALL SSDIFF(IMAX,DELTAX,Q,USCONC,ALPHA,AREA,AREA2,DFACE,DSBOUND,
     &            QLATIN,CLATIN,CONC,HPLUSF,HPLUSB,HMULTF,HMULTB,HDIV,
     &            HDIVF,HADV,NSOLUTE,LAMBDA,LAMBDA2,LAMHAT2,CSBACK)
*
*     compute steady-state storage zone and streambed sediment
*     concentrations
*
      CALL SSCONC(IMAX,CONC,CONC2,ALPHA,AREA,AREA2,LAMBDA2,NSOLUTE,
     &            LAMHAT2,CSBACK,KD,SORB)
      RETURN
      END





