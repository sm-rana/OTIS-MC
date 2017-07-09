************************************************************************
*
*     Subroutine INPUT              Called by: MAININIT
*
*     input physical and chemical parameters
* SM: assigns values to the variables and do not do any sort of processing or calculations
************************************************************************
      SUBROUTINE INPUT(TSTEP,TSTART,TFINAL,IPRINT,XSTART,DSBOUND,PRTOPT,
     &                 IMAX,DISP,AREA2,ALPHA,DELTAX,NSOLUTE,LAMBDA,
     &                 LAMBDA2,NPRINT,PRTLOC,NBOUND,USTIME,USCONC,
     &                 USCONCN,QSTEP,QSTART,QLATIN,QLATOUT,CLATIN,AREA,
     &                 FLOWLOC,QVALUE,AVALUE,NFLOW,QINVAL,CLVAL,PINDEX,
     &                 WT,IBOUND,USBC,CSBACK,RHOLAM,LAMHAT,LAMHAT2,KD,
     &                 ISORB,IDECAY,DIST)
      INCLUDE 'fmodules.inc'
*
*     subroutine arguments
*
      INTEGER*4 IMAX,NPRINT,IPRINT,NBOUND,NSOLUTE,PRTOPT,NFLOW,IBOUND,
     &          ISORB,IDECAY
      INTEGER*4 PINDEX(*)
      DOUBLE PRECISION DSBOUND,TSTART,TFINAL,TSTEP,QSTART,XSTART,QSTEP
      DOUBLE PRECISION PRTLOC(*),AREA2(*),DELTAX(*),USTIME(*),AREA(*),
     &                 ALPHA(*),DISP(*),QLATIN(*),QLATOUT(*),FLOWLOC(*),
     &                 QVALUE(*),AVALUE(*),QINVAL(*),USCONCN(*),WT(*),
     &                 DIST(*)
      DOUBLE PRECISION USCONC(MAXBOUND,*),CLATIN(MAXSEG,*),
     &                 LAMBDA(MAXSEG,*),LAMBDA2(MAXSEG,*),
     &                 CLVAL(MAXFLOWLOC+1,*),USBC(MAXBOUND,*),
     &                 RHOLAM(MAXSEG,*),LAMHAT(MAXSEG,*),
     $                 LAMHAT2(MAXSEG,*),KD(MAXSEG,*),CSBACK(MAXSEG,*)
*
*     local variables
*
      INTEGER*4 NREACH,LASTSEG(0:MAXREACH)
*
*     read the print option, the time parameters, etc.
*
      CALL INPUT1(TSTEP,TSTART,TFINAL,IPRINT,XSTART,DSBOUND,PRTOPT)
*
*     read the time-invariant parameters for each reach
*
      CALL INPUT2(IMAX,DISP,AREA2,ALPHA,XSTART,DELTAX,NREACH,LASTSEG)
*
*     read the number of solutes and the chemistry flags
*
      CALL INPUT3(NSOLUTE,IDECAY,ISORB)
*
*     read the first-order decay rates
*
      IF (IDECAY .EQ. 1)
     &   CALL INPUT4(NREACH,LASTSEG,NSOLUTE,LAMBDA,LAMBDA2)
*
*     read the sorption parameters
*
      IF (ISORB .EQ. 1) CALL INPUT5(NREACH,LASTSEG,NSOLUTE,RHOLAM,
     &                              LAMHAT,LAMHAT2,KD,CSBACK)
*
*     read the print parameters
*
      CALL INPUT6(NPRINT,PRTLOC,IMAX,DELTAX,XSTART,PINDEX,WT,DIST)
*
*     read the flow variables at each reach or flow location
*
      CALL INPUTQ(QSTEP,QSTART,QLATIN,QLATOUT,CLATIN,AREA,NREACH,
     &            LASTSEG,NSOLUTE,FLOWLOC,QVALUE,AVALUE,NFLOW,QINVAL,
     &            CLVAL)
*
*     read the upstream boundary conditions
*
      CALL INPUT7(NBOUND,USTIME,USCONC,USCONCN,NSOLUTE,IBOUND,USBC,
     &            QSTART,TFINAL)

      RETURN
      END
