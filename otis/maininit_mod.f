************************************************************************
*
*     Subroutine MAININIT            Called by: MAINOLD program
*
*     Obtain input parameters, initialize variables, and preprocess.
*     If the steady-state option has been invoked, output results
*
************************************************************************
      !SM: change_area() needs nreach and lastseg
      !SM: input2_mod() needs nseg_tmp, rchlen_tmp

      SUBROUTINE maininit_mod(IMAX,NPRINT,PINDEX,WT,IPRINT,DELTAX,Q,
     &   USTIME,
     &                    USCONC,USCONCN,AREA,DSBOUND,QLATIN,CLATIN,
     &                    TSTART,TFINAL,TSTEP,CONC,CONC2,QSTEP,PRTOPT,
     &                    NSOLUTE,LAMBDA,CHEM,NFLOW,QINDEX,QWT,STOPTYPE,
     &                    JBOUND,TIMEB,TSTOP,QSTOP,DFACE,HPLUSF,HPLUSB,
     &                    HMULTF,HMULTB,HDIV,HDIVF,HADV,GAMMA,BTERMS,
     &                    AFACE,A,C,AWORK,BWORK,LAM2DT,DSDIST,USDIST,
     &                    QINVAL,CLVAL,AVALUE,QVALUE,IBOUND,USBC,NBOUND,
     &                    BCSTOP,SGROUP2,SGROUP,LHATDT,SORB,LHAT2DT,KD,
     &                    ISORB,CLATINN,QN,AREAN,QLATINN,GAMMAN,AFACEN,
     &                    AN,CN,BTERMSN,TWOPLUS,AREA2,ALPHA,BN,TGROUP,
     &                    IGROUP,
     &  IDECAY, XSTART, QSTART,PRTLOC,DISP,FLOWLOC,QLATOUT,DIST,
     &  LAMBDA2, LAMHAT,LAMHAT2,RHOLAM,CSBACK,
     &  nreach,lastseg, nseg, rchlen,
     &  new_area, new_disp, new_area2, new_alpha, nseg_tmp, rchlen_tmp)

      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
*     subroutine arguments
*
      INTEGER*4 IMAX,NPRINT,IPRINT,NSOLUTE,PRTOPT,NFLOW,JBOUND,IBOUND,
     &          NBOUND,ISORB
      INTEGER*4 PINDEX(*),QINDEX(*)
      DOUBLE PRECISION DSBOUND,TSTART,TFINAL,TSTEP,QSTEP,TIMEB,TSTOP,
     &                 QSTOP,BCSTOP
      DOUBLE PRECISION DELTAX(*),Q(*),USTIME(*),AREA(*),QLATIN(*),WT(*),
     &                 QWT(*),DSDIST(*),USDIST(*),QINVAL(*),DFACE(*),
     &                 HPLUSF(*),HPLUSB(*),HMULTF(*),HMULTB(*),HDIV(*),
     &                 HDIVF(*),HADV(*),GAMMA(*),AFACE(*),A(*),C(*),
     &                 AVALUE(*),QVALUE(*),USCONCN(*),QN(*),AREAN(*),
     &                 QLATINN(*),AFACEN(*),GAMMAN(*),AN(*),CN(*),
     &                 AREA2(*),ALPHA(*)
      DOUBLE PRECISION USCONC(MAXBOUND,*),CLATIN(MAXSEG,*),
     &                 CONC(MAXSEG,*),CONC2(MAXSEG,*),LAMBDA(MAXSEG,*),
     &                 LAM2DT(MAXSEG,*),CLVAL(MAXFLOWLOC+1,*),
     &                 AWORK(MAXSEG,*),BWORK(MAXSEG,*),USBC(MAXBOUND,*),
     &                 SGROUP2(MAXSEG,*),SGROUP(MAXSEG,*),
     &                 LHATDT(MAXSEG,*),SORB(MAXSEG,*),
     &                 LHAT2DT(MAXSEG,*),KD(MAXSEG,*),CLATINN(MAXSEG,*),
     &                 TWOPLUS(MAXSEG,*),BN(MAXSEG,*),TGROUP(MAXSEG,*),
     &                 BTERMS(MAXSEG,*),BTERMSN(MAXSEG,*),
     &                 IGROUP(MAXSEG,*)
      CHARACTER*(*) CHEM,STOPTYPE
*
*     local variables
*     sm: I have made the local variables subroutine arguments to pass them upwards
*
      INTEGER*4 IDECAY
      DOUBLE PRECISION XSTART,QSTART
      DOUBLE PRECISION PRTLOC(MAXPRINT),DISP(MAXSEG),
     &                 FLOWLOC(MAXFLOWLOC),QLATOUT(MAXSEG),DIST(MAXSEG)
      DOUBLE PRECISION LAMBDA2(MAXSEG,MAXSOLUTE),
     &                 LAMHAT(MAXSEG,MAXSOLUTE),
     &                 LAMHAT2(MAXSEG,MAXSOLUTE),
     &                 RHOLAM(MAXSEG,MAXSOLUTE),
     &                 CSBACK(MAXSEG,MAXSOLUTE)

*************************************************************************
* SM: new variables intorduced:
      integer*4 nreach,lastseg(0:MAXREACH),nseg_tmp, nseg(*)
      double precision new_area, new_disp, new_area2, new_alpha,
     &                 rchlen_tmp, rchlen(*)
*************************************************************************


* SM: I added the variables nreach, lastseg to the paramlist to access them from here
        ! this will change, DISP, AREA2, ALPHA
      !write(*,*) "maininit_mod(): nreach =", nreach
      CALL input2_mod(IMAX,DISP,AREA2,ALPHA,XSTART,DELTAX,NREACH,
     &                  LASTSEG, new_disp, new_area2, new_alpha,
     &                  nseg_tmp, rchlen_tmp)

      CALL change_area(new_area, AREA, nreach, lastseg)

      !write(*,*) "maininit_mod(): loc 10000. Prgram crashes after this:"
      !SM: -- call this every time A, D, As, Alpha is changed
      CALL PREPROC(IMAX,NBOUND,AREA2,DELTAX,Q,USTIME,AREA,ALPHA,DISP,
     &             QLATIN,TSTART,TSTEP,QSTEP,NSOLUTE,LAMBDA,LAMBDA2,
     &             CHEM,STOPTYPE,JBOUND,TIMEB,TSTOP,QSTOP,DFACE,HPLUSF,
     &             HPLUSB,HMULTF,HMULTB,HDIV,HDIVF,HADV,GAMMA,BTERMS,
     &             AFACE,A,C,AWORK,BWORK,LAM2DT,BCSTOP,IBOUND,RHOLAM,KD,
     &             LAMHAT,LAMHAT2,LHAT2DT,TWOPLUS,LHATDT,SGROUP,SGROUP2,
     &             CSBACK,BN,TGROUP,IGROUP,CLATIN)

      !write(*,*) "maininit_mod(): loc 10001"

*
*     save initial values of flow variables and preprocess
*
* SM: My work at VT is not steady state
      IF (CHEM .NE. 'Steady-State')
     &   CALL SAVEN(IMAX,NSOLUTE,AFACE,AFACEN,AREA,AREAN,GAMMA,GAMMAN,Q,
     &              QLATIN,QLATINN,QN,CLATIN,CLATINN,AN,CN,A,C,BTERMS,
     &              BTERMSN,TIMEB,TGROUP,BN,TWOPLUS,ALPHA,IGROUP)

      !write(*,*) "maininit_mod(): loc 10002"
*
*     Initialize concentrations
*
      !Masud: resetting the concentration variables to zero
      do 873 iseg = 1,MAXSEG
        do 874 isolute = 1,MAXSOLUTE
            CONC(iseg,isolute) = 0.D0
            CONC2(iseg,isolute) = 0.D0
874     continue
873   continue


      !write(*,*) "maininit_mod.f --> CONC2:", CONC(160,1),CONC2(160,1)
      !write(*,*) "maininit_mod.f --> USCONC:", USCONC(1,1), USCONC(2,1)
      CALL STEADY(IMAX,AREA2,DELTAX,Q,USCONC,AREA,ALPHA,DSBOUND,QLATIN,
     &            CLATIN,CONC,CONC2,NSOLUTE,LAMBDA,LAMBDA2,DFACE,HPLUSF,
     &            HPLUSB,HMULTF,HMULTB,HDIV,HDIVF,HADV,LAMHAT2,CSBACK,
     &            KD,SORB)
      !write(*,*) "maininit_mod.f --> CONC2**:",CONC(160,1),CONC2(160,1)
      !write(*,*) "maininit_mod(): loc 10003"

      !write(*,*) "Successful END of mainint_mod()!"

      RETURN
      END
