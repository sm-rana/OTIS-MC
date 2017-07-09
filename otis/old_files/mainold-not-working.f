************************************************************************
*
*     OTIS: One-dimensional Transport with Inflow and Storage
*     -------------------------------------------------------
*
*     OTIS is a mathematical simulation model used to characterize the
*     fate and transport of water-borne solutes in streams and rivers.
*     The governing equation underlying the model is the advection-
*     dispersion equation with additional terms to account for transient
*     storage, lateral inflow, first-order decay and sorption. This
*     equation and the associated equations describing transient
*     storage and sorption are solved using a Crank-Nicolson finite
*     difference solution.
*
*     The OTIS solute transport model, data, and documentation are made
*     available by the U.S. Geological Survey (USGS) to be used in the
*     public interest and the advancement of science. You may, without
*     any fee or cost, use, copy, modify, or distribute this software,
*     and any derivative works thereof, and its supporting
*     documentation, SUBJECT TO the USGS software User's Rights Notice
*     (http://water.usgs.gov/software/software_notice.html)
*
*
************************************************************************
*
*     Version: MOD40         Feb. 1998
*
************************************************************************
*
*     Reference
*     ---------
*     Runkel, R.L., 1998, One dimensional transport with inflow and
*       storage (OTIS): A solute transport model for streams and rivers:
*       U.S. Geological Survey Water-Resources Investigation Report
*       98-4018. xx p.
*
*     Homepage:  http://webserver.cr.usgs.gov/otis/home.html
*     ---------
*
*     Author
*     ------
*     R.L. Runkel
*     U.S. Geological Survey
*     Mail Stop 415
*     Denver Federal Center
*     Lakewood, CO 80225
*     Internet: runkel@usgs.gov
*
************************************************************************
*
*                            SEGMENTATION
*
************************************************************************
*
*
*
*     |<--- Hi-1 --->|<---- Hi ---->|<--- Hi+1 --->|
*     ______________________________________________
*     |              |              |              |
*     |     Ci-1     |      Ci      |     Ci+1     |
*     |     Qi-1     |      Qi      |     Qi+1     |
*     |     Ai-1     |      Ai      |     Ai+1     |
*     |              |              |              |
*     ----------------------------------------------
*                  DFACEi         DFACEi+1
*                  AFACEi         AFACEi+1
*
*     where:
*
*           A = AREA
*           C = CONC, CONC2, or SORB
*           H = DELTAX
*
*     as defined below.
*
*
************************************************************************
*
*                      DICTIONARY - INPUT VARIABLES
*
************************************************************************
*
*     Required Input Data
*     -------------------
*     TITLE      simulation title
*     PRTOPT     print option
*     TSTEP      time step interval [hr]
*     PSTEP      print step interval [hr]
*     TSTART     starting time [hr]
*     TFINAL     final time [hr]
*     XSTART     starting distance at upstream boundary [L]
*     DSBOUND    downstream boundary condition (flux) [mass/L^2-s]
*     NREACH     number of reaches
*
*     Required Reach Parameters (One value for each reach)
*     ----------------------------------------------------
*     NSEG       number of segments
*     RCHLEN     length of reach [L]
*     DISP       dispersion coefficient [L^2/s]
*     AREA2      storage zone cross-sectional area [L^2]
*     ALPHA      exchange coefficient [/s]
*
*     Required Print Information
*     --------------------------
*     NPRINT     number of print locations
*     IOPT       interpolation option
*     PRTLOC     distance at which results are printed [L]
*
*     Required Chemical Information
*     -----------------------------
*     NSOLUTE    number of solutes to simulate
*     IDECAY     decay flag (=0 no decay; =1 first-order decay)
*     ISORB      sorption flag (=0 no sorption; =1 Kd sorption)
*
*     First-Order Decay Rate Information (Required for IDECAY=1)
*     ----------------------------------------------------------
*     LAMBDA     decay coefficient for the main channel [/s]
*     LAMBDA2    decay coefficient for the storage zone [/s]
*
*     Kd Sorption Information (Required for ISORB=1)
*     ----------------------------------------------
*     LAMHAT     sorption rate coefficient for the main channel [/sec]
*     LAMHAT2    sorption rate coefficient for the storage zone [/sec]
*     RHO        mass of accessibile sediment/volume water [mass/L^3]
*     KD         distribution coefficient [L^3/mass]
*     CSBACK     background storage zone solute concentration [mass/L^3]
*
*     Required Chemical and Upstream Boundary Condition Information
*     -------------------------------------------------------------
*     NBOUND     number of different upstream boundary conditions
*     IBOUND     boundary condition option
*     USTIME     time at which upstream boundary condition goes
*                into effect [hr]
*     USBC       upstream boundary condition
*
*     Required Flow and Area Information
*     ----------------------------------
*     QSTEP      time interval between changes in flow (0=steady flow)
*     QSTART     volumetric flowrate at upstream boundary [L^3/s]
*     QLATIN     lateral inflow rate [L^3/s-L]
*     QLATOUT    lateral outflow rate [L^3/s-L]
*     AREA       main channel cross-sectional area [L^2]
*     CLATIN     concentration of lateral inflow [mass/L^3]
*
*     Information for Unsteady Flow Regimes (QSTEP > 0)
*     -------------------------------------------------
*     NFLOW      number of locations at which Q and AREA are specified
*     FLOWLOC    distance at which Q and AREA are specified [L]
*     QVALUE     flowrates at specified distances
*     QINVAL     lateral inflow rate at specified distances
*     AVALUE     areas at specified distances
*     CLVAL      concentration of lateral inflows at specified distances
*
************************************************************************
*
*                      DICTIONARY - PROGRAM VARIABLES
*
************************************************************************
*
*     State Variables
*     ---------------
*     CONC       concentration in main channel [mass/L^3]
*     CONC2      concentration in storage zone [mass/L^3]
*     SORB       streambed sediment concentration [mass/mass]
*
*     System definition and hydrology
*     -------------------------------
*     AFACE      cross sectional area @ interface of segments I,I+1
*     DFACE      dispersion coefficient @ interface of segments I,I+1
*     DELTAX     segment length [L]
*     IMAX       number of segments
*     LASTSEG    last segment in each reach
*     DIST       distance corresponding to segment centroid
*     Q          volumetric flowrate [L^3/s]
*     QINDEX     flow location used for interpolation (when QSTEP > 0)
*     QWT        weight used to interpolate between QINDEX and QINDEX+1
*     DSDIST     distance to the nearest downstream flow location/DELTAX
*     USDIST     distance to the nearest upstream flow location / DELTAX
*
*     Matrix Solution
*     ---------------
*     A          upper diagonal of the tridiagonal matrix
*     B          main diagonal of the tridiagonal matrix
*     C          lower diagonal of the tridiagonal matrix
*     D          constant vector
*     AWORK      work vector corresponding to A
*     BWORK      work vector corresponding to B
*
*     Parameter Groups - Time invariant (H = DELTAX)
*     ----------------------------------------------
*     HPLUSF     H(i) + H(i+1)
*     HPLUSB     H(i) + H(i-1)
*     HMULTF     H(i) [ H(i) + H(i+1) ]
*     HMULTB     H(i) [ H(i) + H(i-1) ]
*     HDIV       H(i) / [ H(i) + H(i+1) ]
*     HDIVF      H(i+1) / [ H(i) + H(i+1) ]
*     HADV       0.5/H(i) { H(i+1)/[H(i+1)+H(i)] - H(i-1)/[H(i-1)+H(i)]}
*     LAM2DT     storage zone decay coefficient times the timestep
*     LHATDT     main channel sorption rate times the timestep
*     LHAT2DT    storage zone sorption rate times the timestep
*     RHOLAM     mass of access. sediment times the m.c. sorption rate
*     SGROUP     sorption parameter group
*     SGROUP2    sorption parameter group
*
*     Parameter Groups - Time invariant given steady flow
*     ---------------------------------------------------
*     BN         group related to the main diagonal vector, B
*     BTERMS     group related to the main diagonal vector, B
*     GAMMA      group introduced by storage zone equation
*     IGROUP     inflow group
*     TGROUP     transient storage group
*     TWOPLUS    2 + .....
*
*     Boundary Conditions & Flow Changes
*     ----------------------------------
*     BCSTOP     time at which the boundary conditions change
*     JBOUND     index of the current boundary condition
*     STOPTYPE   indicates whether the next TSTOP is due to a change
*                in the boundary conditions or the flow variables.
*     TSTOP      time at which the boundary condition or flow variables
*                change
*     USCONC     concentration at the upstream boundary
*     QSTOP      time at which the flow variables change
*
*     Program Variables - Misc.
*     -------------------------
*     TIMEB      inverse of the time step [/sec]
*     TIME       time [hr]
*     CHEM       type of chemistry considered (reactive or conservative)
*     IPRINT     number of iterations between printing of results
*     PINDEX     segments for which results are printed
*     WT         weight used to interpolate between PINDEX and PINDEX+1
*
*     Program Variables - time level 'N'
*     ----------------------------------
*     The following variables contain values of associated with the
*     previous timestep, time level 'N'.  Their counterparts at time
*     level 'N+1' are defined above:
*
*     AFACEN,AN,AREAN,BTERMSN,CLATINN,CN,GAMMAN,QLATINN,QN,USCONCN
*
*     i.e. CLATINN is the lateral inflow concentration at time N, while
*     CLATIN is corresponding concentration at time N+1.
*
************************************************************************
*
*                      INCLUDE FILES
*
*************************************************************************
*
*     Maximum Dimensions (set via PARAMETER statements in fmodules.inc)
*     -----------------------------------------------------------------
*     MAXSEG     =5000  maximum number of segments
*     MAXBOUND   =20000 maximum number of upstream boundary conditions
*     MAXREACH   =30    maximum number of reaches
*     MAXPRINT   =30    maximum number of print locations
*     MAXSOLUTE  =3     maximum number of solutes
*     MAXFLOWLOC =30    maximum number of flow locations (unsteady flow)
*
*     Logical Devices (as defined in lda.inc)
*     ---------------------------------------
*     LDCTRL     input control information (I/O filenames)
*     LDPARAM    input file for chemical and physical parameters
*     LDFLOW     input file for flow, lateral inflow, and areas
*     LDECHO     output data and time, echo input parameters
*     LDFILES    output files, one file per solute
*     LDSORB     sorption output files (ISORB=1), one file per solute
*
************************************************************************
      !SUBROUTINE MAINOLD(mcmc_iter, new_area, new_disp, new_area2,
      !&   new_alpha,nseg_tmp,rchlen_tmp)
      !SUBROUTINE MAINOLD(mcmc_iter, new_area, new_disp, new_area2,
      !&   new_alpha,nreach,nseg_tmp,rchlen_tmp, time_start)

*      SUBROUTINE MAINOLD(IMAX, NPRINT, IPRINT, NSOLUTE, PRTOPT, NFLOW,
*     & JBOUND, IBOUND, NBOUND, ISORB, PINDEX, QINDEX,
*     & DSBOUND, TSTART, TFINAL, TSTEP, QSTEP,
*     & TIMEB, TSTOP, QSTOP, BCSTOP,
*     & DELTAX, Q, USTIME, AREA,
*     & QLATIN, WT, QVALUE, AVALUE,
*     & QWT, DSDIST, USDIST, DFACE, HPLUSF,
*     & HPLUSB, HMULTF, HMULTB, HDIV,
*     & HDIVF, HADV, GAMMA, AFACE, A, C,
*     & QINVAL, USCONCN, QN, AREAN, QLATINN,
*     & AFACEN, GAMMAN, AN, CN,
*     & AREA2, ALPHA,
*     & USCONC, CLATIN, CONC, CONC2,
*     & CLVAL, LAMBDA, LAM2DT, AWORK,
*     & BWORK, USBC,     SGROUP2, SGROUP,
*     & LHATDT,  SORB,     LHAT2DT, KD,
*     & CLATINN, TWOPLUS,  BN, TGROUP,
*     & BTERMS,  BTERMSN,  IGROUP,
*     & CHEM, STOPTYPE,
*     & mcmc_iter, new_area, new_disp, new_area2,
*     & new_alpha, nreach, nseg_tmp, rchlen_tmp, time_start)

      SUBROUTINE MAINOLD(mcmc_iter, new_area, new_disp, new_area2,
     & new_alpha, nreach, nseg_tmp, rchlen_tmp, time_start)
*     PROGRAM MAIN
*
*     dimensional parameters
*
      INCLUDE 'fmodules.inc'
*
*     other variables
*
*************************************************************************************
      integer*4 mcmc_iter, nreach, nseg_tmp
      double precision new_area, new_disp, new_area2, new_alpha
      double precision rchlen_tmp, time_start

*************************************************************************************
* Local variables ***********
************************************************
      !SM: NSEG and RCHLEN gets initialized in input2.f need to trickle it up all the way to main()
      integer*4 lastseg(0:MAXREACH), nseg(MAXREACH)
      double precision rchlen(MAXREACH)
*************************************************************************************

      INTEGER*4 IMAX,NPRINT,IPRINT,NSOLUTE,PRTOPT,NFLOW,JBOUND,IBOUND,
     &          NBOUND,ISORB
      INTEGER*4 PINDEX(MAXPRINT),QINDEX(MAXSEG)
      DOUBLE PRECISION DSBOUND,TSTART,TFINAL,TSTEP,QSTEP,TIMEB,TSTOP,
     &                 QSTOP,BCSTOP
      DOUBLE PRECISION DELTAX(MAXSEG),Q(MAXSEG),USTIME(MAXBOUND+1),
     &                 AREA(MAXSEG),QLATIN(MAXSEG),WT(MAXPRINT),
     &                 QVALUE(MAXFLOWLOC),AVALUE(MAXFLOWLOC),
     &                 QWT(MAXSEG),DSDIST(MAXSEG),USDIST(MAXSEG),
     &                 DFACE(MAXSEG),HPLUSF(MAXSEG),HPLUSB(MAXSEG),
     &                 HMULTF(MAXSEG),HMULTB(MAXSEG),HDIV(MAXSEG),
     &                 HDIVF(MAXSEG),HADV(MAXSEG),GAMMA(MAXSEG),
     &                 AFACE(MAXSEG),A(MAXSEG),C(MAXSEG),
     &                 QINVAL(MAXFLOWLOC+1),USCONCN(MAXSOLUTE),QN(2),
     &                 AREAN(MAXSEG),QLATINN(MAXSEG),AFACEN(2),
     &                 GAMMAN(MAXSEG),AN(MAXSEG),CN(MAXSEG),
     &                 AREA2(MAXSEG),ALPHA(MAXSEG)

      DOUBLE PRECISION USCONC(MAXBOUND,MAXSOLUTE),
     &                 CLATIN(MAXSEG,MAXSOLUTE),CONC(MAXSEG,MAXSOLUTE),
     &                 CONC2(MAXSEG,MAXSOLUTE),
     &                 CLVAL(MAXFLOWLOC+1,MAXSOLUTE),
     &                 LAMBDA(MAXSEG,MAXSOLUTE),
     &                 LAM2DT(MAXSEG,MAXSOLUTE),
     &                 AWORK(MAXSEG,MAXSOLUTE),BWORK(MAXSEG,MAXSOLUTE),
     &                 USBC(MAXBOUND,MAXSOLUTE),
     &                 SGROUP2(MAXSEG,MAXSOLUTE),
     &                 SGROUP(MAXSEG,MAXSOLUTE),
     &                 LHATDT(MAXSEG,MAXSOLUTE),SORB(MAXSEG,MAXSOLUTE),
     &                 LHAT2DT(MAXSEG,MAXSOLUTE),KD(MAXSEG,MAXSOLUTE),
     &                 CLATINN(MAXSEG,MAXSOLUTE),
     &                 TWOPLUS(MAXSEG,MAXSOLUTE),BN(MAXSEG,MAXSOLUTE),
     &                 TGROUP(MAXSEG,MAXSOLUTE),
     &                 BTERMS(MAXSEG,MAXSOLUTE),
     &                 BTERMSN(MAXSEG,MAXSOLUTE),
     &                 IGROUP(MAXSEG,MAXSOLUTE)
      CHARACTER*12 CHEM,STOPTYPE
*
*     Obtain input parameters, initialize variables, and preprocess.
*     If the steady-state option has been invoked, output results
*
* SM: Need to call MAININIT() multiple times to run MCMC simulations with many
* of parameters:
*
      ! Run the original program for the first iteration:

      ! ----------------- IF (mcmc_iter .EQ. 0) BLOCK ----------------------------

      IF (mcmc_iter.EQ.0) THEN
      write(*,*) "In MAINOLD() --> Passed Variables:"
      write(*,*) "A, D, As, Al:", new_area, new_disp, new_area2,
     & new_alpha

      write(*,*) "In MAINOLD() --> Calling MAININIT():"
      ! DEBUG:
      ! nseg(1) = 500
      ! rchlen(1) = 80

      CALL MAININIT(IMAX,NPRINT,PINDEX,WT,IPRINT,DELTAX,Q,USTIME,
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
     &                    IGROUP,nreach,lastseg, nseg, rchlen)
                ! SM: NSEG(MAXREACH), double PRECISION RCHLEN(MAXREACH)
                ! to be tricked up from input()
       write(*,*) "mainold() --> NREACH = ", nreach
       !SUCCESS: write(*,*) "DEBUG: nseg(1),rchlen(1) = ",nseg(1),rchlen(1)
*
*     perform the dynamic simulation
*
      time_start = TSTART ! SM: getting start time for the ELSE block
      nseg_tmp = nseg(1)
      rchlen_tmp = rchlen(1)
      !write(*,*) "Nseg, Rchlen = ", NSEG(1), RCHLEN(1) NREACH

      write(*,*) "Calling MAINRUN()"
      CALL MAINRUN(IMAX,NPRINT,PINDEX,WT,IPRINT,DELTAX,Q,USTIME,USCONC,
     &             USCONCN,AREA,DSBOUND,QLATIN,CLATIN,TSTART,TFINAL,
     &             TSTEP,CONC,CONC2,QSTEP,PRTOPT,NSOLUTE,LAMBDA,CHEM,
     &             NFLOW,QINDEX,QWT,STOPTYPE,JBOUND,TIMEB,TSTOP,QSTOP,
     &             DFACE,HPLUSF,HPLUSB,HMULTF,HMULTB,HDIV,HDIVF,HADV,
     &             GAMMA,BTERMS,AFACE,A,C,AWORK,BWORK,LAM2DT,DSDIST,
     &             USDIST,QINVAL,CLVAL,AVALUE,QVALUE,IBOUND,USBC,NBOUND,
     &             BCSTOP,SGROUP2,SGROUP,LHATDT,SORB,LHAT2DT,KD,ISORB,
     &             CLATINN,QN,AREAN,QLATINN,GAMMAN,AFACEN,AN,CN,BTERMSN,
     &             TWOPLUS,AREA2,ALPHA,BN,TGROUP,IGROUP)

      ELSE


       do 1986 mcmc_iter = 0, mcmc_len
      write(*,*) "Interation = ", mcmc_iter
1986   continue


      write(*,*) "Calling maininit_mod()"
       CALL maininit_mod(IMAX,DELTAX,Q,USTIME,AREA,DSBOUND,QLATIN,
     &                    CLATIN,TSTART,TSTEP,CONC,CONC2,QSTEP,
     &                    NSOLUTE,LAMBDA,CHEM,STOPTYPE,
     &                    JBOUND,TIMEB,TSTOP,QSTOP,DFACE,HPLUSF,HPLUSB,
     &                    HMULTF,HMULTB,HDIV,HDIVF,HADV,GAMMA,BTERMS,
     &                    AFACE,A,C,AWORK,BWORK,LAM2DT,
     &                    IBOUND,NBOUND,BCSTOP,SGROUP2,SGROUP,LHATDT,
     &                    SORB,LHAT2DT,KD,CLATINN,QN,AREAN,QLATINN,
     &                    GAMMAN,AFACEN, AN,CN,BTERMSN,TWOPLUS,AREA2,
     &                    ALPHA,BN,TGROUP, IGROUP, nreach, lastseg,
     &   new_area, new_disp, new_area2, new_alpha, nseg_tmp, rchlen_tmp)

      write(*,*) "Calling mainrun_mod()"
      CALL mainrun_mod(IMAX,NPRINT,PINDEX,WT,IPRINT,DELTAX,Q,
     &                   USTIME,
     &                   USCONC,USCONCN,AREA,DSBOUND,QLATIN,CLATIN,
     &                   time_start,
     &                   TFINAL,TSTEP,CONC,CONC2,QSTEP,PRTOPT,NSOLUTE,
     &                   LAMBDA,CHEM,NFLOW,QINDEX,QWT,STOPTYPE,JBOUND,
     &                   TIMEB,TSTOP,QSTOP,DFACE,HPLUSF,HPLUSB,HMULTF,
     &                   HMULTB,HDIV,HDIVF,HADV,GAMMA,BTERMS,AFACE,A,C,
     &                   AWORK,BWORK,LAM2DT,DSDIST,USDIST,QINVAL,CLVAL,
     &                   AVALUE,QVALUE,IBOUND,USBC,NBOUND,BCSTOP,
     &                   SGROUP2,SGROUP,LHATDT,SORB,LHAT2DT,KD,ISORB,
     &                   CLATINN,QN,AREAN,QLATINN,GAMMAN,AFACEN,AN,CN,
     &                   BTERMSN,TWOPLUS,AREA2,ALPHA,BN,TGROUP,IGROUP)
      ENDIF
      ! -----------END--- IF (mcmc_iter .EQ. 0) BLOCK ----------------------------




*
*     close files
*
      CALL CLOSEF(NSOLUTE,ISORB)


*      STOP
* SM: adding return instead of STOP
      RETURN
      END




