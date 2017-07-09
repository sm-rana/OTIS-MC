************************************************************************
*
*     Subroutine MAINRUN        Called by: MAIN Program
*
*     call the dynamic simulation routine and output results.
*
************************************************************************
       SUBROUTINE mainrun_mod(IMAX,NPRINT,PINDEX,WT,IPRINT,DELTAX,Q,
     &                   USTIME,
     &                   USCONC,USCONCN,AREA,DSBOUND,QLATIN,CLATIN,TIME,
     &                   TFINAL,TSTEP,CONC,CONC2,QSTEP,PRTOPT,NSOLUTE,
     &                   LAMBDA,CHEM,NFLOW,QINDEX,QWT,STOPTYPE,JBOUND,
     &                   TIMEB,TSTOP,QSTOP,DFACE,HPLUSF,HPLUSB,HMULTF,
     &                   HMULTB,HDIV,HDIVF,HADV,GAMMA,BTERMS,AFACE,A,C,
     &                   AWORK,BWORK,LAM2DT,DSDIST,USDIST,QINVAL,CLVAL,
     &                   AVALUE,QVALUE,IBOUND,USBC,NBOUND,BCSTOP,
     &                   SGROUP2,SGROUP,LHATDT,SORB,LHAT2DT,KD,ISORB,
     &                   CLATINN,QN,AREAN,QLATINN,GAMMAN,AFACEN,AN,CN,
     &                   BTERMSN,TWOPLUS,AREA2,ALPHA,BN,TGROUP,IGROUP,
     &                   sim_time, conc_main, conc_ts)
      INCLUDE 'fmodules.inc'
*
*     subroutine arguments
*
      INTEGER*4 IMAX,NPRINT,IPRINT,NSOLUTE,PRTOPT,NFLOW,JBOUND,IBOUND,
     &          NBOUND,ISORB
      INTEGER*4 PINDEX(*),QINDEX(*)
      DOUBLE PRECISION DSBOUND,TIME,TFINAL,TSTEP,QSTEP,TIMEB,TSTOP,
     &                 QSTOP,BCSTOP
      DOUBLE PRECISION DELTAX(*),Q(*),USTIME(*),AREA(*),QLATIN(*),WT(*),
     &                 QWT(*),DSDIST(*),USDIST(*),QINVAL(*),DFACE(*),
     &                 HPLUSF(*),HPLUSB(*),HMULTF(*),HMULTB(*),HDIV(*),
     &                 HDIVF(*),HADV(*),GAMMA(*),AFACE(*),A(*),C(*),
     &                 LAM2DT(*),AVALUE(*),QVALUE(*),USCONCN(*),QN(*),
     &                 AREAN(*),QLATINN(*),AFACEN(*),GAMMAN(*),AN(*),
     &                 CN(*),AREA2(*),ALPHA(*)
      DOUBLE PRECISION USCONC(MAXBOUND,*),CLATIN(MAXSEG,*),
     &                 CONC(MAXSEG,*),CONC2(MAXSEG,*),
     &                 LAMBDA(MAXSEG,*),CLVAL(MAXFLOWLOC+1,*),
     &                 AWORK(MAXSEG,*),BWORK(MAXSEG,*),
     &                 USBC(MAXBOUND,*),SGROUP2(MAXSEG,*),
     &                 SGROUP(MAXSEG,*),LHATDT(MAXSEG,*),SORB(MAXSEG,*),
     &                 LHAT2DT(MAXSEG,*),KD(MAXSEG,*),CLATINN(MAXSEG,*),
     &                 TWOPLUS(MAXSEG,*),BN(MAXSEG,*),TGROUP(MAXSEG,*),
     &                 BTERMS(MAXSEG,*),BTERMSN(MAXSEG,*),
     &                 IGROUP(MAXSEG,*)
      CHARACTER*(*) CHEM,STOPTYPE

      integer print_index, num_solute
      double precision sim_time(MAXBOUND), conc_main(MAXBOUND),
     & conc_ts(MAXBOUND)
*
*     local variables
*
      INTEGER*4 J,NCALLS
*
*     if this is a steady-state run simply return
*
      IF (CHEM .EQ. 'Steady-State') RETURN
*
*     compute the required number of calls to the dynamic solution
*     routine.  Each call advances the solution (in time) to the next
*     output point.
*
      NCALLS = INT((INT((TFINAL - TIME)/TSTEP) + 1) / IPRINT) + 1
      !SM: each call to dynamic updates TIME
*
*     begin the time loop
*
      !write(*,*) "mainrun_mod.f--> Area, Disp, As, Alp"
      !write(*,*) AREA(1), DFACE(1), AREA2(1), ALPHA(1)
      !masud:
      print_index = 1
      num_solute = 1
      !write(*,*) "Printing from inside mainrun_mod.f -->"
      DO 20 J = 1, NCALLS

      CALL DYNAMIC(IMAX,IPRINT,DELTAX,Q,USTIME,USCONC,USCONCN,AREA,
     &                DSBOUND,QLATIN,CLATIN,TIME,TSTEP,CONC,CONC2,QSTEP,
     &                NSOLUTE,LAMBDA,QVALUE,AVALUE,QWT,QINDEX,NFLOW,
     &                QINVAL,CLVAL,CHEM,DSDIST,USDIST,STOPTYPE,JBOUND,
     &                TIMEB,TSTOP,QSTOP,DFACE,HPLUSF,HPLUSB,HMULTF,
     &                HMULTB,HDIV,HDIVF,HADV,GAMMA,BTERMS,AFACE,A,C,
     &                AWORK,BWORK,LAM2DT,IBOUND,USBC,NBOUND,BCSTOP,
     &                SGROUP2,SGROUP,LHATDT,SORB,LHAT2DT,KD,CLATINN,QN,
     &                AREAN,QLATINN,GAMMAN,AFACEN,AN,CN,BTERMSN,TWOPLUS,
     &                AREA2,ALPHA,BN,TGROUP,IGROUP)


         !   write(*,*) "C1, C2 : ",
         !& CONC(PINDEX(print_index),num_solute),
         !& CONC2(PINDEX(print_index),num_solute)

      sim_time(J) = TIME
      conc_main(J) = CONC(PINDEX(print_index),num_solute)
      conc_ts(J) = CONC2(PINDEX(print_index),num_solute)

 20   CONTINUE

       !write(*,*) "Successful END of mainrun_mod()!"
      RETURN
      END


