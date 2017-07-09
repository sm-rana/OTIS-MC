***********************************************************************
*
*     Subroutine DYNAMIC             Called by: MAINRUN
*
*     Simulate dynamic solute transport
*
***********************************************************************
      SUBROUTINE DYNAMIC(IMAX,IPRINT,DELTAX,Q,USTIME,USCONC,USCONCN,
     &                   AREA,DSBOUND,QLATIN,CLATIN,TIME,TSTEP,CONC,
     &                   CONC2,QSTEP,NSOLUTE,LAMBDA,QVALUE,AVALUE,QWT,
     &                   QINDEX,NFLOW,QINVAL,CLVAL,CHEM,DSDIST,USDIST,
     &                   STOPTYPE,JBOUND,TIMEB,TSTOP,QSTOP,DFACE,HPLUSF,
     &                   HPLUSB,HMULTF,HMULTB,HDIV,HDIVF,HADV,GAMMA,
     &                   BTERMS,AFACE,A,C,AWORK,BWORK,LAM2DT,IBOUND,
     &                   USBC,NBOUND,BCSTOP,SGROUP2,SGROUP,LHATDT,SORB,
     &                   LHAT2DT,KD,CLATINN,QN,AREAN,QLATINN,GAMMAN,
     &                   AFACEN,AN,CN,BTERMSN,TWOPLUS,AREA2,ALPHA,BN,
     &                   TGROUP,IGROUP)
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
*     subroutine arguments
*
      CHARACTER*(*) CHEM,STOPTYPE
      INTEGER*4 IMAX,IPRINT,NSOLUTE,NFLOW,JBOUND,IBOUND,NBOUND,QINDEX(*)
      DOUBLE PRECISION DSBOUND,TIME,TSTEP,QSTEP,TIMEB,TSTOP,QSTOP,BCSTOP
      DOUBLE PRECISION DELTAX(*),Q(*),USTIME(*),AREA(*),QLATIN(*),
     &                 QWT(*),QINVAL(*),DSDIST(*),USDIST(*),DFACE(*),
     &                 HPLUSF(*),HPLUSB(*),HMULTF(*),HMULTB(*),HDIV(*),
     &                 HDIVF(*),HADV(*),GAMMA(*),AFACE(*),A(*),C(*),
     &                 QVALUE(*),AVALUE(*),USCONCN(*),QN(*),AREAN(*),
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
*
*     local variables
*
      INTEGER*4 J
*
*     format statement
*
 100  FORMAT(/,10X,'WARNING:'
     &       /,10X,'Your request to change the boundary condition or ',
     &       /,10X,'the flow variables at ',1PE11.5,' hours can not',
     &       /,10X,'be met because it is not aligned with the time ',
     &       /,10X,'step.  The requested change is effective at ',
     &       1PE11.5,' hours.')
*
*     begin the time loop
*
      DO 20 J = 1, IPRINT
*
*        increment TIME, and check to see if the boundary condition
*        and/or flow variables need updating.  The IF statement executes
*        a `dowhile TIME > TSTOP'
*
            ! Masud: TIME starts with TSTART and then it is updated
            ! Masud: TSTOP has been initialized to be TSTART (i.e., = 0)
         TIME = TIME + TSTEP
 10      IF (TIME .GT. TSTOP+1.D-7) THEN
*
*
*           Update the boundary conditions and/or the flow variables.
*
*           If the change occurs w/i the integration inverval (i.e.
*           TIME-TSTEP < TSTOP < TIME), the change will be effective
*           immediatly (at TIME-TSTEP rather than at TSTOP); print a
*           warning message to this effect.  Note that the user can
*           avoid this problem by making the time step compatible with
*           the boundary conditions and the flow routing.
*
            IF (ABS(TSTOP-TIME+TSTEP).GT. 1.D-5)
     &         WRITE(LDECHO,100) TSTOP,TIME-TSTEP
            CALL UPDATE(IMAX,Q,AREA,QLATIN,CLATIN,NSOLUTE,QVALUE,AVALUE,
     &                  QWT,QINDEX,NFLOW,JBOUND,STOPTYPE,QSTOP,QSTEP,
     &                  USTIME,TSTOP,TIMEB,DELTAX,DFACE,HPLUSF,HPLUSB,
     &                  HMULTF,HMULTB,HDIV,HDIVF,HADV,GAMMA,BTERMS,
     &                  AFACE,A,C,AWORK,BWORK,QINVAL,CLVAL,DSDIST,
     &                  USDIST,LAMBDA,LAM2DT,IBOUND,USBC,NBOUND,USCONC,
     &                  BCSTOP,TSTEP,TIME,LHAT2DT,SGROUP,KD,TWOPLUS,
     &                  AREA2,ALPHA,BN,BTERMSN,TGROUP,GAMMAN,IGROUP,
     &                  QLATINN,CLATINN,AREAN)
            GOTO 10
         ENDIF
*
*        compute concentrations
*
         IF (CHEM .EQ. 'Conservative') THEN
            CALL CONSER(IMAX,JBOUND,DELTAX,Q,USCONC,USCONCN,AREA,DFACE,
     &                  DSBOUND,CONC,CONC2,NSOLUTE,GAMMA,AFACE,AWORK,
     &                  BWORK,QN,AREAN,GAMMAN,AFACEN,AN,CN,C,TWOPLUS,BN,
     &                  TGROUP,IGROUP)
         ELSE
            CALL REACT(IMAX,JBOUND,DELTAX,Q,USCONC,USCONCN,AREA,DFACE,
     &                 DSBOUND,CONC,CONC2,NSOLUTE,GAMMA,AFACE,AWORK,
     &                 BWORK,LAM2DT,SGROUP2,SGROUP,LHATDT,SORB,LHAT2DT,
     &                 KD,QN,AREAN,GAMMAN,AFACEN,AN,CN,C,TWOPLUS,ALPHA,
     &                 BN,TGROUP,IGROUP)
         ENDIF
*
*        for unsteady flow, save various parameters at time step 'N' and
*        preprocess
*
         IF (QSTEP .NE. 0.D0)
     &      CALL SAVEN(IMAX,NSOLUTE,AFACE,AFACEN,AREA,AREAN,GAMMA,
     &                 GAMMAN,Q,QLATIN,QLATINN,QN,CLATIN,CLATINN,AN,CN,
     &                 A,C,BTERMS,BTERMSN,TIMEB,TGROUP,BN,TWOPLUS,ALPHA,
     &                 IGROUP)

 20   CONTINUE

      RETURN
      END

