************************************************************************
*
*     Subroutine UPDATE             Called by: DYNAMIC, DYNAMIC2
*
*     Implement a change in the boundary conditions and/or the flow
*     variables.  This routine is called when:
*
*        (TIME - TSTEP)  =   TSTOP
*
*     where:
*
*        TSTOP           time of change
*        TIME            new time (time level N+1)
*        (TIME-TSTEP)    old time (time level N)
*
*     UPDATE the appropriate variables and reset TSTOP.
*
************************************************************************
      SUBROUTINE UPDATE(IMAX,Q,AREA,QLATIN,CLATIN,NSOLUTE,QVALUE,AVALUE,
     &                  QWT,QINDEX,NFLOW,JBOUND,STOPTYPE,QSTOP,QSTEP,
     &                  USTIME,TSTOP,TIMEB,DELTAX,DFACE,HPLUSF,HPLUSB,
     &                  HMULTF,HMULTB,HDIV,HDIVF,HADV,GAMMA,BTERMS,
     &                  AFACE,A,C,AWORK,BWORK,QINVAL,CLVAL,DSDIST,
     &                  USDIST,LAMBDA,LAM2DT,IBOUND,USBC,NBOUND,USCONC,
     &                  BCSTOP,TSTEP,TIME,LHAT2DT,SGROUP,KD,TWOPLUS,
     &                  AREA2,ALPHA,BN,BTERMSN,TGROUP,GAMMAN,IGROUP,
     &                  QLATINN,CLATINN,AREAN)
      INCLUDE 'fmodules.inc'
*
*     subroutine arguments
*
      CHARACTER*(*) STOPTYPE
      INTEGER*4 IMAX,NSOLUTE,NFLOW,JBOUND,IBOUND,NBOUND,QINDEX(*)
      DOUBLE PRECISION QSTOP,QSTEP,TSTOP,TIMEB,BCSTOP,TSTEP,TIME
      DOUBLE PRECISION Q(*),AREA(*),QLATIN(*),QVALUE(*),AVALUE(*),
     &                 QWT(*),USTIME(*),DELTAX(*),DFACE(*),HPLUSF(*),
     &                 HPLUSB(*),HMULTF(*),HMULTB(*),HDIV(*),HDIVF(*),
     &                 HADV(*),GAMMA(*),AFACE(*),A(*),C(*),QINVAL(*),
     &                 DSDIST(*),USDIST(*),AREA2(*),ALPHA(*),GAMMAN(*),
     &                 QLATINN(*),AREAN(*)
      DOUBLE PRECISION CLATIN(MAXSEG,*),CLVAL(MAXFLOWLOC+1,*),
     &                 LAMBDA(MAXSEG,*),LAM2DT(MAXSEG,*),
     &                 AWORK(MAXSEG,*),BWORK(MAXSEG,*),
     &                 USBC(MAXBOUND,*),USCONC(MAXBOUND,*),
     &                 LHAT2DT(MAXSEG,*),SGROUP(MAXSEG,*),KD(MAXSEG,*),
     &                 TWOPLUS(MAXSEG,*),BN(MAXSEG,*),TGROUP(MAXSEG,*),
     &                 BTERMS(MAXSEG,*),BTERMSN(MAXSEG,*),
     &                 IGROUP(MAXSEG,*),CLATINN(MAXSEG,*)
*
*     update the flow variables and/or the boundary condition
*
      IF (STOPTYPE .EQ. 'QChange') THEN
*
*        the flow has changed.  Read/Update flow variables and recompute
*        the matrix coefficients
*
         CALL QCHANGE(IMAX,Q,AREA,QLATIN,CLATIN,NSOLUTE,QVALUE,AVALUE,
     &                QWT,QINDEX,NFLOW,JBOUND,QSTOP,QSTEP,TSTOP,USCONC,
     &                QINVAL,CLVAL,DSDIST,USDIST,IBOUND,USBC,NBOUND)
         CALL PREPROC3(IMAX,DELTAX,Q,AREA,DFACE,QLATIN,HPLUSF,HPLUSB,
     &                 HMULTF,HMULTB,HDIV,HDIVF,HADV,GAMMA,BTERMS,AFACE,
     &                 A,C,AWORK,BWORK,TIMEB,NSOLUTE,LAM2DT,LHAT2DT,KD,
     &                 TWOPLUS,AREA2,ALPHA,LAMBDA,SGROUP,BTERMSN,BN,
     &                 TGROUP,GAMMAN,IGROUP,QLATINN,CLATINN,AREAN,
     &                 CLATIN)

      ELSEIF (STOPTYPE .EQ. 'BCondition') THEN
*
*        the boundary condition has changed
*
         CALL BCCHANGE(NSOLUTE,JBOUND,IBOUND,TSTOP,BCSTOP,TSTEP,TIME,
     &                 USTIME,USCONC,USBC)

      ELSE
*
*        both the b.c. and the flow variables have changed
*
         CALL QCHANGE(IMAX,Q,AREA,QLATIN,CLATIN,NSOLUTE,QVALUE,AVALUE,
     &                QWT,QINDEX,NFLOW,JBOUND,QSTOP,QSTEP,TSTOP,USCONC,
     &                QINVAL,CLVAL,DSDIST,USDIST,IBOUND,USBC,NBOUND)
         CALL PREPROC3(IMAX,DELTAX,Q,AREA,DFACE,QLATIN,HPLUSF,HPLUSB,
     &                 HMULTF,HMULTB,HDIV,HDIVF,HADV,GAMMA,BTERMS,AFACE,
     &                 A,C,AWORK,BWORK,TIMEB,NSOLUTE,LAM2DT,LHAT2DT,KD,
     &                 TWOPLUS,AREA2,ALPHA,LAMBDA,SGROUP,BTERMSN,BN,
     &                 TGROUP,GAMMAN,IGROUP,QLATINN,CLATINN,AREAN,
     &                 CLATIN)
         CALL BCCHANGE(NSOLUTE,JBOUND,IBOUND,TSTOP,BCSTOP,TSTEP,TIME,
     &                 USTIME,USCONC,USBC)

      ENDIF
*
*     reset TSTOP
*
      CALL SETSTOP(QSTOP,TSTOP,STOPTYPE,BCSTOP,IBOUND)

      RETURN
      END
