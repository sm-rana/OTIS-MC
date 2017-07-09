************************************************************************
*
*     Subroutine PREPROC                Called by: MAININIT
*
*     preprocess parameter groups and initialize variables.  Most of the
*     preprocessing is completed in the following three routines:
*
*     PREPROC1 - develop parameter groups that are a function of the
*                segment length, and the interfacial dispersion coeff.
*
*     PREPROC2 - develop parameter groups that are functions of the
*                time-invariant model parameters.
*
*     PREPROC3 - develop parameter groups that are functions of the flow
*                variables and decompose the tridiagonal matrix.
*
************************************************************************
      SUBROUTINE PREPROC(IMAX,NBOUND,AREA2,DELTAX,Q,USTIME,AREA,ALPHA,
     &                   DISP,QLATIN,TSTART,TSTEP,QSTEP,NSOLUTE,LAMBDA,
     &                   LAMBDA2,CHEM,STOPTYPE,JBOUND,TIMEB,TSTOP,QSTOP,
     &                   DFACE,HPLUSF,HPLUSB,HMULTF,HMULTB,HDIV,HDIVF,
     &                   HADV,GAMMA,BTERMS,AFACE,A,C,AWORK,BWORK,LAM2DT,
     &                   BCSTOP,IBOUND,RHOLAM,KD,LAMHAT,LAMHAT2,LHAT2DT,
     &                   TWOPLUS,LHATDT,SGROUP,SGROUP2,CSBACK,BN,TGROUP,
     &                   IGROUP,CLATIN)
      INCLUDE 'fmodules.inc'
*
*     subroutine arguments
*
      CHARACTER*(*) CHEM,STOPTYPE
      INTEGER*4 IMAX,NBOUND,NSOLUTE,JBOUND,IBOUND
      DOUBLE PRECISION TSTART,TSTEP,QSTEP,TIMEB,TSTOP,QSTOP,BCSTOP
      DOUBLE PRECISION AREA2(*),DELTAX(*),Q(*),USTIME(*),AREA(*),
     &                 ALPHA(*),DISP(*),QLATIN(*),DFACE(*),HPLUSF(*),
     &                 HPLUSB(*),HMULTF(*),HMULTB(*),HDIV(*),HDIVF(*),
     &                 HADV(*),GAMMA(*),AFACE(*),A(*),C(*)
      DOUBLE PRECISION LAMBDA(MAXSEG,*),LAMBDA2(MAXSEG,*),
     &                 LAM2DT(MAXSEG,*),AWORK(MAXSEG,*),BWORK(MAXSEG,*),
     &                 RHOLAM(MAXSEG,*),KD(MAXSEG,*),LAMHAT(MAXSEG,*),
     &                 LAMHAT2(MAXSEG,*),LHAT2DT(MAXSEG,*),
     &                 TWOPLUS(MAXSEG,*),LHATDT(MAXSEG,*),
     &                 SGROUP(MAXSEG,*),SGROUP2(MAXSEG,*),
     &                 CSBACK(MAXSEG,*),BN(MAXSEG,*),TGROUP(MAXSEG,*),
     &                 BTERMS(MAXSEG,*),IGROUP(MAXSEG,*),
     &                 CLATIN(MAXSEG,*)
*
*     compute the inverse of the time step [/sec]
*

      IF (CHEM .NE. 'Steady-State') TIMEB = 1.D0 / (TSTEP * 3600.D0)

*
*     initialize the upstream boundary condition variables.  For a step
*     boundary condition (IBOUND = 1 or 2), the boundary condition will
*     change at user-specified times indicated by USTIME.  For a
*     continuous boundary condition (IBOUND =3), the b.c. will change
*     every time step.  In either case, set the time at which the
*     boundary condition changes, BCSTOP.
*
      USTIME(NBOUND+1) = 999.D99
      JBOUND = 1

      IF (IBOUND .EQ. 1 .OR. IBOUND .EQ. 2) THEN
         BCSTOP = USTIME(JBOUND+1)
      ELSEIF (IBOUND .EQ. 3) THEN
         BCSTOP = TSTART
      ENDIF
*
*     initialize the time at which a change in flow occurs.  if the flow
*     is steady, set QSTOP to an arbitrarily large number.
*
      IF (QSTEP .EQ. 0.D0) THEN
         QSTOP = 999.D99
      ELSE
         QSTOP = TSTART + QSTEP
      ENDIF
*
*     compute parameter groups used within CONSER/REACT/SSDIFF.  If the
*     run is 'Steady-State' call only PREPROC1.  Note that in the call
*     to PREPROC3, GAMMA, BTERMS, QLATIN, CLATIN and AREA are passed in
*     twice as GAMMAN, BTERMSN, QLATINN, CLATINN and AREAN are
*     undefined.
*

      CALL PREPROC1(DISP,DFACE,DELTAX,HPLUSF,HPLUSB,HMULTF,HMULTB,HDIV,
     &              HDIVF,HADV,IMAX)
      IF (CHEM .NE. 'Steady-State') THEN
         CALL PREPROC2(IMAX,NSOLUTE,TIMEB,LAMBDA2,LAM2DT,LAMHAT2,
     &                 LHAT2DT,LAMHAT,LHATDT,SGROUP,SGROUP2,RHOLAM,
     &                 CSBACK)
         CALL PREPROC3(IMAX,DELTAX,Q,AREA,DFACE,QLATIN,HPLUSF,HPLUSB,
     &                 HMULTF,HMULTB,HDIV,HDIVF,HADV,GAMMA,BTERMS,AFACE,
     &                 A,C,AWORK,BWORK,TIMEB,NSOLUTE,LAM2DT,LHAT2DT,KD,
     &                 TWOPLUS,AREA2,ALPHA,LAMBDA,SGROUP,BTERMS,BN,
     &                 TGROUP,GAMMA,IGROUP,QLATIN,CLATIN,AREA,CLATIN)
      ENDIF
*
*     determine TSTOP, the time at which the upstream boundary condition
*     or flow variables change.
*
      CALL SETSTOP(QSTOP,TSTOP,STOPTYPE,BCSTOP,IBOUND)

      RETURN
      END
