************************************************************************
*     
*     Subroutine REACT                Called by: DYNAMIC, DYNAMIC2
*
*     develop the finite difference equation for each segment.  The
*     resulting set of equations is solved to determine the main channel
*     concentrations at time 'N+1'.  The storage zone and streambed
*     sediment concentrations are then determined.
*
*     The solutes undergo first-order decay and/or sorption.
*
*     Local Variable
*     --------------
*     OLDCONC    concentration at time step N
* 
************************************************************************
      SUBROUTINE REACT(IMAX,JBOUND,DELTAX,Q,USCONC,USCONCN,AREA,DFACE,
     &                 DSBOUND,CONC,CONC2,NSOLUTE,GAMMA,AFACE,AWORK,
     &                 BWORK,LAM2DT,SGROUP2,SGROUP,LHATDT,SORB,LHAT2DT,
     &                 KD,QN,AREAN,GAMMAN,AFACEN,AN,CN,C,TWOPLUS,ALPHA,
     &                 BN,TGROUP,IGROUP)
      INCLUDE 'fmodules.inc'
*
*     subroutine arguments
*
      INTEGER*4 IMAX,JBOUND,NSOLUTE
      DOUBLE PRECISION DSBOUND
      DOUBLE PRECISION DELTAX(*),Q(*),AREA(*),DFACE(*),GAMMA(*),
     &                 AFACE(*),USCONCN(*),QN(*),AREAN(*),AFACEN(*),
     &                 GAMMAN(*),AN(*),CN(*),C(*),ALPHA(*)
      DOUBLE PRECISION USCONC(MAXBOUND,*),CONC(MAXSEG,*),
     &                 CONC2(MAXSEG,*),LAM2DT(MAXSEG,*),AWORK(MAXSEG,*),
     &                 BWORK(MAXSEG,*),SGROUP2(MAXSEG,*),
     &                 SGROUP(MAXSEG,*),LHATDT(MAXSEG,*),SORB(MAXSEG,*),
     &                 LHAT2DT(MAXSEG,*),KD(MAXSEG,*),TWOPLUS(MAXSEG,*),
     &                 BN(MAXSEG,*),TGROUP(MAXSEG,*),IGROUP(MAXSEG,*)
*
*     Local variables
*
      INTEGER*4 I,J
      DOUBLE PRECISION OLDCONC(MAXSEG)
*
*     arguments passed to SUBSTIT
*
      DOUBLE PRECISION D(MAXSEG)
*
*     begin loop for each solute
*
      DO 40 J = 1,NSOLUTE
*
*        save concentration from previous timestep
*        and compute D vector.
*
         OLDCONC(1) = CONC(1,J)
         D(1) = BN(1,J) * CONC(1,J)
     &        - CN(1) * CONC(2,J)
*storage
     &        + TGROUP(1,J)*CONC2(1,J)
     &        + ALPHA(1)*SGROUP2(1,J)/ TWOPLUS(1,J)
*inflow
     &        + IGROUP(1,J)
*sorp
     &        + 2.D0 * SORB(1,J) * SGROUP(1,J)
*adv term
     &        + .5d0*(QN(1)*USCONCN(J)/(AREAN(1)*DELTAX(1)) +
     &                Q(1)*USCONC(JBOUND,J)/(AREA(1)*DELTAX(1)))
*disp term
     &        + DFACE(1) * (AFACEN(1)*USCONCN(J)/AREAN(1)
     &           + AFACE(1)*USCONC(JBOUND,J)/AREA(1))
     &           / (DELTAX(1)*DELTAX(1))

         DO 20 I = 2,IMAX-1
            OLDCONC(I) = CONC(I,J)
            D(I) = - AN(I) * CONC(I-1,J)
     &           + BN(I,J) * CONC(I,J) 
     &           - CN(I) * CONC(I+1,J)
*storage
     &           + TGROUP(I,J)*CONC2(I,J)
     &           + ALPHA(I)*SGROUP2(I,J)/ TWOPLUS(I,J)
*inflow
     &           + IGROUP(I,J)
*sorp
     &           + 2.D0 * SORB(I,J) * SGROUP(I,J)
 20      CONTINUE
         I = IMAX  
         OLDCONC(I) = CONC(I,J)
         D(I) = -AN(I) * CONC(I-1,J) 
     &        + BN(I,J) * CONC(I,J)
*storage
     &        + TGROUP(I,J)*CONC2(I,J)
     &        + ALPHA(I)*SGROUP2(I,J)/ TWOPLUS(I,J)
*inflow
     &        + IGROUP(I,J)
*sorp
     &        + 2.D0 * SORB(I,J) * SGROUP(I,J)
*adv
     &         - 0.25D0*DSBOUND/DFACE(I)*(QN(2)/AREAN(I)+Q(I)/AREA(I))
*disp
     &        + 0.5D0*DSBOUND/DELTAX(I)
     &               *(AFACE(I)/AREA(I)+AFACEN(2)/AREAN(I))
*
*        Determine the main channel concentration by solving the system
*        of equations via the Thomas Algorithm (substitution step).
*
         CALL SUBSTIT(IMAX,CONC,AWORK,BWORK,C,D,J)
*
*        Compute the storage and sorbed concentrations
*
         DO 30 I = 1,IMAX
            CONC2(I,J) = ( (2.0D0-GAMMAN(I)-LAM2DT(I,J)-LHAT2DT(I,J))
     &                        *CONC2(I,J) 
     &                      + GAMMAN(I)*OLDCONC(I) + GAMMA(I)*CONC(I,J)
     &                      + 2.D0* SGROUP2(I,J)  )
     &                      / TWOPLUS(I,J) 
            SORB(I,J) = (  (2.D0-LHATDT(I,J))*SORB(I,J)
     &                   + KD(I,J)*LHATDT(I,J)
     &                   * (CONC(I,J)+OLDCONC(I))) / (2.D0+LHATDT(I,J))
 30      CONTINUE
*
*        save the upstream boundary condition for use next timestep
*
         USCONCN(J) = USCONC(JBOUND,J)

 40   CONTINUE

      RETURN
      END
