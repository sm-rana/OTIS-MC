************************************************************************
*
*     Subroutine CONSER                Called by: DYNAMIC, DYNAMIC2
*
*     develop the finite difference equation for each segment.  The
*     resulting set of equations is solved to determine the main channel
*     concentration at time 'N+1'.  The storage zone concentration is
*     then determined.
*
*     This routine is for the case of conservative transport.
*
*     Local Variable:
*     ---------------
*     OLDCONC - main channel concentration at time step N
*
************************************************************************
      SUBROUTINE CONSER(IMAX,JBOUND,DELTAX,Q,USCONC,USCONCN,AREA,DFACE,
     &                  DSBOUND,CONC,CONC2,NSOLUTE,GAMMA,AFACE,AWORK,
     &                  BWORK,QN,AREAN,GAMMAN,AFACEN,AN,CN,C,TWOPLUS,BN,
     &                  TGROUP,IGROUP)
      INCLUDE 'fmodules.inc'
*
*     subroutine arguments
*
      INTEGER*4 IMAX,JBOUND,NSOLUTE
      DOUBLE PRECISION DSBOUND
      DOUBLE PRECISION DELTAX(*),Q(*),AREA(*),DFACE(*),GAMMA(*),
     &                 AFACE(*),USCONCN(*),QN(*),AREAN(*),AFACEN(*),
     &                 GAMMAN(*),AN(*),CN(*),C(*)
      DOUBLE PRECISION USCONC(MAXBOUND,*),CONC(MAXSEG,*),
     &                 CONC2(MAXSEG,*),AWORK(MAXSEG,*),BWORK(MAXSEG,*),
     &                 TWOPLUS(MAXSEG,*),BN(MAXSEG,*),TGROUP(MAXSEG,*),
     &                 IGROUP(MAXSEG,*)
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
         D(1) = BN(1,1)*CONC(1,J) - CN(1)*CONC(2,J)
*storage
     &          + TGROUP(1,1)*CONC2(1,J)
*inflow
     &          + IGROUP(1,J)
*adv term
     &          + .5d0*(QN(1)*USCONCN(J)/(AREAN(1)*DELTAX(1)) +
     &             Q(1)*USCONC(JBOUND,J)/(AREA(1)*DELTAX(1)))
*disp term
     &          + DFACE(1)*(AFACEN(1)*USCONCN(J)/AREAN(1) +
     &                      AFACE(1)*USCONC(JBOUND,J)/AREA(1))
     &            / (DELTAX(1)*DELTAX(1))
         DO 20 I = 2,IMAX-1
            OLDCONC(I) = CONC(I,J)
            D(I) = - AN(I)*CONC(I-1,J) + BN(I,1)*CONC(I,J)
     &             - CN(I)*CONC(I+1,J)
*storage
     &             + TGROUP(I,1)*CONC2(I,J)
*inflow
     &             + IGROUP(I,J)
 20      CONTINUE
         I = IMAX
         OLDCONC(I) = CONC(I,J)
         D(I) = - AN(I)*CONC(I-1,J) + BN(I,1)*CONC(I,J)
*storage
     &          + TGROUP(I,1)*CONC2(I,J)
*inflow
     &          + IGROUP(I,J)
*adv
     &          - 0.25D0*DSBOUND/DFACE(I)*(QN(2)/AREAN(I)+Q(I)/AREA(I))
*disp
     &          + 0.5D0*DSBOUND/DELTAX(I)
     &                 *(AFACE(I)/AREA(I)+AFACEN(2)/AREAN(I))
*
*        Determine the main channel concentration by solving the system
*        of equations via the Thomas Algorithm (substitution step).
*
        !write(*,*) "SM: conser.f D[] --->"
        !write(*,*) C(1),C(2),C(3),C(4)
        !write(*,*) C(17),C(18),C(19),C(20)
        !write(*,*) D(1),D(2),D(3),D(4)
        !write(*,*) D(17),D(18),D(19),D(20)
        !write(*,*) AN(1),AN(2),AN(3),AN(4),AN(5),AN(6),AN(7)
        !write(*,*)BN(1,1),BN(2,1),BN(3,1),BN(4,1),BN(5,1),BN(6,1),BN(7,1)
        !write(*,*)BWORK(1,1),BWORK(2,1),BWORK(3,1),BWORK(4,1),BWORK(5,1)
        !write(*,*) "END of conser.f --->"

         CALL SUBSTIT(IMAX,CONC,AWORK,BWORK,C,D,J)
        !write(*,*) "SM: conser.f D[] --->"
        !write(*,*) C(1),C(2),C(3),C(4),C(5)
        !write(*,*) D(1),D(2),D(3),D(4)
        !write(*,*) D(17),D(18),D(19),D(20)
        !write(*,*) AN(1),AN(2),AN(3),AN(4),AN(5),AN(6),AN(7)
        !write(*,*)BN(1,1),BN(2,1),BN(3,1),BN(4,1),BN(5,1),BN(6,1),BN(7,1)
        !write(*,*)BWORK(1,1),BWORK(2,1),BWORK(3,1),BWORK(4,1),BWORK(5,1)
        !write(*,*) "END of conser.f --->"
*
*        Compute the storage zone concentration
*
         DO 30 I = 1,IMAX
            CONC2(I,J) = ( (2.D0-GAMMAN(I)) * CONC2(I,J)
     &                   + GAMMAN(I)*OLDCONC(I) + GAMMA(I)*CONC(I,J))
     &                      / TWOPLUS(I,1)
 30      CONTINUE
*
*        save the upstream boundary condition for use next timestep
*
         USCONCN(J) = USCONC(JBOUND,J)

 40   CONTINUE

      RETURN
      END
