************************************************************************
*
*     Subroutine SSDIFF               Called by: STEADY
*
*     Develop the finite difference equation for each segment.  The
*     resulting equation set is solved to determine the steady-state
*     concentrations.
*
************************************************************************
      SUBROUTINE SSDIFF(IMAX,DELTAX,Q,USCONC,ALPHA,AREA,AREA2,DFACE,
     &                  DSBOUND,QLATIN,CLATIN,CONC,HPLUSF,HPLUSB,HMULTF,
     &                  HMULTB,HDIV,HDIVF,HADV,NSOLUTE,LAMBDA,LAMBDA2,
     &                  LAMHAT2,CSBACK)
      INCLUDE 'fmodules.inc'
*
*     subroutine arguments
*
      INTEGER*4 IMAX,NSOLUTE
      DOUBLE PRECISION DSBOUND,DELTAX(*),Q(*),ALPHA(*),AREA(*),AREA2(*),
     &                 DFACE(*),QLATIN(*),HPLUSF(*),HPLUSB(*),HMULTF(*),
     &                 HMULTB(*),HDIV(*),HDIVF(*),HADV(*),
     &                 USCONC(MAXBOUND,*),CLATIN(MAXSEG,*),
     &                 CONC(MAXSEG,*),LAMBDA(MAXSEG,*),
     &                 LAMBDA2(MAXSEG,*),LAMHAT2(MAXSEG,*),
     &                 CSBACK(MAXSEG,*)
*
*     Local variables
*
      INTEGER*4 I,J
      DOUBLE PRECISION AFACE(MAXSEG),AREADISP(0:MAXSEG)
*
*     arguments passed to MATRIX
*
      DOUBLE PRECISION A(MAXSEG),B(MAXSEG),C(MAXSEG),D(MAXSEG)
*
*     begin loop for each solute
*
      DO 20 J = 1,NSOLUTE
*
*        Compute the cross sectional area at the interface between
*        segments one and two.
*
         AFACE(1) = HDIV(1) * AREA(2) + HDIVF(1) * AREA(1)
*
*        set coefficients for segment 1.  The coefficients are developed
*        on a term by term basis, so that changes to an individual term
*        of the differential equation can be easily made.
*
         A(1) = 0.D0
*
*        contributions from the advection term,
*
*                    Q   dC
*                -  --- ----
*                    A   dx
*
*
         B(1) = Q(1) * HDIVF(1) / (AREA(1) * DELTAX(1))
         C(1) = Q(1) / (AREA(1) * HPLUSF(1))
         D(1) = Q(1) * USCONC(1,J) / (AREA(1) * DELTAX(1))
*
*        contributions from the dispersion term,
*
*                    1   d        dc
*                +  --- ---- (AD ----)
*                    A   dx       dx
*
*
         AREADISP(1) = AFACE(1) * DFACE(1)
         AREADISP(0) = AREADISP(1)

         B(1) = B(1) + 2.D0 *
     &          (AREADISP(1)/HPLUSF(1) + AREADISP(0)/DELTAX(1))
     &          / (AREA(1)*DELTAX(1))
         C(1) = C(1) - 2.D0 * AREADISP(1) / (AREA(1) * HMULTF(1))
         D(1) = D(1) + 2.D0 * AREADISP(0) * USCONC(1,J) /
     &                 (AREA(1) * DELTAX(1) * DELTAX(1))
*
*        contributions from the lateral inflow term,
*
*                    q
*                +  --- (Cl - C)
*                    A
*
*
         B(1) = B(1) + QLATIN(1) / AREA(1)
         D(1) = D(1) + QLATIN(1) * CLATIN(1,J) / AREA(1)
*
*        contributions from the transient storage term,
*
*                 +  alpha (Cs - C)
*
*        note that if alpha equals zero, there is no contribution
*
         IF (ALPHA(1) .NE. 0.D0) THEN
            B(1) = B(1)
     &           + ALPHA(1)*AREA2(1)*(LAMBDA2(1,J)+LAMHAT2(1,J)) /
     &        (ALPHA(1)*AREA(1)+AREA2(1)*(LAMBDA2(1,J)+LAMHAT2(1,J)))
            D(1) = D(1) + ALPHA(1)*AREA2(1)*LAMHAT2(1,J)*CSBACK(1,J)
     &           / (ALPHA(1)*AREA(1) +
     &                       AREA2(1)*(LAMBDA2(1,J)+LAMHAT2(1,J)))
         ENDIF
*
*        contributions from the decay term,
*
*                 -  lambda C
*
         B(1) = B(1) + LAMBDA(1,J)
*
*        set coefficients for segments 2 through IMAX-1.  The coefficients
*        are developed on a term by term basis, so that changes to an
*        individual term of the differential equation can be easily made.
*
         DO 10 I = 2,IMAX-1
*
*           compute the cross-sectional area at the interface
*           between segments I and I+1.
*
            AFACE(I) = HDIV(I) * AREA(I+1) + HDIVF(I) * AREA(I)
*
*           contributions from the advection term,
*
*                       Q   dC
*                   -  --- ----
*                       A   dx
*
*
            A(I) = - Q(I) / (AREA(I) * HPLUSB(I))
            B(I) = Q(I) * 2.D0 * HADV(I) / AREA(I)
            C(I) = Q(I) / (AREA(I) * HPLUSF(I))
            D(I) = 0.D0
*
*           contributions from the dispersion term,
*
*                       1   d        dc
*                   +  --- ---- (AD ----)
*                       A   dx       dx
*
*
            AREADISP(I) = AFACE(I) * DFACE(I)
            A(I) = A(I) - 2.D0 * AREADISP(I-1) / (AREA(I) * HMULTB(I))
            B(I) = B(I) + 2.D0 *
     &             (AREADISP(I)/HPLUSF(I) + AREADISP(I-1)/HPLUSB(I))
     &             / (AREA(I)*DELTAX(I))
            C(I) = C(I) - 2.D0 * AREADISP(I) / (AREA(I) * HMULTF(I))
*
*           contributions from the lateral inflow term,
*
*                       q
*                   +  --- (Cl - C)
*                       A
*
*
            B(I) = B(I) + QLATIN(I) / AREA(I)
            D(I) = D(I) + QLATIN(I) * CLATIN(I,J) / AREA(I)
*
*           contributions from the transient storage term,
*
*                    +  alpha (Cs - C)
*
*           note that if alpha equals zero, there is no contribution
*
            IF (ALPHA(I) .NE. 0.D0) THEN
               B(I) = B(I)
     &              + ALPHA(I)*AREA2(I)*(LAMBDA2(I,J)+LAMHAT2(I,J)) /
     &              (ALPHA(I)*AREA(I)
     &                      + AREA2(I)*(LAMBDA2(I,J)+LAMHAT2(I,J)))
               D(I) = D(I)
     &              + ALPHA(I)*AREA2(I)*LAMHAT2(I,J)*CSBACK(I,J)
     &              / (ALPHA(I)*AREA(I) +
     &                          AREA2(I)*(LAMBDA2(I,J)+LAMHAT2(I,J)))
            ENDIF
*
*           contributions from the decay term,
*
*                    -  lambda C
*
            B(I) = B(I) + LAMBDA(I,J)

 10      CONTINUE
*
*        Set coefficients for the last segment
*
         I = IMAX
         C(I) = 0.D0
*
*        compute the cross-sectional area at the last interface
*
         AFACE(I) =   AREA(I)
*
*        contributions from the advection term,
*
*                    Q   dC
*                -  --- ----
*                    A   dx
*
*
         A(I) = - Q(I) / (AREA(I) * HPLUSB(I))
         B(I) =   Q(I) * (1.D0 - DELTAX(I-1)/HPLUSB(I))
     &            / (AREA(I)*DELTAX(I))
         D(I) = - Q(I) * DSBOUND / (2.D0 * AREA(I) * DFACE(I))
*
*        contributions from the dispersion term,
*
*                    1   d        dc
*                +  --- ---- (AD ----)
*                    A   dx       dx
*
*
         A(I) = A(I) - 2.D0 * AREADISP(I-1) / (AREA(I) * HMULTB(I))
         B(I) = B(I) + 2.D0 * AREADISP(I-1) / (AREA(I) * HMULTB(I))
         D(I) = D(I) + AFACE(I) * DSBOUND / (AREA(I) * DELTAX(I))
*
*        contributions from the lateral inflow term,
*
*                    q
*                +  --- (Cl - C)
*                    A
*
*
         B(I) = B(I) + QLATIN(I) / AREA(I)
         D(I) = D(I) + QLATIN(I) * CLATIN(I,J) / AREA(I)
*
*        contributions from the transient storage term,
*
*                 +  alpha (Cs - C)
*
*        note that if alpha equals zero, there is no contribution
*
         IF (ALPHA(I) .NE. 0.D0) THEN
            B(I) = B(I)
     &           + ALPHA(I)*AREA2(I)*(LAMBDA2(I,J)+LAMHAT2(I,J)) /
     &           (ALPHA(I)*AREA(I)
     &                   + AREA2(I)*(LAMBDA2(I,J)+LAMHAT2(I,J)))
            D(I) = D(I)
     &           + ALPHA(I)*AREA2(I)*LAMHAT2(I,J)*CSBACK(I,J)
     &           / (ALPHA(I)*AREA(I) +
     &                       AREA2(I)*(LAMBDA2(I,J)+LAMHAT2(I,J)))
         ENDIF
*
*        contributions from the decay term,
*
*                 -  lambda C
*
         B(I) = B(I) + LAMBDA(I,J)
*
*        Determine the in-stream concentration by solving the system of
*        equations via the Thomas Algorithm.
*
         CALL MATRIX(IMAX,CONC,A,B,C,D,J)

 20   CONTINUE

      RETURN
      END
