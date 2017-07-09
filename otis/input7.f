************************************************************************
*
*     Subroutine INPUT7               Called by: INPUT
*
*     read the upstream boundary conditions
* SM: intializes two important variables: USTIME and USCONC(NBOUND,1) with time and
*                                           concentration of the upstream boundary
************************************************************************
      SUBROUTINE INPUT7(NBOUND,USTIME,USCONC,USCONCN,NSOLUTE,IBOUND,
     &                  USBC,QSTART,TFINAL)
*
*     dimensional parameters and logical devices
*
      INCLUDE 'fmodules.inc'
      INCLUDE 'lda.inc'
*
*     subroutine arguments
*
      INTEGER*4 IBOUND,NBOUND,NSOLUTE
      DOUBLE PRECISION QSTART,TFINAL
      DOUBLE PRECISION USTIME(*),USCONCN(*)
      DOUBLE PRECISION USCONC(MAXBOUND,*),USBC(MAXBOUND,*)
*
*     local variables
*
      INTEGER*4 I,J
      CHARACTER*500 BUFFER
*
*     Input Format Statements
*
 1000 FORMAT(2I5)
 1020 FORMAT(70D13.5)
*
*     Output Format Statements
*
 2000 FORMAT(//,10X,'Number of Upstream Boundary Conditions: ',I5,/,10X,
     &       'Boundary Condition Option             : ',I5,/)
 2020 FORMAT(/,10X,'Upstream Boundary Conditions, Solute #',I5,/,10X,
     &43('-'),/,10X,'Begin                Concentration',/,10X,
     &              'Time [hour]          [mass/L^3]',/,10X,34('-'))
 2040 FORMAT(/,10X,'Upstream Boundary Conditions, Solute #',I3,/,10X,
     &43('-'),/,10X,'Begin                Flux',/,10X,
     &              'Time [hour]          [mass/second]',/,10X,34('-'))
 2060 FORMAT(10X,1P,E11.4,9X,E12.4)
*
*     read the number of upstream boundary conditions and the boundary
*     condition option. Test to see if the number of boundary conditions
*     exceeds the maximum and for an invalid boundary condition option.
*
      CALL GETLINE(LDPARAM,BUFFER)
      READ (BUFFER,1000) NBOUND,IBOUND
      WRITE (LDECHO,2000) NBOUND,IBOUND
      IF (NBOUND .GT. MAXBOUND) CALL ERROR(5,NBOUND,MAXBOUND)
      IF (IBOUND .NE. 1 .AND. IBOUND .NE. 2 .AND. IBOUND .NE. 3)
     &   CALL ERROR3(3,IBOUND)
*
*     Read the upstream boundary conditions for each solute.
*
      DO 10 I = 1,NBOUND
         CALL GETLINE(LDPARAM,BUFFER)
         READ (BUFFER,1020) USTIME(I),(USBC(I,J),J=1,NSOLUTE)
 10   CONTINUE
*
*     for IBOUND=3, test to make sure that the time of the last b.c.
*     is at or after the simulation end time.  This is required for
*     determining a boundary concentration via interpolation.
*
      IF ((IBOUND .EQ. 3).AND.(USTIME(NBOUND) .LT. TFINAL))
     &   CALL ERROR2(6,USTIME(NBOUND),TFINAL)

*
*     write the boundary conditions for each solute
*
      DO 30 J = 1,NSOLUTE
         IF (IBOUND .EQ. 1 .OR. IBOUND .EQ. 3) THEN
            WRITE (LDECHO,2020) J
         ELSE
            WRITE (LDECHO,2040) J
         ENDIF
         DO 20 I = 1,NBOUND
            WRITE (LDECHO,2060) USTIME(I),USBC(I,J)
 20      CONTINUE
 30   CONTINUE
*
*     set the upstream boundary concentration (USCONC) based on the
*     specified boundary condition (USBC).  If the boundary condition
*     is specified as a flux (IBOUND=2) divide the boundary condition
*     by the flow at the upstream boundary to obtain the concentration.
*
      !Masud: IBOUND = 3 for continuous profile

      DO 50 J = 1,NSOLUTE
         DO 40 I = 1,NBOUND
            IF (IBOUND .EQ. 1 .OR. IBOUND .EQ. 3) THEN
               USCONC(I,J) = USBC(I,J)
            ELSEIF (IBOUND .EQ. 2) THEN
               USCONC(I,J) = USBC(I,J) / QSTART
            ENDIF
 40      CONTINUE
 50   CONTINUE
*
*     initialize the boundary conditions for time step 'N', USCONCN.
*
      DO 60 J = 1,NSOLUTE
         USCONCN(J) = USCONC(1,J)
 60   CONTINUE

      RETURN
      END
