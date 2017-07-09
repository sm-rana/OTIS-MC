************************************************************************
*
*     Subroutine change_area
*     It replicates the behavior of QSTEADY to change X-Sectional Area for
*     Monte Carlo Simulation
*
*     Called by: INPUTQ
*
*
************************************************************************
      SUBROUTINE get_fmodules_inc(MAXREACH1, MAXPRINT1, MAXBOUND1,
     &   MAXSEG1, MAXSOLUTE1, MAXFLOWLOC1)

*     MASUD:
*     Importing the constant values to the main() function in the C interface
*     from the file 'fmodules.inc'
*     Feb 20, 2017
*
*
      INCLUDE 'fmodules.inc'
*
*     subroutine arguments
*

      integer*4 MAXREACH1, MAXPRINT1, MAXBOUND1, MAXSEG1, MAXSOLUTE1,
     & MAXFLOWLOC1

      MAXREACH1 = MAXREACH
      MAXPRINT1 = MAXPRINT
      MAXBOUND1 = MAXBOUND
      MAXSEG1 = MAXSEG
      MAXSOLUTE1 = MAXSOLUTE
      MAXFLOWLOC1 = MAXFLOWLOC


      RETURN
      END
