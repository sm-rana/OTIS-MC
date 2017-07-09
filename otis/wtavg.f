************************************************************************
*
*     Function WTAVG            Called by: QAINIT
*
*     compute the weighted average value of QLATIN, QLATOUT, and CLATIN,
*     when a flow location is within a given segment. (Unsteady Flow
*     Conditions only)
*
************************************************************************
      FUNCTION WTAVG(QDIST,USVAL,DSVAL)
*
*     function arguments
*
      DOUBLE PRECISION WTAVG,QDIST,USVAL,DSVAL

      WTAVG = (0.5D0 - QDIST) * USVAL + (0.5D0 + QDIST) * DSVAL

      RETURN
      END
