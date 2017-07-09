************************************************************************
*
*     Function YINTERP            Called by: BCCHANGE
*
*     Use linear interpolation to determine y(X), given X1, X2, y(X1)
*     and y(X2):
*
*     y(x) = y(X1) + {X - X1}{y(X2) - y(X1)}
*                    -----------------------
*                            X2 - X1
*
*     Variables
*     ---------
*     Y1       -   y(X1)
*     Y2       -   y(X2)
*     YINTERP  -   y
*
************************************************************************
      FUNCTION YINTERP(X,X1,X2,Y1,Y2)
*
*     function definition
*
      DOUBLE PRECISION YINTERP
*
*     subroutine arguments
*
      DOUBLE PRECISION X,X1,X2,Y1,Y2

      YINTERP = Y1 + (X - X1)/(X2-X1) * (Y2 - Y1)

      RETURN
      END
