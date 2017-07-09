###########################################################################
# 
#  Makefile for OTIS - MOD40
# 
#  This makefile is accessed by ../../makefile to compile OTIS.
#
#############################################################################
#
#  Notes:
#
#
#  I. Fortran Definitions
#
#      1) Debug Flags
#
#         FFLAGS = -C -g -u
#
#      2) Sample Optimization Flags
#
#         OS/Hardware       Compiler                    FFLAGS
#         -------------------------------------------------------
#         SunOS/SUN         SPARCompiler Fortran        -fast -O3
#         Solaris/SUN       SPARCompiler Fortran        -fast -O3
#         Linux/Pentium II  GNU F77                     -O
#	  DG-UX/DG          Green Hills Fortran-88000   -OLM -u
#	  AIX/IBM           AIX XL FORTRAN Comp/6000    -O3
#	  Digital UNIX/DEC  Digital Fortran 77          -O4	 
#
#
#  II. OTIS Object files
#
#        OBJ  - source code that is local to OTIS
#        SHARE - source code that is shared with OTIS-P (../share)
#
###########################################################################

F77 = f77 $(FFLAGS)
FFLAGS = -O -u

#FC = gfortran
#F77 = gfortran $(FFLAGS)
#FFLAGS = -O -u -fno-underscoring

SHARE = ../share/share.a
OBJ = closef.o input.o ldainit.o main.o maininit.o mainrun.o openfile.o
###########################################################################
#
#   OTIS
#
###########################################################################
otis:	$(OBJ)
	$(F77) -o $@ $(OBJ) $(SHARE)
#	$(F77) $(FFLAGS) -o $@ $(OBJ) $(SHARE)    
#	$(F77) $(FFLAGS) -o $@ $(OBJ) $(SHARE)
#	mv -f $@ ../../$@

