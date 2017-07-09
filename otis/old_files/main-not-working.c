#include <stdio.h>
#include <time.h>
// #include "iniparser.h"

// MASUD:
// These values are defined here on the basis of fmodumes.inc
// MAKE SURE the values match in the file and here!!
#define MAXREACH    30
#define MAXPRINT    30
#define MAXBOUND    20000
#define MAXSEG      5000
#define MAXSOLUTE   3
#define MAXFLOWLOC  30

//int mainold_(int*,double*,double*,double*,double*,int*,double*);
/*
int mainold_(int* mcmc_iter, double *new_area, double *new_disp, double* new_area2,
             double *new_alpha, int* nseg_tmp, double* rchlen_tmp, double *time_start
             double *sim_time, double *conc_main, double *conc_ts);
*/

int mainold_(int* mcmc_len, double* new_area, double* new_disp, double* new_area2,
             double* new_alpha, int* nreach, int* nseg_tmp, double* rchlen_tmp, double*time_start);
/*
int mainold_( int* IMAX, int* NPRINT, int* IPRINT, int* NSOLUTE, int* PRTOPT, int* NFLOW,
		 int* JBOUND, int* IBOUND, int* NBOUND, int* ISORB, int* PINDEX, int* QINDEX,
		 double* DSBOUND,   double* TSTART,     double* TFINAL, double* TSTEP, double* QSTEP,
		 double* TIMEB,     double* TSTOP,      double* QSTOP, double* BCSTOP,
    	 double* DELTAX,    double* Q,          double* USTIME, double* AREA,
		 double* QLATIN,   double* WT,        double* QVALUE, double* AVALUE,
         double* QWT,      double* DSDIST,    double* USDIST, double* DFACE, double* HPLUSF,
		 double* HPLUSB,   double* HMULTF,    double* HMULTB, double* HDIV,
         double* HDIVF,    double* HADV,      double* GAMMA, double* AFACE, double* A, double* C,
         double* QINVAL,   double* USCONCN,   double* QN, double* AREAN, double* QLATINN,
		 double* AFACEN,   double* GAMMAN,    double* AN, double* CN,
         double* AREA2,    double* ALPHA,
    	 double* USCONC[],  double* CLATIN[],   double* CONC[], double* CONC2[],
         double* CLVAL[],   double* LAMBDA[],   double* LAM2DT[], double* AWORK[],
		 double* BWORK[],   double* USBC[],     double* SGROUP2[], double* SGROUP[],
		 double* LHATDT[],  double* SORB[],     double* LHAT2DT[], double* KD[],
		 double* CLATINN[], double* TWOPLUS[],  double* BN[], double* TGROUP[],
		 double* BTERMS[],  double* BTERMSN[],  double* IGROUP[],
		 char* CHEM, char* STOPTYPE,
		 int* mcmc_iter, double* new_area, double* new_disp, double* new_area2,
         double* new_alpha, int* nreach, int* nseg_tmp, double* rchlen_tmp, double*time_start);
*/
int main(int argc, char** argv )
{
    printf("Beginning running main()\n");

    int mcmc_iter, mcmc_len;
    int nreach, nseg_tmp;
    double rchlen_tmp, time_start;
    double new_area, new_disp, new_area2, new_alpha;

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
    /*
    int lastseg[MAXREACH];
    int nseg[MAXREACH];
    double rchlen[MAXREACH];
    //SM: NSEG and RCHLEN gets initialized in input2.f need to trickle it up all the way to main()
    //double sim_time[MAXBOUND], conc_main[MAXBOUND],conc_ts[MAXBOUND];
    */

    int IMAX,NPRINT,IPRINT,NSOLUTE,PRTOPT,NFLOW,JBOUND,IBOUND,NBOUND,ISORB;
    int PINDEX[MAXPRINT], QINDEX[MAXSEG];

    double DSBOUND, TSTART, TFINAL,TSTEP, QSTEP, TIMEB, TSTOP, QSTOP, BCSTOP;

    double  DELTAX[MAXSEG],Q[MAXSEG],USTIME[MAXBOUND+1],
            AREA[MAXSEG],QLATIN[MAXSEG],WT[MAXPRINT],
            QVALUE[MAXFLOWLOC],AVALUE[MAXFLOWLOC],
            QWT[MAXSEG],DSDIST[MAXSEG],USDIST[MAXSEG],
            DFACE[MAXSEG],HPLUSF[MAXSEG],HPLUSB[MAXSEG],
            HMULTF[MAXSEG],HMULTB[MAXSEG],HDIV[MAXSEG],
            HDIVF[MAXSEG],HADV[MAXSEG],GAMMA[MAXSEG],
            AFACE[MAXSEG],A[MAXSEG],C[MAXSEG],
            QINVAL[MAXFLOWLOC+1],USCONCN[MAXSOLUTE],QN[2],
            AREAN[MAXSEG],QLATINN[MAXSEG],AFACEN[2],
            GAMMAN[MAXSEG],AN[MAXSEG],CN[MAXSEG],
            AREA2[MAXSEG],ALPHA[MAXSEG];

    double  USCONC[MAXSOLUTE][MAXBOUND],
            CLATIN[MAXSOLUTE][MAXSEG],CONC[MAXSOLUTE][MAXSEG],
            CONC2[MAXSOLUTE][MAXSEG],
            CLVAL[MAXSOLUTE][MAXFLOWLOC+1],
            LAMBDA[MAXSOLUTE][MAXSEG],
            LAM2DT[MAXSOLUTE][MAXSEG],
            AWORK[MAXSOLUTE][MAXSEG],BWORK[MAXSOLUTE][MAXSEG],
            USBC[MAXSOLUTE][MAXBOUND],
            SGROUP2[MAXSOLUTE][MAXSEG],
            SGROUP[MAXSOLUTE][MAXSEG],
            LHATDT[MAXSOLUTE][MAXSEG],SORB[MAXSOLUTE][MAXSEG],
            LHAT2DT[MAXSOLUTE][MAXSEG],KD[MAXSOLUTE][MAXSEG],
            CLATINN[MAXSOLUTE][MAXSEG],
            TWOPLUS[MAXSOLUTE][MAXSEG],BN[MAXSOLUTE][MAXSEG],
            TGROUP[MAXSOLUTE][MAXSEG],
            BTERMS[MAXSOLUTE][MAXSEG],
            BTERMSN[MAXSOLUTE][MAXSEG],
            IGROUP[MAXSOLUTE][MAXSEG];

    char CHEM[12], STOPTYPE[12];


////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
    printf("End of variable declaration\n");

    mcmc_len = 2;
    mcmc_iter = 0;
    new_area = 6.00E-2;
    new_disp = 3.20E-1;
    new_area2 = 2.60E-2;
    new_alpha = 2.90E-4;


    printf("MCMC iter= %d\n----------------------------------------\n", mcmc_iter);
    /*
    mainold_(&IMAX, &NPRINT, &IPRINT, &NSOLUTE, &PRTOPT, &NFLOW,
		 &JBOUND, &IBOUND, &NBOUND, &ISORB, PINDEX, QINDEX,
		 &DSBOUND,   &TSTART,     &TFINAL, &TSTEP, &QSTEP,
		 &TIMEB,     &TSTOP,      &QSTOP, &BCSTOP,
    	 DELTAX,   Q,         USTIME, AREA,
		 QLATIN,   WT,        QVALUE, AVALUE,
         QWT,      DSDIST,    USDIST, DFACE, HPLUSF,
		 HPLUSB,   HMULTF,    HMULTB, HDIV,
         HDIVF,    HADV,      GAMMA, AFACE, A, C,
         QINVAL,   USCONCN,   QN, AREAN, QLATINN,
		 AFACEN,  GAMMAN,    AN, CN,
         AREA2,    ALPHA,
    	 USCONC,  CLATIN,   CONC, CONC2,
         CLVAL,   LAMBDA,   LAM2DT, AWORK,
		 BWORK,   USBC,     SGROUP2, SGROUP,
		 LHATDT,  SORB,     LHAT2DT, KD,
		 CLATINN, TWOPLUS,  BN, TGROUP,
		 BTERMS,  BTERMSN,  IGROUP,
		 CHEM, STOPTYPE,
		 &mcmc_iter, &new_area, &new_disp, &new_area2,
         &new_alpha, &nreach, &nseg_tmp, &rchlen_tmp, &time_start);
    */

    mainold_(&mcmc_iter, &new_area, &new_disp, &new_area2,
             &new_alpha, &nreach, &nseg_tmp, &rchlen_tmp, &time_start);

    printf("USCONC = %f, %f\n",USCONC[0][0], USCONC[1][0]);
    printf("USCONC = %f, %f\n",USCONC[0][1], USCONC[1][1]);
    printf("USCONC = %f, %f\n",USCONC[0][2], USCONC[1][2]);
    printf("USCONC = %f, %f\n",USCONC[0][3], USCONC[1][3]);
    printf("USCONC = %f, %f\n",USCONC[0][4], USCONC[1][4]);
    printf("USCONC = %f, %f\n",USCONC[0][5], USCONC[1][5]);
    printf("USCONC = %f, %f\n",USCONC[0][6], USCONC[1][6]);
    printf("main(): nreach = %d\n",nreach);
    printf("------------------ X -------------------\n\n");

    //MAINOLD(mcmc_iter, new_area, new_disp, new_area2, new_alpha, CONC,CONC2)
    /*
        Plan for the MCMC run
        maininit_mod()
            input_mod()
                input2()    <-- perturb parameters
            preproc()
        mainrun_mod()   <-- get CONC1 and CONC2 from it
                        <-- interporlate at the observation time
                        <-- compute likelihood
        Trickle back:
            from input2.f --> input.f --> :
                NREACH
                NSEG
                RCHLEN
    */
    // maininit_mod_()
    // printf("In (myc.c) The value of x obtained from fortran = %d",x);

    printf("main() --> nseg_tmp = %d, rchlen_tmp = %.1f, time_start = %.1f \n\n",
           nseg_tmp, rchlen_tmp, time_start);

    return 0;
}
