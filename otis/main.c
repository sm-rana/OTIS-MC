#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h> //to get DBL_MAX

#include "rvgs.h"
#include "rngs.h"
#include "iniparser.h"

// MASUD:
// These values are defined here on the basis of fmodumes.inc
// MAKE SURE the values match in the file and here!!
#define MAXREACH    1
#define MAXPRINT    1
#define MAXBOUND    10000
#define MAXSEG      5000
#define MAXSOLUTE   1
#define MAXFLOWLOC  1

/* Original fmodules.inc
      PARAMETER (MAXREACH = 30)
      PARAMETER (MAXPRINT = 30)
      PARAMETER (MAXBOUND = 20000)
      PARAMETER (MAXSEG = 5000)
      PARAMETER (MAXSOLUTE = 3)
      PARAMETER (MAXFLOWLOC = 30)
*/
//int mainold_(int*,double*,double*,double*,double*,int*,double*);
/*
int mainold_(int* mcmc_iter, double *new_area, double *new_disp, double* new_area2,
             double *new_alpha, int* nseg_tmp, double* rchlen_tmp, double *time_start
             double *sim_time, double *conc_main, double *conc_ts);
*/

double likelihood(double *obs, double *mod, int nobs, double sigma, double *ll);

int mainold_(int* IMAX, int* NPRINT, int* IPRINT, int* NSOLUTE, int* PRTOPT, int* NFLOW, int* JBOUND, int* IBOUND,
             int* NBOUND, int* ISORB, int* PINDEX, int* QINDEX,
             double* DSBOUND, double* TSTART, double* TFINAL, double* TSTEP,
             double* QSTEP, double* TIMEB, double* TSTOP, double* QSTOP, double* BCSTOP,
             double* DELTAX, double* Q, double* USTIME, double* AREA, double* QLATIN, double* WT,
             double* QVALUE, double* AVALUE, double* QWT, double* DSDIST, double* USDIST,
             double* DFACE, double* HPLUSF, double* HPLUSB, double* HMULTF, double* HMULTB, double* HDIV,
             double* HDIVF, double* HADV, double* GAMMA, double* AFACE, double* A, double* C,
             double* QINVAL, double* USCONCN, double* QN, double* AREAN, double* QLATINN, double* AFACEN,
             double* GAMMAN, double* AN, double* CN, double* AREA2, double* ALPHA,
             double* USCONC, double* CLATIN, double* CONC, double* CONC2, double* CLVAL, double* LAMBDA, double* LAM2DT, double* AWORK,
             double* BWORK, double* USBC, double* SGROUP2, double* SGROUP, double* LHATDT, double* SORB, double* LHAT2DT, double* KD, double* CLATINN,
             double* TWOPLUS, double* BN, double* TGROUP, double* BTERMS, double* BTERMSN, double* IGROUP,
             char *CHEM, char *STOPTYPE,
             int* IDECAY, double* XSTART, double* QSTART, // <-- local variables maininit.f
             double* PRTLOC, double*  DISP,  double* FLOWLOC, double* QLATOUT, double* DIST, // <-- local variables maininit.f
             double* LAMBDA2, double* LAMHAT, double* LAMHAT2, double* RHOLAM, double* CSBACK, // <-- local variables maininit.f
             int* ncout, int* mcmc_len, double* new_area, double* new_disp, double* new_area2,
             double* new_alpha, int* nreach, int* nseg_tmp, double* rchlen_tmp, double*time_start,
             double* sim_time, double*conc_main, double* conc_ts, double* test);
             //double* USCONC[], double* CLATIN[], double* CONC[], double* CONC2[], double* CLVAL[], double* LAMBDA[], double* LAM2DT[],
             //double* AWORK[], double* BWORK[], double* USBC[], double* SGROUP2[], double* SGROUP[], double* LHATDT[], double* SORB[],
             //double* LHAT2DT[], double* KD[], double* CLATINN[], double* TWOPLUS[], double* BN[], double* TGROUP[], double* BTERMS[],
             //double* BTERMSN[], double* IGROUP[],

//int get_fmodules_inc_(int* MAXREACH, int* MAXPRINT, int* MAXBOUND, int* MAXSEG, int* MAXSOLUTE, int* MAXFLOWLOC);

int main(int argc, char** argv ){

    int i, j, read_flag, count, ind_tmp, ind_cur, count_start, count_final;
    int mcmc_iter, mcmc_len, burn_in, ncout, nobs, nreach, nseg_tmp, accept, reject, flag_run;
    int *tind, nline;

    double rchlen_tmp, time_start, printStep;
    double new_area, new_disp, new_area2, new_alpha, time_tmp, conc1_tmp, conc2_tmp, param_tmp;
    double *sim_time, *obs_time, *conc_main, *obs_conc1, *conc_ts, *conc_interp, *tfrac;
    double *chain_A, *chain_As, *chain_D, *chain_Al, *chain_LL;

    double tmp_ll, cur_ll, last_ll, pstd_A, pstd_D, pstd_As, pstd_Al, std_obs, time_mcmc;  // pstd --> proposal standard deviation

    long seed;

    double *test;   // This variable is only for debug only.
    clock_t time_beginning, time_end;

//    int tmp1=5;
//    double tmpmod[]={10,12,14,15,16};
//    double tmpobs[]={10.5,12.9,13.5,11.1,14.9};

    char *fname, *fname2;
    FILE *fp_obs, *fp_out, *fp_debug;
    dictionary *dict;

    sim_time = (double*)calloc(MAXBOUND,sizeof(double));
    conc_main = (double*)calloc(MAXBOUND,sizeof(double));
    conc_ts = (double*)calloc(MAXBOUND,sizeof(double));

    test = (double*)calloc(6,sizeof(double));

    //fname = "observation01.inp";


////////////////////////////////////////////////////////////////////////////////////////
///////////// Following contains propagated variables from Fortran MAIN ////////////////
////////////////////////////////////////////////////////////////////////////////////////
/*
    IMAX    equals to NSEG in params file
*/
    int IMAX,NPRINT,IPRINT,NSOLUTE,PRTOPT,NFLOW,JBOUND,IBOUND,NBOUND,ISORB;
    int PINDEX[MAXPRINT], QINDEX[MAXSEG];
    int debug_flag[10];

    double DSBOUND,TSTART,TFINAL,TSTEP,QSTEP,TIMEB,TSTOP,QSTOP,BCSTOP;
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
    /*

    double *DELTAX,*Q,*USTIME,*AREA,*QLATIN,*WT,*QVALUE,*AVALUE,*QWT,*DSDIST,*USDIST,*DFACE,

    *DELTAX,*Q,*USTIME,*AREA,*QLATIN,*WT,*QVALUE,*AVALUE,
        *QWT,*DSDIST,USDIST[MAXSEG],
        DFACE[MAXSEG],HPLUSF[MAXSEG],HPLUSB[MAXSEG],
        HMULTF[MAXSEG],HMULTB[MAXSEG],HDIV[MAXSEG],
        HDIVF[MAXSEG],HADV[MAXSEG],GAMMA[MAXSEG],
        AFACE[MAXSEG],A[MAXSEG],C[MAXSEG],
        QINVAL[MAXFLOWLOC+1],USCONCN[MAXSOLUTE],QN[2],
        AREAN[MAXSEG],QLATINN[MAXSEG],AFACEN[2],
        GAMMAN[MAXSEG],AN[MAXSEG],CN[MAXSEG],
        AREA2[MAXSEG],ALPHA[MAXSEG];

    */

    char CHEM[12], STOPTYPE[12];
    // masud: 2D arrays in Fortran are actually flat arrays to be traversed as column major
    // declaring arrays with indexes is crashing the program -- because of stack overflow!
    // dynamic declaration (on the heap) will prevent crashes.
    double *USCONC, *CLATIN, *CONC, *CONC2, *CLVAL, *LAMBDA, *LAM2DT, *AWORK,
           *BWORK, *USBC, *SGROUP2, *SGROUP, *LHATDT, *SORB, *LHAT2DT, *KD, *CLATINN,
           *TWOPLUS, *BN, *TGROUP, *BTERMS, *BTERMSN, *IGROUP;

    int IDECAY; // <-- local variables maininit.f
    double XSTART,QSTART; //<-- local variables maininit.f
    double *PRTLOC, *DISP,   *FLOWLOC,  *QLATOUT,  *DIST,
           *LAMBDA2,  *LAMHAT,  *LAMHAT2,  *RHOLAM,  *CSBACK;  // <-- local variables from maininit.f

    // calloc dimension --> COL*ROW
    USCONC = (double*)calloc(MAXSOLUTE*MAXBOUND, sizeof(double));
    CLATIN = (double*)calloc(MAXSOLUTE*MAXSEG, sizeof(double));
    CONC = (double*)calloc(MAXSOLUTE*MAXSEG, sizeof(double));
    CONC2 = (double*)calloc(MAXSOLUTE*MAXSEG, sizeof(double));
    CLVAL = (double*)calloc(MAXSOLUTE*(MAXFLOWLOC+1), sizeof(double));
    LAMBDA = (double*)calloc(MAXSOLUTE*MAXSEG, sizeof(double));
    LAM2DT = (double*)calloc(MAXSOLUTE*MAXSEG, sizeof(double));
    AWORK = (double*)calloc(MAXSOLUTE*MAXSEG, sizeof(double));
    BWORK = (double*)calloc(MAXSOLUTE*MAXSEG, sizeof(double));
    USBC = (double*)calloc(MAXSOLUTE*MAXBOUND, sizeof(double));
    SGROUP2 = (double*)calloc(MAXSOLUTE*MAXSEG, sizeof(double));
    SGROUP = (double*)calloc(MAXSOLUTE*MAXSEG, sizeof(double));
    LHATDT = (double*)calloc(MAXSOLUTE*MAXSEG, sizeof(double));
    SORB = (double*)calloc(MAXSOLUTE*MAXSEG, sizeof(double));
    LHAT2DT = (double*)calloc(MAXSOLUTE*MAXSEG, sizeof(double));
    KD = (double*)calloc(MAXSOLUTE*MAXSEG, sizeof(double));
    CLATINN = (double*)calloc(MAXSOLUTE*MAXSEG, sizeof(double));
    TWOPLUS = (double*)calloc(MAXSOLUTE*MAXSEG, sizeof(double));
    BN = (double*)calloc(MAXSOLUTE*MAXSEG, sizeof(double));
    TGROUP = (double*)calloc(MAXSOLUTE*MAXSEG, sizeof(double));
    BTERMS = (double*)calloc(MAXSOLUTE*MAXSEG, sizeof(double));
    BTERMSN = (double*)calloc(MAXSOLUTE*MAXSEG, sizeof(double));
    IGROUP = (double*)calloc(MAXSOLUTE*MAXSEG,sizeof(double));

    // ---------------- Local variables from maininit.f --------------------
    PRTLOC = (double*)calloc(MAXPRINT, sizeof(double));
    DISP = (double*)calloc(MAXSEG, sizeof(double));
    FLOWLOC = (double*)calloc(MAXFLOWLOC, sizeof(double));
    QLATOUT = (double*)calloc(MAXSEG, sizeof(double));
    DIST = (double*)calloc(MAXSEG, sizeof(double));
    LAMBDA2 = (double*)calloc(MAXSEG*MAXSOLUTE, sizeof(double));
    LAMHAT = (double*)calloc(MAXSEG*MAXSOLUTE, sizeof(double));
    LAMHAT2 = (double*)calloc(MAXSEG*MAXSOLUTE, sizeof(double));
    RHOLAM = (double*)calloc(MAXSEG*MAXSOLUTE, sizeof(double));
    CSBACK = (double*)calloc(MAXSEG*MAXSOLUTE, sizeof(double));
    // ---------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
    /*
    printf("-----------  main.C  ---------------------------\n");

    printf("%d arguments passed: -->\n",argc);
    for (i = 0; i<argc; i++){
        printf("\t %s\n",argv[i]);
    }
    */

    /**********************************************************
                    Get parameters from ini file
    **********************************************************/
    dict = iniparser_load(argv[1]);
    if (dict == NULL) {
        printf("Error reading ini file\n");
		return 2; // ERROR_FILE_READ
	}
	//[OTIS]
    new_area = iniparser_getdouble(dict, "OTIS:AREA", -1.5);
    new_disp = iniparser_getdouble(dict, "OTIS:DISP", -1.5);
    new_area2 = iniparser_getdouble(dict, "OTIS:AREA2", -1.5);
    new_alpha = iniparser_getdouble(dict, "OTIS:ALPHA", -1.5);
    //[MCMC]
    mcmc_len = iniparser_getint(dict, "MCMC:mcmc_len", -1);
    burn_in = (int)(0.01*(float)mcmc_len*iniparser_getdouble(dict, "MCMC:burn_in(prct)", -1));
    mcmc_len++;
    std_obs = iniparser_getdouble(dict, "MCMC:obs_std", -1.5);
    pstd_A = iniparser_getdouble(dict, "MCMC:proposal_std_A", -1.5);
    pstd_D = iniparser_getdouble(dict, "MCMC:proposal_std_D", -1.5);
    pstd_As = iniparser_getdouble(dict, "MCMC:proposal_std_As", -1.5);
    pstd_Al = iniparser_getdouble(dict, "MCMC:proposal_std_Al", -1.5);
    seed = iniparser_getint(dict, "MCMC:rng_seed", -1);
    //[FILES]
    fname = iniparser_getstring(dict, "FILES:observation_file", "");
    fname2 = iniparser_getstring(dict, "FILES:mcmc_output_file", "mcmc_out.csv");
    //[DEBUG]
    debug_flag[0] = iniparser_getint(dict, "DEBUG:print_sample", 0);

    if (new_area <0 || new_disp<0 || new_area2<0 || new_alpha<0){
        printf("Error in initial parameter values provided in the ini file.\n");
        return 2; // ERROR_FILE_READ
    }
    //printf("INI DEBUG: %f, %f, %f, %f\n", new_area, new_disp, new_area2,new_alpha);
    //printf("INI DEBUG: %d, %d, %f, %f, %f, %f, %f\n", mcmc_len, burn_in, std_obs, pstd_A, pstd_D, pstd_As, pstd_Al);
    //printf("INI DEBUG: Obs_std = %f \n", std_obs);
    //printf("INI DEBUG: Seed = %d\n",seed);

/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
    chain_A = (double*)calloc(mcmc_len, sizeof(double));    chain_A[0] = new_area;
    chain_D = (double*)calloc(mcmc_len, sizeof(double));    chain_D[0] = new_disp;
    chain_As = (double*)calloc(mcmc_len, sizeof(double));   chain_As[0] = new_area2;
    chain_Al = (double*)calloc(mcmc_len, sizeof(double));   chain_Al[0] = new_alpha;
    chain_LL = (double*)calloc(mcmc_len, sizeof(double));

    flag_run = 0;
    mainold_(&IMAX,&NPRINT,&IPRINT,&NSOLUTE,&PRTOPT,&NFLOW,&JBOUND,&IBOUND,&NBOUND,&ISORB,
                 PINDEX, QINDEX,
                 &DSBOUND, &TSTART, &TFINAL, &TSTEP, &QSTEP, &TIMEB, &TSTOP,&QSTOP, &BCSTOP,
                 DELTAX,Q,USTIME,AREA,QLATIN,WT,
                 QVALUE,AVALUE,QWT,DSDIST,USDIST,
                 DFACE,HPLUSF,HPLUSB,HMULTF,HMULTB,HDIV,
                 HDIVF,HADV,GAMMA, AFACE,A,C,
                 QINVAL,USCONCN,QN,AREAN,QLATINN,AFACEN,
                 GAMMAN,AN,CN, AREA2,ALPHA,
                 USCONC, CLATIN, CONC, CONC2, CLVAL, LAMBDA, LAM2DT, AWORK,
                 BWORK, USBC, SGROUP2, SGROUP, LHATDT, SORB, LHAT2DT, KD, CLATINN,
                 TWOPLUS, BN, TGROUP, BTERMS, BTERMSN, IGROUP, CHEM, STOPTYPE,
                 &IDECAY, &XSTART, &QSTART, PRTLOC,  DISP,  FLOWLOC, QLATOUT, DIST, //<-- local variables maininit.f
                 LAMBDA2, LAMHAT, LAMHAT2, RHOLAM, CSBACK, //<-- local variables maininit.f
                 &ncout, &flag_run, &new_area, &new_disp, &new_area2, &new_alpha, &nreach, &nseg_tmp, &rchlen_tmp,
                 &time_start, sim_time, conc_main, conc_ts, test);

    if (ncout > MAXBOUND){
        printf("Print time step too small, output time points exceeding maximum\n specified limit (MAXBOUND) of %d\n",MAXBOUND);
        return 3;
    }
    printStep = (double)IPRINT * TSTEP;
    // printf("Start and End, printStep = %f, %f, %f\n", time_start, TFINAL, printStep);

   /**********************************************************
                Get observation from the observation file
    **********************************************************/
    printf("Observation File : %s\n",fname);
    printf("MCMC Output File : %s\n\n",fname2);
    fp_obs = fopen(fname,"r");
    fp_out = fopen(fname2,"w");

    read_flag = 1;
    count = 0;
    count_start = 0;
    count_final = 0;

    while(read_flag > 0){
        count++;
        read_flag = fscanf(fp_obs, "%lf %lf", &time_tmp, &conc1_tmp);
        //printf("read_flag, time, conc = %d, %f, %f\n", read_flag, time_tmp, conc1_tmp);
        if(time_tmp > TFINAL)
            break;
        if (time_tmp < time_start + printStep)
            count_start++;
    }
    count--;
    count_final = count;
    nobs = count_final - count_start;
    if((nobs) <= 0){
        printf("nobs, count_final, count_start = %d, %d, %d\n", nobs, count_final, count_start);
        printf("There is no observation in the simulation time range\n");
        printf("Qutting...\n");
        exit(1);    // NO_OBS --> enum ErrCode
    }
    printf("There are %d observations within the specified simulation time.\n\n", nobs);
    //printf("Obs starting index = %d\n", count_start);

    // check if observation time range matches the range of simulation time start and end
    if(nobs > MAXBOUND){
        printf("Number of observation greater than the specified limit (MAXBOUND = %d). \n\n", MAXBOUND);
        exit(1); // LIMIT_EXCEEDED --> enum ErrCode
    }

    obs_time = (double*)calloc(nobs+1,sizeof(double));
    obs_conc1 = (double*)calloc(nobs+1,sizeof(double));
    conc_interp = (double*)calloc(nobs+1,sizeof(double));   // will contain interpolated modeled concentration
    tind = (int*)calloc(nobs+1, sizeof(int));
    tfrac = (double*)calloc(nobs+1,sizeof(double));         // interpolation fraction

    rewind(fp_obs);
    ind_tmp = 0;
    for(i=0; i< count_final; i++){
        fscanf(fp_obs, "%lf %lf", &time_tmp, &conc1_tmp);
        if(i >= count_start){
            obs_time[ind_tmp] = time_tmp;
            obs_conc1[ind_tmp] = conc1_tmp;
            //printf("ind, time, conc = %d, %f, %f\n", ind_tmp, obs_time[ind_tmp], obs_conc1[ind_tmp]);
            ind_tmp++;
        }
    }
    //printf("%f \t %f\n\n", sim_time[0],sim_time[1]);

    ind_cur = 1;
    for (i = 0; i<nobs; i++){
        for(j=ind_cur; j<ncout; j++){
            if(obs_time[i] < sim_time[j]){
                //printf("obs, sim --- %f, %f\n",obs_time[i], sim_time[j]);
                ind_cur = j-1;
                tind[i] = j;
                tfrac[i] = (obs_time[i] - sim_time[j-1])/(sim_time[j] - sim_time[j-1]);
                //if(i<10)
                //    printf("time, conc, tind, tfrac = %f, %f, %d, %f\n", obs_time[i], obs_conc1[i], tind[i], tfrac[i]);
                break;
            }
        }
    }

    // -----------------------------------------------------------------------//
    // -----------------------------------------------------------------------//
    //seed = 2;   // set seed to -1 for seed = clock time
                // TODO: some better algorithm for random number generator should be used

    PlantSeeds(seed);
    //printf("double number = %f", exp(2.5));
    last_ll = -DBL_MAX;
    cur_ll = 0;
    tmp_ll = 0;
    accept = 0;
    reject = 0;
    flag_run = 1;

    printf("MCMC Iterations requested = %d\n",mcmc_len-1);
    printf("Iteration:  ");

    //FILE *fpDebugConc;
    //fpDebugConc = fopen("Debug_tfrac.csv","w");
    time_beginning = clock();

    for(mcmc_iter = 0; mcmc_iter < mcmc_len -1; mcmc_iter++){
        new_area = chain_A[mcmc_iter];
        new_disp = chain_D[mcmc_iter];
        new_area2 = chain_As[mcmc_iter];
        new_alpha = chain_Al[mcmc_iter];

        mainold_(&IMAX,&NPRINT,&IPRINT,&NSOLUTE,&PRTOPT,&NFLOW,&JBOUND,&IBOUND,&NBOUND,&ISORB,
                 PINDEX, QINDEX,
                 &DSBOUND, &TSTART, &TFINAL, &TSTEP, &QSTEP, &TIMEB, &TSTOP,&QSTOP, &BCSTOP,
                 DELTAX,Q,USTIME,AREA,QLATIN,WT,
                 QVALUE,AVALUE,QWT,DSDIST,USDIST,
                 DFACE,HPLUSF,HPLUSB,HMULTF,HMULTB,HDIV,
                 HDIVF,HADV,GAMMA, AFACE,A,C,
                 QINVAL,USCONCN,QN,AREAN,QLATINN,AFACEN,
                 GAMMAN,AN,CN, AREA2,ALPHA,
                 USCONC, CLATIN, CONC, CONC2, CLVAL, LAMBDA, LAM2DT, AWORK,
                 BWORK, USBC, SGROUP2, SGROUP, LHATDT, SORB, LHAT2DT, KD, CLATINN,
                 TWOPLUS, BN, TGROUP, BTERMS, BTERMSN, IGROUP, CHEM, STOPTYPE,
                 &IDECAY, &XSTART, &QSTART, PRTLOC,  DISP,  FLOWLOC, QLATOUT, DIST, //<-- local variables maininit.f
                 LAMBDA2, LAMHAT, LAMHAT2, RHOLAM, CSBACK, //<-- local variables maininit.f
                 &ncout, &flag_run, &new_area, &new_disp, &new_area2, &new_alpha, &nreach, &nseg_tmp, &rchlen_tmp,
                 &time_start, sim_time, conc_main, conc_ts, test);

        // interpolate simulated concentration in conc_interp[N] vector
        //printf("\ni, tind[i], tfrac[i], C1, C2, conc_interp[i]\n");

        for(i=0; i<nobs; i++){
            //c2 = conc_main[tind[i]]
            //c1 = conc_main[tind[i]-1]
            conc_interp[i] = conc_main[tind[i]]*tfrac[i] + conc_main[tind[i]-1]*(1-tfrac[i]);
            //printf("%d, %d, %f, %f, %f, %f\n",i, tind[i], tfrac[tind[i]], conc_main[tind[i]-1], conc_main[tind[i]], conc_interp[i]);
            //fprintf(fpDebugConc,"%f\n",tfrac[i]);
        }

        //fclose(fpDebugConc);

        likelihood(obs_conc1,conc_interp,nobs,std_obs,&tmp_ll);
        cur_ll = tmp_ll;

        //printf("main.c: curLL, lasLL = %f, %f, %f\n",cur_ll,last_ll,exp(cur_ll-last_ll));

        if(Uniform(0,1) < exp(cur_ll - last_ll)){
            /************ grow chain with newly accepted value *************/
            if(mcmc_iter > burn_in)
                accept++;
            last_ll = cur_ll;
            // store accepted modeled values -- not really
        }
        else{// roll-back to previous sample
            /************ grow chain with old values *************/
            if(mcmc_iter == 0){
                printf("Error in likelihood calculation\nQutting...");
                goto CLEANUP; // only for debug purpose;
            }
            if(mcmc_iter > burn_in)
                reject++;
            // For loop should crash if the very first sample is not accepted -- which also should not happen!!
            cur_ll = last_ll;
            //for i in range(n_cluster):
            //    patt_chain[chain_start[i] + i_mc] = patt_chain[chain_start[i] + i_mc - 1]
            chain_A[mcmc_iter] = chain_A[mcmc_iter-1];
            chain_D[mcmc_iter] = chain_D[mcmc_iter-1];
            chain_As[mcmc_iter] = chain_As[mcmc_iter-1];
            chain_Al[mcmc_iter] = chain_Al[mcmc_iter-1];
        }

        /************ grow likelihood chain *************/
        chain_LL[mcmc_iter] = cur_ll;

        /************ grow parameter chain with normal() *************/
        do{
            param_tmp = chain_A[mcmc_iter] + Normal(0,pstd_A);
        }while(param_tmp < 0);
        new_area = param_tmp;
        //new_area = (param_tmp > 0 ? param_tmp : chain_A[mcmc_iter]);

        do{
            param_tmp = chain_D[mcmc_iter] + Normal(0,pstd_D);
        }while(param_tmp < 0);
        new_disp = param_tmp;
        //new_disp = (param_tmp > 0 ? param_tmp : chain_D[mcmc_iter]);
        do{
            param_tmp = chain_As[mcmc_iter] + Normal(0,pstd_As);
        }while(param_tmp < 0);
        new_area2 = param_tmp;
        //new_area2 = (param_tmp > 0 ? param_tmp : chain_As[mcmc_iter]);

        do{
            param_tmp = chain_Al[mcmc_iter] + Normal(0,pstd_Al);
        }while(param_tmp < 0);
        new_alpha = param_tmp;
        //new_alpha = (param_tmp > 0 ? param_tmp : chain_Al[mcmc_iter]);

        //Grow Chains:
        chain_A[mcmc_iter + 1] = new_area;
        chain_D[mcmc_iter + 1] = new_disp;
        chain_As[mcmc_iter + 1] = new_area2;
        chain_Al[mcmc_iter + 1] = new_alpha;

        if((mcmc_iter + 1) % 500 == 0)
            printf("%d  ",(int)((mcmc_iter+1)));
    } //for(mcmc_iter = 0; mcmc_iter<mcmc_len; mcmc_iter++)
    printf("\n\n");

    time_end = clock();
    time_mcmc = (double)(time_end - time_beginning) / CLOCKS_PER_SEC;
    printf("Total time required for MCMC = %.2f hrs\n",time_mcmc/3600);

    printf("Acceptance rate: %.2f\n\n",((float)accept/(float)mcmc_len*100));

    fprintf(fp_out, "Likelihood, Area, Disp, Area2, Alpha\n");
    for(i=0; i<mcmc_len -1; i++)
        fprintf(fp_out,"%f,%f,%f,%f,%f\n",chain_LL[i],chain_A[i],chain_D[i],chain_As[i],chain_Al[i]);

    // ---------------DEBUG BLOCK-------------------------------------------------
    if(debug_flag[0]){
        fp_debug = fopen("Debug_ConcOut.csv","w");
        nline = (int)((TFINAL - 0)/TSTEP/IPRINT) + 10;
        fprintf(fp_debug,"Printing iteration %d\n",mcmc_iter-2);
        fprintf(fp_debug,"A,D,As,Al\n");
        fprintf(fp_debug,"%f,%f,%f,%f\n",chain_A[mcmc_iter-2],chain_D[mcmc_iter-2],chain_As[mcmc_iter-2],chain_Al[mcmc_iter-2]);
        fprintf(fp_debug,"\nTime,Conc\n");
        for(i=0; i<nline; i++){
            fprintf(fp_debug,"%f, %f\n",sim_time[i],conc_main[i]);
        }
        fclose(fp_debug);
    }
    // ---------------end of DEBUG BLOCK -------------------------------------------------

    CLEANUP:
    fclose(fp_obs); fclose(fp_out);
    free(test); free(sim_time); free(conc_main);  free(conc_ts);
    free(obs_time); free(obs_conc1); free(conc_interp); free(tind); free(tfrac);

    printf("Freeing up dynamic variables...");

    free(PRTLOC); free(DISP);   free(FLOWLOC);  free(QLATOUT);  free(DIST);
    free(LAMBDA2);  free(LAMHAT);  free(LAMHAT2);  free(RHOLAM);  free(CSBACK);

    free(chain_A); free(chain_D); free(chain_As); free(chain_Al); free(chain_LL);

    free(USCONC); free(CLATIN); free(CONC); free(CONC2); free(CLVAL); free(LAMBDA);
    free(LAM2DT); free(AWORK); free(BWORK); free(USBC); free(SGROUP2); free(SGROUP);
    free(LHATDT); free(SORB); free(LHAT2DT); free(KD); free(CLATINN); free(TWOPLUS);
    free(BN); free(TGROUP); free(BTERMS); free(BTERMSN); free(IGROUP);

    printf(" Success!\n\nExit.\n");
    exit(0);

    return 0;
}


double likelihood(double *obs, double *mod, int nobs, double sigma, double *ll)
{
    //FILE * fptmp;
    //fptmp = fopen("Debug_LL.csv","w");
    //fprintf(fptmp,"Sigma,%f,nobs,%d\n",sigma,nobs);
    int i;
    *ll = 0;
    for(i=0; i<nobs; i++){
        *ll += -0.5*pow(((obs[i] - mod[i])/sigma),2);
        //fprintf(fptmp,"%f,%f,%f\n",obs[i],mod[i],*ll);
    }

    //fclose(fptmp);
    *ll = *ll/(float)nobs;
    //printf("DEBUG: Likelihood = %f\n", *ll);
    return 0;
}
