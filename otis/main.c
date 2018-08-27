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
#define MAXBOUND    20000
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

double getNormal(double mu, double sigma);
double likelihood(double *obs, double *mod, int nobs, double sigma, double *ll);

// count frequency to draw histogram of modeled output
int updateFreq(int nobs, int *bin_indexes, int *freq);
int countFrequency(double *obs, double *mod, int nobs, float bin_width, int nbins, int *bin_indexes, int *freq);
//int printModelConc(FILE* fpModel, double* conc, int ncout, int nchains);
int printModelConc2(FILE* fpModel, double* conc, int ncout);
int makeEcho(FILE *fp, char* config_file, int mcmcLen, int burnIn, int rngSeed, int nobservation, char* concChainFile,
             float pAccept, float elapsedTime, int* debugMode, double pstdD, double pstdA, double pstdAs, double pstdAl,
             double stdObs, int nobs, double* obsT, double* obsC);

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
int read_obsFile(char* fname, int *nobs, float* time, float* conc);

int main(int argc, char** argv ){
    int i, j, count, ind_tmp, ind_cur, count_start, count_final;
    int mcmc_iter, mcmc_len, burn_in, ncout, nobs, nreach, nseg_tmp, accept, reject, flag_run;
    int nline, nbins;
    int *tind, *freq, *bin_indexes;
    int saveModelConc;      // flag: save chains of models for creating confidence bounds
    int lineLength = 100, nobs_tmp;   // for reading observation file

    //float binWidth;         // used to make histogram. pretty useless, need to remove it.

    double rchlen_tmp, time_start, printStep;
    double new_area, new_disp, new_area2, new_alpha, time_tmp, conc1_tmp;
    double *sim_time, *obs_time, *conc_main, *obs_conc1, *conc_ts, *conc_interp, *tfrac;
    double *chain_A, *chain_As, *chain_D, *chain_Al, *chain_LL;
    double tmp_ll, cur_ll, last_ll, pstd_A, pstd_D, pstd_As, pstd_Al, std_obs, time_mcmc;  // pstd --> proposal standard deviation
    double *lastConc;
    double *test;   // This variable is only for debug only.

    long seed;

    char buffer[101];   // for reading observation file.
    char *fname, *fname2, *fModConc;
    clock_t time_beginning, time_end;



    FILE *fp_obs, *fp_out, *fp_debug, *fpEcho, *fpModels;
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
    printf("\n\n\n\t\t------------------------------\n"
           "\t\t    OTIS-MC  (version 0.1)  \n"
           "\t\t     (ranasd@mail.uc.edu)   \n"
           "\t\t         Oct 06, 2017   \n"
           "\t\t------------------------------\n\n");
    //printf("=======================================================================\n");

    dict = iniparser_load(argv[1]);
    if (dict == NULL) {
        printf("  Error reading ini file\n");
        printf("  Run OTIS-MC.exe with the configuration file.\n");
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
    nbins = iniparser_getint(dict, "MCMC:nbins", 13);                       // number of bins for generating the histogram of modeled concentration (around the observation concentration)
    //binWidth = (float)iniparser_getdouble(dict, "MCMC:binWidth", 0.1);      // width of bin (in units of concentration of solute being modeled)
    //[FILES]
    fname = iniparser_getstring(dict, "FILES:observation_file", "");
    fname2 = iniparser_getstring(dict, "FILES:mcmc_output_file", "mcmc_out.csv");
    //[DEBUG]
    debug_flag[0] = iniparser_getint(dict, "DEBUG:print_sample", 0);
    saveModelConc = iniparser_getint(dict, "DEBUG:saveModelConc", 0);

    if (new_area <0 || new_disp<0 || new_area2<0 || new_alpha<0){
        printf("  Error in initial parameter values provided in the ini file.\n");
        return 2; // ERROR_FILE_READ
    }
    //printf("INI DEBUG: %f, %f, %f, %f\n", new_area, new_disp, new_area2,new_alpha);
    //printf("INI DEBUG: %d, %d, %f, %f, %f, %f, %f\n", mcmc_len, burn_in, std_obs, pstd_A, pstd_D, pstd_As, pstd_Al);
    //printf("INI DEBUG: Obs_std = %f \n", std_obs);
    //printf("INI DEBUG: Seed = %d\n",seed);
    debug_flag[1] = saveModelConc;
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
        printf("  Print time step too small: No. of output time points = %d\n  Exceeding maximum specified limit (MAXBOUND) of %d\n",ncout, MAXBOUND);
        return 3;
    }
    printStep = (double)IPRINT * TSTEP;
    // printf("Start and End, printStep = %f, %f, %f\n", time_start, TFINAL, printStep);

   /**********************************************************
                Get observation from the observation file
    **********************************************************/
    printf("\tObservation File : %s\n",fname);
    printf("\tMCMC Output File : %s\n\n",fname2);
    fp_obs = fopen(fname,"r");
    fp_out = fopen(fname2,"w");

    if(fp_obs == NULL){
        printf("Error opening observation file: %s\nUnsuccessful Exit...\n\n",fname);
        return -1;
    }

    while(fgets(buffer, lineLength, fp_obs) != NULL){
        if(buffer[0] == '#'){
            continue;
        }
        else{   // read the number of observations
            if(sscanf(buffer, "%d", &nobs_tmp) != 1){
                printf("Error reading NUMBER of OBSERVATIONS from file: %s\nUnsuccessful Exit...\n\n",fname);
                return -1;
            }
            printf("\tNumber of Observations: %d\n", nobs_tmp);
            break;
        }
    }

    count = 0;
    count_start = 0;
    count_final = 0;
    /* ------------------------------------------------------------------- */
    // Masud (8-28-2017): Improved version of observation file reading:
    for(i=0; i<nobs_tmp; i++){
        count++;
        if(fgets(buffer, lineLength, fp_obs) == NULL){
            printf("Error reading TIME -- CONC pairs from file: %s\nUnsuccessful Exit...\n\n",fname);
            printf("Number of lines read: %d\n", i);
            exit(1);
        }
        if(sscanf(buffer, "%lf %lf", &time_tmp, &conc1_tmp) != 2){
            printf("Error reading TIME -- CONC pairs from file: %s\n"
                   "Number of lines read: %d\n"
                   "Unsuccessful Exit...\n\n",fname,i);
            exit(1);
        }
        if(time_tmp > TFINAL)
            break;
        if (time_tmp < time_start + printStep)
            count_start++;
    }
    /* ------------------------------------------------------------------- /
    // Old version of observation file reading:
    while(read_flag > 0){
        count++;
        read_flag = fscanf(fp_obs, "%lf %lf", &time_tmp, &conc1_tmp);
        if(time_tmp > TFINAL)
            break;
        if (time_tmp < time_start + printStep)
            count_start++;
    }
    / ------------------------------------------------------------------- */
    count--;
    count_final = count;
    nobs = count_final - count_start;       // actual number of observation within the simulation range
    /* ------------------------------------------------------------------- */
    if((nobs) <= 0){
        printf("nobs, count_final, count_start = %d, %d, %d\n", nobs, count_final, count_start);
        printf("There is no observation in the simulation time range\n");
        printf("Qutting...\n");
        exit(1);    // NO_OBS --> enum ErrCode
    }
    printf("\tThere are %d observations within the specified simulation time.\n\n", nobs);
    //printf("Obs starting index = %d\n", count_start);

    // check if observation time range matches the range of simulation time start and end
    if(nobs > MAXBOUND){
        printf("\tNumber of observation greater than the specified limit (MAXBOUND = mcmc_len%d). \n\n", MAXBOUND);
        exit(1); // LIMIT_EXCEEDED --> enum ErrCode
    }

    obs_time = (double*)calloc(nobs+1,sizeof(double));
    obs_conc1 = (double*)calloc(nobs+1,sizeof(double));
    conc_interp = (double*)calloc(nobs+1,sizeof(double));   // will contain interpolated modeled concentration
    tind = (int*)calloc(nobs+1, sizeof(int));
    tfrac = (double*)calloc(nobs+1,sizeof(double));         // interpolation fraction

    freq = (int*)calloc(nbins*nobs, sizeof(int));       // will contain histogram frequency of the modeled concentration
    bin_indexes = (int*)calloc(nobs, sizeof(int));      // will contain the indexes of the previously updated bins

    fModConc = "(Not Saved)";
    if(saveModelConc){
        // This block is for saving accepted MCMC modeled concentration chains for plotting purpose
        fModConc = "modelConc.csv";
        //fprintf(fpEcho, "\nModel concentration output written to: \"%s\"\n", "modelConc.csv");

        fpModels = fopen(fModConc,"w");

        fprintf(fpModels, "chain number, %d\n", mcmc_len - 1 - burn_in);
        fprintf(fpModels, "%f", sim_time[0]);
        for(i=1; i<ncout; i++)
            fprintf(fpModels, ",%f", sim_time[i]);
        fprintf(fpModels, "\n");
    }

    lastConc = (double*) calloc(ncout, sizeof(double));
    for(i=0; i<ncout; i++)
        lastConc[i] = conc_main[i];

    /* ----- Read and Store the Observation from file ----------------------------------- */
    rewind(fp_obs);
    // go through all the comments at the head of the observation file and read the first
    // line that has the number of observations:
    while(fgets(buffer, lineLength, fp_obs) != NULL){
        if(buffer[0] == '#'){
            continue;
        }
        else{   // read the number of observations and then break out of loop
            break;
        }
    }

    // Read Time-Conc pairs:
    ind_tmp = 0;
    for(i=0; i< count_final; i++){
        fscanf(fp_obs, "%lf %lf", &time_tmp, &conc1_tmp);
        if(i >= count_start){
            obs_time[ind_tmp] = time_tmp;
            obs_conc1[ind_tmp] = conc1_tmp;
            ind_tmp++;
        }
    }
    /* ------------------------------------------------------------------- */
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
    //seed = 2;   // set seed to -1 for seed = clock time
                // TODO: some better algorithm for random number generator should be used

    PlantSeeds(seed);   // a function from rngs.c
    //printf("double number = %f", exp(2.5));
    last_ll = -DBL_MAX;
    cur_ll = 0;
    tmp_ll = 0;
    accept = 0;
    reject = 0;
    flag_run = 1;

    printf("\tMCMC Iterations requested = %d\n",mcmc_len-1);
    printf("\tIteration:  ");
    time_beginning = clock();
    /* ------------------------------------------------------------------------------------------------------ */
    /* ------------------------------------------ MCMC Begins ----------------------------------------------- */
    /* ------------------------------------------------------------------------------------------------------ */
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
        // ********************************************************************************************** //
        // Masud (9-25-2017):
        // For GLUE implementation the if--else block should be turned off
        // and all samples should be accepted
        // ********************************************************************************************** //

        if(Uniform(0,1) < exp(cur_ll - last_ll)){
            /************ grow chain with newly accepted value *************/
            if(mcmc_iter >= burn_in){
                accept++;
                if(saveModelConc){
                    for(i=0; i<ncout; i++)
                        lastConc[i] = conc_main[i];     // save last accepted model concentrations for rejected parameters
                    printModelConc2(fpModels,lastConc,ncout);
                }

            }
            last_ll = cur_ll;
        }
        else{// roll-back to previous sample
            /************ grow chain with old values for rejected parameter samples *************/
            if(mcmc_iter >= burn_in){
                reject++;
                if(saveModelConc){
                    printModelConc2(fpModels,lastConc,ncout);
                }
            }

            // For loop should crash if the very first sample is not accepted -- which also should not happen!!
            cur_ll = last_ll;

            chain_A[mcmc_iter] = chain_A[mcmc_iter-1];
            chain_D[mcmc_iter] = chain_D[mcmc_iter-1];
            chain_As[mcmc_iter] = chain_As[mcmc_iter-1];
            chain_Al[mcmc_iter] = chain_Al[mcmc_iter-1];
        }

        /************ grow likelihood chain *************/
        chain_LL[mcmc_iter] = cur_ll;

        /************ grow parameter chain with normal() *************/
        // Masud (9-25-2017):
        // For GLUE implementation the next four DO loop should be replaced with
        // Uniform distribution instead of Normal

        /*
        // Original code:
        new_area = chain_A[mcmc_iter]   + Normal(0,pstd_A);
        new_disp = chain_D[mcmc_iter]   + Normal(0,pstd_D);
        new_area2 = chain_As[mcmc_iter] + Normal(0,pstd_As);
        new_alpha = chain_Al[mcmc_iter] + Normal(0,pstd_Al);
        */
        // if proposal std is zero (or less) then that parameter is held constant
        new_area = chain_A[mcmc_iter]   + getNormal(0.0, pstd_A);
        new_disp = chain_D[mcmc_iter]   + getNormal(0.0, pstd_D);
        new_area2 = chain_As[mcmc_iter] + getNormal(0.0, pstd_As);
        new_alpha = chain_Al[mcmc_iter] + getNormal(0.0, pstd_Al);

        while (new_area<0.0 || new_disp<0.0 || new_area2<0.0 || new_alpha<0.0){
            mcmc_iter++;    // move the chain iterator forward
            if(mcmc_iter >= burn_in){
                reject++;   // reject the negative sample
                if(saveModelConc){
                    printModelConc2(fpModels,lastConc,ncout);  // print concentrations for the last accepted parameters
                }
            }
            // roll back parameter samples and likelihood chain
            cur_ll = last_ll;
            chain_A[mcmc_iter] = chain_A[mcmc_iter-1];
            chain_D[mcmc_iter] = chain_D[mcmc_iter-1];
            chain_As[mcmc_iter] = chain_As[mcmc_iter-1];
            chain_Al[mcmc_iter] = chain_Al[mcmc_iter-1];
            chain_LL[mcmc_iter] = cur_ll;

             /************ Try sampling again *************/
            new_area = chain_A[mcmc_iter]   + getNormal(0.0, pstd_A);
            new_disp = chain_D[mcmc_iter]   + getNormal(0.0, pstd_D);
            new_area2 = chain_As[mcmc_iter] + getNormal(0.0, pstd_As);
            new_alpha = chain_Al[mcmc_iter] + getNormal(0.0, pstd_Al);
        }

        //Propose new parameter values:
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

    printf("\tTotal time required for MCMC = %.2f hrs\n",time_mcmc/3600);
    printf("\tAcceptance rate: %.2f\n\n",((float)accept/(float)mcmc_len*100));


    // print parameter chains to file
    fprintf(fp_out, "Likelihood, Area, Disp, Area2, Alpha\n");
    for(i=0; i<mcmc_len -1; i++)
        fprintf(fp_out,"%f,%f,%f,%f,%f\n",chain_LL[i],chain_A[i],chain_D[i],chain_As[i],chain_Al[i]);


    fpEcho = fopen("echo-mc.out", "w");

    makeEcho(fpEcho, argv[1], mcmc_len-1, (int)burn_in, (int)seed, nobs,fModConc,((float)accept/(float)mcmc_len*100),
             time_mcmc/3600.0, debug_flag, pstd_D, pstd_A, pstd_As, pstd_Al, std_obs, nobs, obs_time, obs_conc1);

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

    printf("\tFreeing allocated memory ...");

    //CLEANUP:
    fclose(fp_obs); fclose(fp_out); fclose(fpEcho);

    free(test); free(sim_time); free(conc_main);  free(conc_ts);
    free(obs_time); free(obs_conc1); free(conc_interp); free(tind); free(tfrac);

    if(saveModelConc){
        fclose(fpModels);
    }
    free(lastConc);


    free(PRTLOC); free(DISP);   free(FLOWLOC);  free(QLATOUT);  free(DIST);
    free(LAMBDA2);  free(LAMHAT);  free(LAMHAT2);  free(RHOLAM);  free(CSBACK);

    free(chain_A); free(chain_D); free(chain_As); free(chain_Al); free(chain_LL);
    free(freq); free(bin_indexes);

    free(USCONC); free(CLATIN); free(CONC); free(CONC2); free(CLVAL); free(LAMBDA);
    free(LAM2DT); free(AWORK); free(BWORK); free(USBC); free(SGROUP2); free(SGROUP);
    free(LHATDT); free(SORB); free(LHAT2DT); free(KD); free(CLATINN); free(TWOPLUS);
    free(BN); free(TGROUP); free(BTERMS); free(BTERMSN); free(IGROUP);

    printf(" Success!\n");
    printf("  ==================================================================================\n");
    exit(0);

    return 0;
}

double getNormal(double mu, double sigma){
    if (sigma<=0.0)
        return 0.0;
    return Normal(mu, sigma);
}

/*
int printModelConc(FILE* fpModel, double* conc, int ncout, int nchains)
{
    int i, chain, istart;
    for(chain = 0; chain < nchains; chain++){
        istart = chain * ncout;
        for(i=0; i<ncout; i++)
            fprintf(fpModel, "%f, ",conc[istart + i]);
        fprintf(fpModel, "\n");
    }
    return 0;
}
*/

int printModelConc2(FILE* fpModel, double* conc, int ncout)
{
    int i;

    fprintf(fpModel, "%f",conc[0]);     // First concentration is printed separately to prevent the last comma
    for(i=1; i<ncout; i++)
            fprintf(fpModel, ",%f",conc[i]);
    fprintf(fpModel, "\n");

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

int updateFreq(int nobs, int *bin_indexes, int *freq)
{
    // this function will be used to update bin frequency when a particular parameter set is rejected
    //      and a rollback to the previous parameter set is performed.
    int i;

    for (i=0; i<nobs; i++)
    {
        //printf("updateFreq(): bin_indexes[%d] = %d", i, bin_indexes[i]);
        freq[bin_indexes[i]]++;
    }


    return 0;
}

int countFrequency(double *obs, double *mod, int nobs, float bin_width, int nbins, int *bin_indexes, int *freq){
/*
    Counts the frequency of model values in bins around the observed value for each accepted parameter sets
    from the MCMC simulations so that model prediction bounds can be drawn.
*/
    int iobs;
    int binInd, bin_upper, bin_lower;         // index of the bin, the center being zero
    int binCenter_t, bin_center;
    double obs_t;       // value for the bin containing the obs
    const float TOLERANCE = 0.00001;

    bin_center = (nbins - 1)/2;        // index of the mid point of freq[] which contains the observation
    for(iobs = 0; iobs < nobs; iobs++)
    {
        bin_lower = nbins * iobs;
        binCenter_t = bin_center + bin_lower;
        bin_upper = nbins * (iobs + 1) - 1;

        //printf("bin_lower, binCenter_t, bin_upper = %d, %d, %d\n", bin_lower, binCenter_t, bin_upper);

        if(mod[iobs] > obs[iobs])
        {
            obs_t = obs[iobs] - bin_width/2.0;
            binInd = (int)ceil((mod[iobs] - obs_t)/bin_width) - 1;
            //printf("Mod = %f, binInd = %d\n", mod[i], binInd);
            binInd += binCenter_t;                              // actual index in freq[]
            if (binInd > bin_upper)                         // if modeled value is too high than observed
            {
                freq[bin_upper] += 1;
                bin_indexes[iobs] = bin_upper;
            }
            else
            {
                freq[binInd] += 1;
                bin_indexes[iobs] = binInd;
            }

        }
        else
        {
            obs_t = obs[iobs] + bin_width/2.0;
            binInd = -1 * (int)floor((obs_t - mod[iobs])/bin_width + TOLERANCE) ;  // bindInd will be negative
            //printf("Mod = %f, binInd = %d\n", mod[i], binInd);
            binInd += binCenter_t;                              // actual index in freq[]
            if (binInd < bin_lower)                                 // if modeled value is too low than observed
            {
                freq[bin_lower] += 1;
                bin_indexes[iobs] = bin_lower;
            }
            else
            {
                freq[binInd] += 1;
                bin_indexes[iobs] = binInd;
            }
        }
    }
    return 0;
}


int makeEcho(FILE *fp, char* config_file, int mcmcLen, int burnIn, int rngSeed, int nobservation, char* concChainFile,
             float pAccept, float elapsedTime, int* debugMode, double pstdD, double pstdA, double pstdAs, double pstdAl,
             double stdObs, int nobs, double* obsT, double* obsC)
{
    // -- Local Variables -------------------------------------- //
    int i;
    int debug = 0;

    time_t timer;
    char buffer[26];
    struct tm* currentTime;
    // --------------------------------------------------------- //

    time(&timer);
    currentTime = localtime(&timer);
    strftime(buffer, 26, "%m-%d-%Y (%H:%M:%S)", currentTime);

    for (i=0; i<10; i++){
        if(debugMode[i]){
            debug = 1;
            break;
        }
    }

    fprintf(fp, "\n\n");
    fprintf(fp, "\t\t\t\t\t  MARKOV CHAIN MONTE CARLO PARAMETER ESTIMATION WITH\n");
    fprintf(fp, "\t\t\t\t\t  OTIS SOLUTE TRANSPORT MODEL (Runkel 1998)\n\n");
    fprintf(fp, "\t\t\t\t\t  VERSION: 0.1 (2017)\n\n");
    fprintf(fp, "\t\t\t\t\t  S.M. Masud Rana (sm.masudrana@gmail.com)\n");
    fprintf(fp, "\t\t\t\t\t  Aug 21, 2017\n");
    fprintf(fp, "\t\t\t\t\t  --------------------------------------------------\n\n");

    fprintf(fp, "\n\t\t  DATE: %s\n", buffer);

    fprintf(fp,"\n\n");
    fprintf(fp, "\t\t  INPUT\n\t\t  --------------------------------------------------------------\n");
    fprintf(fp, "\t\t  MCMC CONFIGURATION FILE                  :   %s\n", config_file);
    fprintf(fp, "\t\t  MCMC ITERATIONS                          :   %d\n", mcmcLen);
    fprintf(fp, "\t\t  BURN-IN                                  :   %d\n", burnIn);
    fprintf(fp, "\t\t  RANDOM NO GENERATOR SEED (-1 CLOCK TIME) :   %d\n", rngSeed);
    fprintf(fp, "\t\t  NO. D/S OBSERVATION                      :   %d\n", nobservation);
    fprintf(fp, "\t\t  MEASUREMENT STD (NORMAL DIST)            :   %.3f\n", stdObs);
    fprintf(fp, "\t\t  PROPOSAL STD (NORMAL DIST)\n");
    fprintf(fp, "\t\t  (Dispersion, Area, TS Area, Alpha)       :   %.2E, %.2E, %.2E, %.2E\n", pstdD, pstdA, pstdAs, pstdAl);
    if(debug)
        fprintf(fp, "\n\t\t  DEBUG MODE                               :   %s\n", "ON");

    fprintf(fp, "\n\n");
    fprintf(fp, "\t\t  OUTPUT\n\t\t  --------------------------------------------------------------\n");
    fprintf(fp, "\t\t  ACCEPTANCE RATE (%%)                      :   %.2f\n", pAccept);
    fprintf(fp, "\t\t  TOTAL TIME ELAPSED (Hours)               :   %.3f\n", elapsedTime);

    if(debugMode[0])
        fprintf(fp, "\t\t  CONCENTRATIONS FOR LAST ACCEPTED PARAMS  :   %s\n", "Debug_ConcOut.csv");
    if(debugMode[1])
        fprintf(fp, "\t\t  CONCENTRATION CHAINS WRITTEN TO          :   %s\n", concChainFile);

    fprintf(fp, "\n\n");
    fprintf(fp, "\t\t  OBSERVATIONS USED\n\t\t  --------------------------------------------------------------\n");
    fprintf(fp, "\t\t  TIME            CONC\n\t\t  --------        --------\n");

    for(i=0; i<nobs; i++)
        fprintf(fp, "\t\t   %7.4f         %5.2f\n", obsT[i], obsC[i]);

    return 0;
}

int read_obsFile(char* fname, int *nobs, float* time, float* conc){
    FILE * fp;
    //char* fname = "obsE12.inp";
    char buffer[201];
    int lineLength = 200, i;
    //float time[500], conc[500];

    fp = fopen(fname,"r");

    if(fp == NULL){
        printf("Error opening observation file: %s\nUnsuccessful Exit...\n\n",fname);
        return -1;
    }

    while(fgets(buffer,lineLength,fp) != NULL){
        if(buffer[0] == '#'){
            continue;
        }
        else{   // read the number of observations
            if(sscanf(buffer, "%d", nobs) != 1){
                printf("Error reading NUMBER of OBSERVATIONS from file: %s\nUnsuccessful Exit...\n\n",fname);
                return -1;
            }
            break;
        }
    }

    //printf("\n\nNumber of observations to read: %d\n", *nobs);

    for(i=0; i<(*nobs); i++){
        if(fgets(buffer,lineLength,fp) == NULL){
            printf("Error reading TIME -- CONC pairs from file: %s\nUnsuccessful Exit...\n\n",fname);
            return -1;
        }
        if(sscanf(buffer, "%f %f", &time[i], &conc[i]) != 2){
            printf("Error reading TIME -- CONC pairs from file: %s\n"
                   "Number of lines read: %d\n"
                   "Unsuccessful Exit...\n\n",fname,i);
            return -1;
        }
    }

    fclose(fp);
    return 0;
}


