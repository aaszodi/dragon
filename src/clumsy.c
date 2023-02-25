/* ==== PROGRAM clumsy.c ==== */

/* Clusters a set of PDB structures together.
 * All chains must have the same length. All atoms or C-alphas
 * are considered. Uniform weighting only
 */

/* ANSI C, 30. Apr. 1998. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* ---- INCLUDE FILES ---- */

#include "cmdopt.h" /* command line parameter parser */
#include "rotpdb.h"    /* PDB I/O and rotation */
#include "dslclu.h"    /* metric distance single-linkage clustering */

/* ---- DEFINITIONS ---- */

#define WINLEN 3    /* default C-alpha trace smoothing window length */
#define SMCYC 5	    /* default no. of smoothing cycles */

/* ---- PROTOTYPES ---- */

static void avg_str(const char *Outfnm, const Dslclu_ *Clu, double ***Structs, 
	int Size, char **Snames, const Pdbentry_ *Pdbdescr);
static double **cluster_avg(const Dslclu_ *Clu, double ***Structs,
	int Size);

static void avg_sd(double Data[], int Datno, double *Avg, double *Sd);
static void print_help(const char *Progname);

/* ==== MAIN ==== */

int main(int argc, char *argv[])
{
    double ***Structs=NULL;
    Trimat_ Rms=NULL;
    int Winlen=WINLEN, Smcyc=SMCYC, Fileno, Fx, 
	Size=0, Structsize=0, i, j, Fno, Sno, 
	Nok=0, Nbd=0, Nxx=0;
    Dslclu_ *Clus=NULL;
    double Rmsavg=0.0, Rmssd=0.0, Rmsmin=HUGE_VAL, Rmsmax=0.0;
    char **Snames=NULL, **Remarks=NULL;
    char *Name, *Outfnm=NULL, *Outfnm2=NULL;
    char Smooth=0, Allatoms=0, Rch;
    Pdbentry_ *Pdbdescr=NULL, *Pdbout=NULL;
    
    /* parse the command line */
    parse_optstr("as w%d<window_len> c%d<smooth_cycno> o%s<output>");
    Fx=get_options(argc, argv);
    if (Fx>=argc-1)
    {
	fprintf(stderr, "! %s: Please specify at least two structures\n", argv[0]);
	print_help(argv[0]);
	exit(EXIT_FAILURE);
    }
    
    /* decide if all-atom or C-alpha superimpositions will be made:
     * no smoothing for all-atom comparisons
     */
    Allatoms=optval_bool('a');
    Smooth=Allatoms? 0: optval_bool('s');
    
    /* set up storage for the structures and their names */
    Structs=(double ***) calloc(argc-Fx, sizeof(double **));
    Snames=(char **) calloc(argc-Fx, sizeof(char *));
    
    /* Get the structures. Lengths must be equal */
    for (Sno=Fileno=0; Fileno<argc-Fx; Fileno++, Sno++)
    {
	/* read the 1st chain from the next file
	 * and store the original PDB layout in Pdbdescr
	 * from the first structure (when Pdbdescr==NULL)
	 */
	Structs[Sno]=get_vectors(argv[Fx+Fileno], &Pdbdescr, Allatoms, NULL, &Structsize);
	if (Structs[Sno]==NULL || Structsize<=0)
	{
	    fprintf(stderr, "\n? %s: Cannot process %s\n", 
		argv[0], argv[Fx+Fileno]);
	    Sno--; continue;
	}
	
	if (!Sno) Size=Structsize;  /* save first length */
	else if (Size!=Structsize)
	{
	    fprintf(stderr, "\n? %s: Structure from \"%s\": size mismatch (%d!=%d)\n", 
		    argv[0], argv[Fx+Fileno], Structsize, Size);
	    Sno--; continue;
	}
	
	/* store filename if struct is OK */
	Snames[Sno]=(char *) calloc(strlen(argv[Fx+Fileno])+1, sizeof(char));
	strcpy(Snames[Sno], argv[Fx+Fileno]);
    }
    
    if (Sno<2)
    {
	fprintf(stderr, "\n! %s: Too few (%d) valid input structures, exiting...\n",
	    argv[0], Sno);
	exit(EXIT_FAILURE);
    }
    
    /* get the remaining options */
    if (optval_str('o', &Outfnm))
	Outfnm2=(char *) calloc(strlen(Outfnm)+20, sizeof(char));
    
    /* set smoothing parameters if smoothing is on: replace by
     * default values if not specified or the values are silly
     */
    if (Smooth)
    {
	if (!optval_int('w', &Winlen) ||
		Winlen<=0 || Winlen>Size/8)
	    Winlen=WINLEN;
	if (!optval_int('c', &Smcyc) ||
		Smcyc<=0 || Smcyc>Size-2)
	    Smcyc=SMCYC;
    }
    
    printf("# PDB %s%s rigid body superposition: metric single-linkage clustering: %s\n", 
	    (Smooth)? "smoothed ":"", 
	    (Allatoms)? "all-atom":"C-alpha", 
	    argv[0]);
	    
    /* smooth structs and weights */
    if (Smooth)
    {
	printf("# Window length=%d, no. of smooth cycles=%d\n", 
		    Winlen, Smcyc);
	for (i=0; i<Sno; i++)
	    Structs[i]=smooth_chains(Structs[i], &Size, Winlen, Smcyc);
    }
    
    /* generate pairwise RMS values */
    Rms=alloc_trimat(Sno);
    for (i=0; i<Sno; i++)
	for (j=0; j<i; j++)
	    Rms[i][j]=rotate_vectors(Structs[i], Structs[j], NULL, Size);
    puts("# The RMS matrix:");
    list_trimat(Rms, Sno, 80, 4, 1);
    puts("# List of structures:");
    for (i=0; i<Sno; i++)
	printf("[%d] %s\n", i, Snames[i]);
    
    /* do the clustering */
    Clus=make_dslclus(Rms, Sno);
    if (Clus==NULL)
    {
	fprintf(stderr, "\n! %s: Could not perform clustering, exiting...\n", argv[0]);
	exit(EXIT_FAILURE);
    }
    
    print_dslclus(Clus, stdout);
    if (Outfnm!=NULL && Clus->No>=2)   /* struct output requested */
    {
	printf("# Saved to \"%s\"\n\n", Outfnm);
	avg_str(Outfnm, Clus, Structs, Size, Snames, Pdbdescr);
    }
    
    exit(EXIT_SUCCESS);
}

/* ==== FUNCTIONS ==== */

/* ---- Cluster averaging ---- */

/* avg_str(): creates an average structure from the structures
 * contained in the cluster pointed to by Clu. The structures
 * are supplied by Structs (all having size Size), the RMS rotations
 * will be performed with uniform weighting. The average
 * is generated by a bottom-up traversal of the cluster tree.
 * Once it is there, the original structures are superposed
 * onto it and the whole bunch is written to Outfnm in PDB format
 * if Outfnm!=NULL.
 * Snames contains the original file names of the structures, 
 * Pdbdescr is the target PDB record (NULL for smoothed C-alpha
 * comparisons) whose 1st chain is the same for all compared entities.
 */
static void avg_str(const char *Outfnm, const Dslclu_ *Clu, double ***Structs, 
	int Size, char **Snames, const Pdbentry_ *Pdbdescr)
{
    const int REMLEN=60;
    
    double **Avg=NULL, *Arms=NULL;
    double Rmsmin=DBL_MAX, Rmsmax=-9999.9, Rmsavg, Rmssd;
    char **Remarks=NULL;
    char *Name=NULL;
    Pdbentry_ *Pdbout;
    int i, k;
    
    /* don't do anything if Clu is empty */
    if (Clu==NULL) return;
    
    /* leaf */
    if (Clu->No<2) return;
    
    /* get the average structure recursively, store as first chain */
    Avg=cluster_avg(Clu, Structs, Size);
    if (Outfnm!=NULL)
	Pdbout=start_struct(Avg, Size, Pdbdescr);
    
    /* rotate all structs onto the average, save RMSs */
    Arms=(double *) calloc(Clu->No, sizeof(double));
    for (k=0; k<Clu->No; k++)
    {
	i=Clu->Members[k];
	Arms[k]=rotate_vectors(Avg, Structs[i], NULL, Size);
	if (Arms[k]<Rmsmin) Rmsmin=Arms[k];
	if (Arms[k]>Rmsmax) Rmsmax=Arms[k];
	
	/* store as k-th chain */
	if (Outfnm!=NULL)
	    add_struct(Pdbout, Structs[i], Size, Pdbdescr, 'A'+k);
    }
    avg_sd(Arms, Clu->No, &Rmsavg, &Rmssd);
    puts("# RMS deviation from cluster average:");
    for (k=0; k<Clu->No; k++)
    {
	i=Clu->Members[k];
	printf("[%d]: RMS=%.3f A\n", i, Arms[k]);
    }
    printf("\n#    Best   |    Avg    +/-    SD     |   Worst\n  %.3e | %.3e +/- %.3e | %.3e\n\n", 
	    Rmsmin, Rmsavg, Rmssd, Rmsmax);
    
    /* Get the S.D. of each average position and put into average's B-factor field */
    target_sd(Pdbout);
        
    /* finish up PDB output */
    if (Outfnm!=NULL)
    {
	/* set up remarks */
	Remarks=(char **) calloc(Clu->No+6, sizeof(char *));
	Remarks[0]=(char *) calloc(REMLEN, sizeof(char));
	sprintf(Remarks[0], "CLUSTER OF %d STRUCTURE%s, FIRST CHAIN (0) IS THE AVERAGE", 
	    Clu->No, (Clu->No==1)? "": "S");
	for (k=0; k<Clu->No; k++)
	{
	    i=Clu->Members[k];
	    Remarks[k+1]=(char *) calloc(REMLEN, sizeof(char));
	    Name=strrchr(Snames[i], '/');
	    Name=(Name==NULL)? Snames[i]: Name+1;
	    sprintf(Remarks[k+1], "%s %c %.4e A", Name, 'A'+k, Arms[k]);
	}
	
	Remarks[Clu->No+1]=(char *) calloc(REMLEN, sizeof(char));
	sprintf(Remarks[Clu->No+1], "BEST RMS=%.4e A", Rmsmin);
	Remarks[Clu->No+2]=(char *) calloc(REMLEN, sizeof(char));
	sprintf(Remarks[Clu->No+2], "WORST RMS=%.4e A", Rmsmax);
	Remarks[Clu->No+3]=(char *) calloc(REMLEN, sizeof(char));
	sprintf(Remarks[Clu->No+3], "AVERAGE RMS=%.4e +/- %.4e A", Rmsavg, Rmssd);
	Remarks[Clu->No+4]=(char *) calloc(REMLEN, sizeof(char));
	sprintf(Remarks[Clu->No+4], "FIRST CHAIN: B-FACTOR IS S.D. OF DISTANCES FROM AVERAGE");
	Remarks[Clu->No+5]=(char *) calloc(REMLEN, sizeof(char));
	sprintf(Remarks[Clu->No+5], "OTHER CHAINS: B-FACTOR IS DISTANCE FROM AVERAGE");
		
	/* list to PDB file */
	put_pdb(Outfnm, Pdbout, Remarks, Clu->No+6);
	
	/* cleanup */
	free_pdb(Pdbout);
	for (i=0; i<Clu->No+6; i++) free(Remarks[i]);
	free(Remarks);
    }
    
    /* more cleanup */
    for (i=0; i<Size; i++) free(Avg[i]);
    free(Avg);
    free(Arms);
}
/* END of avg_str() */

/* cluster_avg(): if this cluster is a leaf, then a copy of the
 * structure it contains is returned. Otherwise, two recursive
 * calls are made on the subclusters, the averages returned, 
 * then subcluster 2's average is superimposed onto subcluster 1, 
 * the two are averaged in 1, 2's average is deleted and 1 is
 * returned. Logic based on the property of the Dslclu_ trees that there
 * are either leaves or 2-nodes and the leftmost branches are
 * the deepest.
 */
static double **cluster_avg(const Dslclu_ *Clu, double ***Structs,
	int Size)
{
    double **Substr1=NULL, **Substr2=NULL;
    double Rms;
    int i, j, N1, N2;
    
    /* leaf */
    if (Clu->Sub1==NULL && Clu->Sub2==NULL)
    {
	/* make copy of the only struct it contains */
	Substr1=(double **) calloc(Size, sizeof(double *));
	i=Clu->Members[0];
	for (j=0; j<Size; j++)
	{
	    Substr1[j]=(double *) calloc(3, sizeof(double));	/* one vector */
	    memcpy(Substr1[j], Structs[i][j], 3*sizeof(double));
	}
	return(Substr1);
    }
    
    /* not leaf, has two sub-nodes. Get the (average)
     * structures from both, rotate the 2nd onto the 1st, 
     * weighted average into the 1st, remove 2nd and return 1st
     */
    Substr1=cluster_avg(Clu->Sub1, Structs, Size);
    Substr2=cluster_avg(Clu->Sub2, Structs, Size);
    
    Rms=rotate_vectors(Substr1, Substr2, NULL, Size);
    
    N1=Clu->Sub1->No; N2=Clu->Sub2->No;
    for (j=0; j<Size; j++)
    {
	for (i=0; i<3; i++) /* average vectors in the 1st set */
	{
	    Substr1[j][i]*=N1;
	    Substr1[j][i]+=N2*Substr2[j][i];
	    Substr1[j][i]/=(N1+N2);
	}
	free(Substr2[j]);
    }
    free(Substr2);
    return(Substr1);
}
/* END of cluster_avg() */

/* avg_sd: returns the average and SD (uncorrected, N-weighted) of a
 * data vector Data[] that is Datno long. Avg, Sd are not modified if
 * Datno<=0. Sd==1.0 if Datno=1.
 */
static void avg_sd(double Data[], int Datno, double *Avg, double *Sd)
{
    double Sx, Sx2;
    int d;

    *Avg= *Sd=0.0;
    if (Datno<=0) return;

    Sx=Sx2=0.0;
    for (d=0; d<Datno; d++)
    {
	Sx+=Data[d];
	Sx2+=(Data[d]*Data[d]);
    }
    *Avg=Sx/Datno;
    *Sd=(Datno==1)? 0.0: sqrt(fabs((Sx2-Sx*Sx/Datno)/Datno));
}
/* END of avg_sd */

/* print_help: lists the available options to stderr. */
static void print_help(const char *Progname)
{
    fprintf(stderr, "Usage: %s %s PDB_files... \n", Progname, opt_helpstr());
    fputs("\t the PDB_file(s) will be aligned to each other and clustered\nOptions:-\n", stderr);
    fputs("\t-o outfile: save aligned structures to \"outfile\" in PDB format\n", stderr);
    fputs("\t-a: superimpose all atoms (default C-alpha only)\n", stderr);
    fputs("\t-s: smooth the C-alpha trace (default off, ignored with \"-a\" option)\n", stderr);
    fprintf(stderr, "\t-w <int>, default=%d: smoothing window length\n",
	WINLEN);
    fprintf(stderr, "\t-c <int>, default=%d: no. of smoothing cycles\n",
	SMCYC);
}
/* END of print_help() */

/* ==== END OF PROGRAM clumsy.c ==== */
