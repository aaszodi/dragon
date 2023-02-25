/* ==== PROGRAM rank.c ==== */

/* For ranking DRAGON output PDB files. */

/* ANSI C, 21. June 1998. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* ---- MODULES ---- */

#include "cmdopt.h"

/* ---- TYPEDEFS ---- */

/* Sco_: stores the name of a DRAGON output file in PDB format,
 * the scores (bond, nonbond, external restraint, secstr)
 * and the rankings for each score.
 */
typedef struct
{
	char *Name;	/* PDB file name */
	float Bn, Nb, Rs, Sc;	/* score values */
	int Brank, Nrank, Rrank, Srank;	/* individual score ranks */
} Sco_;

/* Scoflag_: flags for choosing scores to be taken into account
 * when performing composite scoring. The flags may be ORed together.
 */
typedef enum {BN=1, NB=2, RS=4, SC=8} Scoflag_;

/* ---- GLOBAL VARIABLES ---- */

int Scoflags=0;	/* holds the score flags for composite ranking */

/* ---- PROTOTYPES ---- */

int cmp_bn(const void* p1, const void *p2);
int cmp_nb(const void* p1, const void *p2);
int cmp_rs(const void* p1, const void *p2);
int cmp_sc(const void* p1, const void *p2);
int cmp_composite(const void* p1, const void *p2);

/* ==== MAIN ==== */

/* The program takes one or more PDB filenames as arguments.
 * These are supposed to have been created by DRAGON Version 4.16
 * or above and must contain score lines among their REMARK cards.
 */
int main(int argc, char *argv[])
{
	#define LINELEN 120

	Sco_ *Scos=NULL;
	char Line[LINELEN], *Sp;
	FILE *Inf=NULL;
	int Fileno,i,Firstfile,Rdnb,Rdbn,Rdrs, Rdsc;
	
	/* get options and file parameters */
	parse_optstr("bnrs");
	Firstfile=get_options(argc, argv);
	if (Firstfile<0 || argc-Firstfile<1)
	{
		fprintf(stderr,"\n! Usage: %s %s DRAGON_PDB_file(s)\n",
			argv[0], opt_helpstr());
		fputs("\t-b: sort on bond score\n", stderr);
		fputs("\t-n: sort on non-bond score\n", stderr);
		fputs("\t-r: sort on restraint score\n", stderr);
		fputs("\t-s: sort on secondary structure score\n", stderr);
		fputs("\tFlags may be combined, default: -bnrs\n", stderr);
		exit(1);
	}

	/* set up the score flags */
	if (optval_bool('b')) Scoflags|=BN;
	if (optval_bool('n')) Scoflags|=NB;
	if (optval_bool('r')) Scoflags|=RS;
	if (optval_bool('s')) Scoflags|=SC;
	if (!Scoflags) Scoflags=(BN|NB|RS|SC); /* default all */
	
	/* get the score values from the file(s) */
	Scos=(Sco_*)calloc(argc-1,sizeof(Sco_));
	for (i=Firstfile, Fileno=0; i<argc; i++)
	{
		Inf=fopen(argv[i],"r");
		if (Inf==NULL)
		{
			fprintf(stderr,"\n? %s: Cannot open \"%s\",skipped\n",argv[0],argv[i]);
			continue;
		}
		Rdbn=Rdnb=Rdrs=Rdsc=0;
		while(NULL!=fgets(Line,LINELEN,Inf))
		{
			if (!Rdbn && NULL!=(Sp=strstr(Line,"BOND SCORE:")))
			{
				Rdbn=sscanf(Sp+11,"%f",&(Scos[Fileno].Bn));
				continue;
			}
			if (!Rdnb && NULL!=(Sp=strstr(Line,"BUMP SCORE:")))
			{
				Rdnb=sscanf(Sp+11,"%e",&(Scos[Fileno].Nb));
				continue;
			}
			if (!Rdrs && NULL!=(Sp=strstr(Line,"RESTRAINT SCORE:")))
			{
				Rdrs=sscanf(Sp+16,"%e",&(Scos[Fileno].Rs));
				continue;
			}
			if (!Rdsc && NULL!=(Sp=strstr(Line,"SECONDARY STRUCTURE SCORE:")))
			{
				Rdsc=sscanf(Sp+26,"%e",&(Scos[Fileno].Sc));
				continue;
			}
		}
		fclose(Inf);
		
		if (!Rdbn || !Rdnb || !Rdrs || !Rdsc)
		{
			fprintf(stderr,"\n? %s: Scores(s) missing from \"%s\",skipped\n",argv[0],argv[i]);
			continue;
		}

		/* looks OK, store filename */
		Scos[Fileno].Name=(char *)calloc(strlen(argv[i])+1,sizeof(char));
		Scos[Fileno].Name[strlen(argv[i])]='\0';
		strcpy(Scos[Fileno].Name,argv[i]);
		Fileno++;
	}
	
	if(!Fileno)
	{
		fprintf(stderr,"\n? %s: No valid files\n",argv[0]);
		exit(2);
	}

	/* rank on each score */
	if (Fileno>=2)
	{
	    if (Scoflags & BN)
	    {
		qsort(Scos, Fileno, sizeof(Sco_),cmp_bn);
		for (i=0; i<Fileno; i++) 
		    Scos[i].Brank=(i && Scos[i].Bn==Scos[i-1].Bn)? Scos[i-1].Brank: i+1;
	    }
	    if (Scoflags & NB)
	    {
		qsort(Scos, Fileno, sizeof(Sco_),cmp_nb);
		for (i=0; i<Fileno; i++)
		    Scos[i].Nrank=(i && Scos[i].Nb==Scos[i-1].Nb)? Scos[i-1].Nrank: i+1;
	    }
	    if (Scoflags & RS)
	    {
		qsort(Scos, Fileno, sizeof(Sco_),cmp_rs);
		for (i=0; i<Fileno; i++)
		    Scos[i].Rrank=(i && Scos[i].Rs==Scos[i-1].Rs)? Scos[i-1].Rrank: i+1;
	    }
	    if (Scoflags & SC)
	    {
		qsort(Scos, Fileno, sizeof(Sco_),cmp_sc);
		for (i=0; i<Fileno; i++)
		    Scos[i].Srank=(i && Scos[i].Sc==Scos[i-1].Sc)? Scos[i-1].Srank: i+1;
	    }
	    
	    /* sort on composite rank if required */
	    qsort(Scos,Fileno,sizeof(Sco_),cmp_composite);
	}
	
	/* listing */
	for (i=0; i<Fileno; i++)
	{
		printf("%d %s", i+1, Scos[i].Name);
		if (Scoflags & BN)
		    printf(" Bn=%9.3e (%d)", Scos[i].Bn,Scos[i].Brank);
		if (Scoflags & NB)
		    printf(" Nb=%9.3e (%d)", Scos[i].Nb,Scos[i].Nrank);
		if (Scoflags & RS)
		    printf(" Rs=%9.3e (%d)", Scos[i].Rs,Scos[i].Rrank);
		if (Scoflags & SC)
		    printf(" Sc=%9.3e (%d)", Scos[i].Sc,Scos[i].Srank);
		putchar('\n');
	}
	exit(0);
}

/* ==== FUNCTIONS ==== */

int cmp_bn(const void* p1, const void *p2)
{
	Sco_ *s1=(Sco_*)p1, *s2=(Sco_*)p2;
	if (s1->Bn==s2->Bn) return(0);
	return((s1->Bn<s2->Bn)? -1: 1);
}
int cmp_nb(const void* p1, const void *p2)
{
	Sco_ *s1=(Sco_*)p1, *s2=(Sco_*)p2;
	if (s1->Nb==s2->Nb) return(0);
	return((s1->Nb<s2->Nb)? -1: 1);
}
int cmp_rs(const void* p1, const void *p2)
{
	Sco_ *s1=(Sco_*)p1, *s2=(Sco_*)p2;
	if (s1->Rs==s2->Rs) return(0);
	return((s1->Rs<s2->Rs)? -1: 1);
}
int cmp_sc(const void* p1, const void *p2)
{
	Sco_ *s1=(Sco_*)p1, *s2=(Sco_*)p2;
	if (s1->Sc==s2->Sc) return(0);
	return((s1->Sc<s2->Sc)? -1: 1);
}
int cmp_composite(const void* p1, const void *p2)
{
	extern int Scoflags;	/* says which score to use: global */
	Sco_ *s1=(Sco_*)p1, *s2=(Sco_*)p2;
	int c1=0, c2=0;
	
	if (Scoflags & BN)
	{
	    c1+=s1->Brank; c2+=s2->Brank;
	}
	if (Scoflags & NB)
	{
	    c1+=s1->Nrank; c2+=s2->Nrank;
	}
	if (Scoflags & RS)
	{
	    c1+=s1->Rrank; c2+=s2->Rrank;
	}
	if (Scoflags & SC)
	{
	    c1+=s1->Srank; c2+=s2->Srank;
	}
	return(c1-c2);
}

/* ==== END OF PROGRAM rank.c ==== */

