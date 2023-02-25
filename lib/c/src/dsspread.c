/* ==== FUNCTIONS dsspread.c ==== */

/* DSSP reader routine. */

/* ANSI C, IRIX 5.3, 30. Apr. 1996. Andris */

/* ---- MODULE HEADER ---- */

#include "dsspread.h"

/* ==== FUNCTIONS ==== */

/* dssp_read: reads the text file Dsspfnm produced
 * by DSSP, Version Oct. 1988 (ref: W. Kabsch, C. Sander, Biopolymers
 * 22:2577-2637 (1983)). Returns an array of size Size, containing
 * chain ID, residue no, 1-letter AA code, Kabsch/Sander secondary
 * structure code and solvent accessibility for each residue
 * or NULL on error. Sets the number of chains to Chainno.
 */
Dssprec_ *dssp_read(const char *Dsspfnm, unsigned int *Size, unsigned int *Chainno)
{
    #define LINELEN 132
    
    FILE *Dssp=NULL;
    char Line[LINELEN],Dummy[18];  /* read-in buffers */
    Dssprec_ *Entry=NULL, *Cur;            /* temp pointer */
    unsigned int i, Entryj, Nres;
    char c;
    
    /* init reading, skip lines */
    Nres=*Chainno=0;
    
    /* opening ceremony */
    if (NULL==(Dssp=fopen(Dsspfnm, "r")))
    {
	fprintf(stderr, "\n? dssp_read(%s): Cannot open\n", Dsspfnm);
	return(NULL);
    }

    /* read total no. of residues, no. of chains */
    do
    {
	if(NULL==fgets(Line,LINELEN,Dssp))
	{
	    fprintf(stderr, "\n? dssp_read(%s): Cannot read\n", Dsspfnm);
	    fclose(Dssp); return(NULL);
	}
	if (NULL!=strstr(Line, "TOTAL NUMBER OF RESIDUES"))
	    break;
    } while(1); 
    if (2>sscanf(Line,"%d %d",&Nres,Chainno))
    {
	fprintf(stderr, "\n? dssp_read(%s): Cannot parse:\n%s\n", Dsspfnm, Line);
	fclose(Dssp); return(NULL);
    }
    
    /* skip all info until the last text line is encountered */
    do 
    {
	if(NULL==fgets(Line,LINELEN,Dssp))
	{
	    fprintf(stderr, "\n? dssp_read(%s): Cannot read\n", Dsspfnm);
	    fclose(Dssp); return(NULL);
	}
	if(NULL!=strstr(Line, "#  RESIDUE AA STRUCTURE"))
	    break;
    }
    while (1);
        
    /* adjust storage, allocate */
    Nres+=*Chainno-1;
    Entry= (Dssprec_ *) calloc(Nres,sizeof(Dssprec_));

    while (NULL!=fgets(Line,LINELEN,Dssp)) 
    {
	sscanf(Line,"%d",&Entryj);  /* reads pos */
	Entryj--;   /* indexes res position now */
        Cur=Entry+Entryj;  /* current pos in array */
	if (Line[13]=='!')
	{                   /* chain break */
	    Cur->Resno=0; Cur->Res='!';
	    Cur->Disulf=Cur->Chain=Cur->Secstruct=Cur->Turns3=Cur->Turns4=Cur->Turns5=' ';
	    Cur->Bend=Cur->Chir=Cur->Bridge1=Cur->Bridge2=Cur->Sheet=' ';
	    Cur->Access=0;
	    Cur->Phi=Cur->Psi=0.0;
	    Cur->Ca[0]=Cur->Ca[1]=Cur->Ca[2]=0.0;
            continue;              /* nothing else to be done */
	}

	/* normal entry */
	sscanf(Line, "%*d%d", &(Cur->Resno));
	Cur->Chain=Line[11];	/* chain ID or ' ' if single chain */
	Cur->Res=Line[13];	/* AA code */
	
	/* transfer the funny char descriptors */
	Cur->Secstruct=Line[16];
	Cur->Turns3=Line[18]; Cur->Turns4=Line[19]; Cur->Turns5=Line[20];
	Cur->Bend=Line[21]; Cur->Chir=Line[22];
	Cur->Bridge1=Line[23]; Cur->Bridge2=Line[24];
	sscanf(Line+25, "%d%d", &(Cur->Beta1), &(Cur->Beta2));
	Cur->Sheet=Line[33];
	
	/* read accessibility,phi,psi,and C-alpha coordinates */
	sscanf(Line+34,
	  "%d%d,%lf%d,%lf%d,%lf%d,%lf%lf%lf%lf%lf%lf%lf%lf%lf",
	  &(Cur->Access),
	  &(Cur->Nho[0].Offs), &(Cur->Nho[0].En), 
	  &(Cur->Ohn[0].Offs), &(Cur->Ohn[0].En), 
	  &(Cur->Nho[1].Offs), &(Cur->Nho[1].En), 
	  &(Cur->Ohn[1].Offs), &(Cur->Ohn[1].En), 
	  &(Cur->Tco), &(Cur->Kappa), &(Cur->Alpha), 
	  &(Cur->Phi), &(Cur->Psi), &(Cur->Ca[0]),&(Cur->Ca[1]),&(Cur->Ca[2]));
	
	/* Rename half-cystines. These are coded as 'a'..'z' */
	if (Cur->Res>='a' && Cur->Res<='z')
	{
	    Cur->Disulf=Cur->Res; Cur->Res='C';
	}
	else Cur->Disulf=' ';
    }        /* while */
    
    fclose(Dssp);
    *Size=Nres; return(Entry);
    #undef LINELEN    
}     
/* END of dssp_read() */

/* dssp_cadist(): returns the CA:CA distance between two DSSP records
 * pointed to by Dp1 and Dp2. Checks if they're NULL.
 */
double dssp_cadist(const Dssprec_ *Dp1, const Dssprec_ *Dp2)
{
    register double Dist=0.0, D;
    register unsigned int i;
    
    if (Dp1==NULL || Dp2==NULL)
    {
	fprintf(stderr, "\n? dssp_cadist(%x, %x): NULL argument\n", Dp1, Dp2);
	return(0.0);
    }
    
    for (i=0; i<3; i++)
    {
	D=Dp1->Ca[i]-Dp2->Ca[i]; 
	Dist+=D*D;
    }
    return(sqrt(Dist));
}
/* END of dssp_cadist() */

/* ==== END OF FUNCTIONS dsspread.c ==== */
