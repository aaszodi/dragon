/* ==== FUNCTIONS rotpdb.c ==== */

/* Contains routines for PDB I/O and for optimal superposition
 * of PDB-derived C-alpha chains.
 */

/* ANSI C, 30. Apr. 1998. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <time.h>

/* ---- MODULE HEADER ---- */

#include "rotpdb.h"

/* NOTE: SGI provides single-precision floating point functions
 * such as sqrtf() etc. Some machines (SUNs in particular) don't
 * know about this: Define NO_MATHFLOATFUNC on the command line
 */
#ifdef NO_MATHFLOATFUNC
#define sqrtf sqrt
#endif

/* ==== FUNCTIONS ==== */

/* ---- Rotation ---- */

/* get_vectors: transfers the coordinates of the
 * first chain in the PDB file Pdbfnm to an array of 3D
 * vectors suitable for processing by the Bestrot routines.
 * If *Newentry==NULL, then it will be set to point to
 * the PDB entry record constructed from the contents of
 * the PDB file (used for reading the first entry). If *Newentry!=NULL,
 * then the first chain in it is compared to the first chain
 * just read in and NULL is returned if they don't match.
 * Ignored if Newentry==NULL.
 * If Allatoms!=0, then all atoms are read, otherwise C-alphas only.
 * The vectors are centred on the centroid.
 * The no. of vectors is returned in Size, the return value
 * is the [Size][3] array of the vectors or NULL if the file
 * could not be read or was not a protein.....etc.
 * Also returns the sequence of the chain in 1-letter codes in Seq.
 * Also creates and returns a weight vector in *Wgt if Wgt!=NULL;
 * the weights will be "reciprocal" to the B-factors in the
 * structure (if meaningful B-factors are specified) or
 * uniform otherwise.
 */
double **get_vectors(const char *Pdbfnm, Pdbentry_ **Newentry, 
	char Allatoms, double **Wgt, int *Size)
{
    Pdbentry_ *Entry=NULL;
    double **Vectors=NULL;
    double *Ctr=NULL, *W=NULL;
    double Bfact, Bmin=HUGE_VAL, Bmax=-HUGE_VAL;
    Chain_ *Chain;
    int i, j, c, s;
    
    /* read the PDB entry */
    Entry=get_pdb(Pdbfnm, Allatoms? ALLATOMS: CALPHA, STRICT);
    if (Entry==NULL || Entry->Chainno<=0)
    {
	free_pdb(Entry); *Size=0; return(NULL);
    }
    
    /* if *Newentry!=NULL, then compare the first chains
     * between the two (by comparing the sequences and number of
     * atoms): return an error on mismatch
     */
    if (Newentry!=NULL && *Newentry!=NULL)
    {
	if (Entry->Chains->Atomno!=(*Newentry)->Chains->Atomno ||
	    strcmp(Entry->Chains->Seq, (*Newentry)->Chains->Seq))
	{
	    fputs("\n? get_vectors(): Current chain does not match target\n", stderr);
	    free_pdb(Entry); *Size=0; return(NULL);
	}
    }
    
    /* process the first chain only (monomers...) */
    Chain=Entry->Chains;
    s=Allatoms? Chain->Atomno: Chain->Aano;
    Vectors=(double **) calloc(s, sizeof(double *));
    
    /* copy the coordinates */
    for (i=0; i<s; i++)
    {
	Vectors[i]=(double *) calloc(3, sizeof(double));
	Vectors[i][0]=Chain->Atoms[i].X;
	Vectors[i][1]=Chain->Atoms[i].Y;
	Vectors[i][2]=Chain->Atoms[i].Z;
    }
    
    /* construct weight vector if needed */
    if (Wgt!=NULL)
    {
	/* get minimal and maximal B-factor */
	for (i=0; i<s; i++)
	{
	    Bfact=Chain->Atoms[i].Bfact;
	    if (Bfact<Bmin) Bmin=Bfact;
	    if (Bfact>Bmax) Bmax=Bfact;
	}
	
	W=(double *) calloc(s, sizeof(double));
	if (Bmax==Bmin || Bmin>Bmax)	/* make uniform wgts */
	    for (i=0; i<s; i++) W[i]=1.0;
	else
	{	/* inverse: Bmin->1, Bmax->0 */
	    Bmin=Bmax-Bmin;
	    for (i=0; i<s; i++)
		W[i]=(Bmax-Chain->Atoms[i].Bfact)/Bmin;
	}
	*Wgt=W;
    }
    
    /* shift to centroid */
    Ctr=center_vectors(Vectors, NULL, (Wgt!=NULL)? W: NULL, s);
    free(Ctr);
    
    
    *Size=s;
    
    /* If *Newentry was NULL, this means that it should be
     * set to point to where Entry pointed. Otherwise, 
     * remove Entry
     */
    if (Newentry!=NULL && *Newentry==NULL) *Newentry=Entry;
    else free_pdb(Entry);
    return(Vectors);
}
/* END of get_vectors */

/* rotate_vectors: performs the McLachlan rotation on the vector sets
 * Target and Struct (both [Size][3]) so that Struct will be
 * rotated to Target. Weights are supplied in W[].
 * Return value: the RMS difference.
 */
double rotate_vectors(double **Target, double **Struct, 
	const double W[], int Size)
{
    Sqmat_ Transform;
    int i, j, k;
    double Rms, Temp[3];
    
    Transform=alloc_sqmat(3);   /* the rotation matrix */
    
    /* get the RMS value and the transform matrix */
    Rms=best_rot(Struct, Target, W, Size, Transform);
    
    /* perform the rotation in place */
    for (i=0; i<Size; i++)
    {
	for (j=0; j<3; j++)
	{
	    Temp[j]=0.0;
	    for (k=0; k<3; k++)
		Temp[j]+=Transform[j][k]*Struct[i][k];
	}
	for (j=0; j<3; j++) Struct[i][j]=Temp[j];
    }
    
    free_matrix(Transform);    
    return(Rms);
}
/* END of rotate_vectors */

/* ---- Smoothing ---- */

/* smooth_chains: do a moving-average smoothing of the
 * coordinate vectors stored in Struct (Size x 3) Cycno times
 * with a Wlen-wide window. The smoothed results will be kept in
 * Struct which will be truncated appropriately on exit.
 * Return value: a pointer to the realloc()-d Struct whose
 * new size is adjusted in Size.
 */
double **smooth_chains(double **Struct, int *Size, int Wlen, int Cycno)
{
    int c, i, j, w, Len;
    double **Newstruct;
    
    /* do smoothing Cycno times */
    for (Len=*Size, c=0; c<Cycno; c++)
    {
	Len-=(Wlen-1);	/* averaged length will shrink by Wlen-1 */
	/* average Wlen neighbours in place in the first element
	 * of the window (indexed by i)
	 */
	for (i=0; i<Len; i++)
	{
	    for (w=i+1; w<i+Wlen; w++)	/* averaging in window */
		for (j=0; j<3; j++)
		    Struct[i][j]+=Struct[w][j];
	    for (j=0; j<3; j++)
		Struct[i][j]/=Wlen; /* result kept in i-th */
	}   /* for i */
    }	/* for c */
    
    /* size adjustment: chop off tail */
    for (i=Len; i<*Size; i++) free(Struct[i]);
    Newstruct=(double **) realloc(Struct, Len*sizeof(double *));
    *Size=Len;
    return(Newstruct);
}
/* END of smooth_chains */

/* smooth_wgt: smooths the weight vector Wgt. The algorithm and
 * the parameters used are the same as in smooth_chains().
 * Note that Size is not modified and Wgt is not-reallocated.
 */
void smooth_wgt(double Wgt[], int Size, int Wlen, int Cycno)
{
    int c, i, j, w, Len;
    
    /* do smoothing Cycno times */
    for (Len=Size, c=0; c<Cycno; c++)
    {
	Len-=(Wlen-1);	/* averaged length will shrink by Wlen-1 */
	/* average Wlen neighbours in place in the first element
	 * of the window (indexed by i)
	 */
	for (i=0; i<Len; i++)
	{
	    for (w=i+1; w<i+Wlen; w++)	/* averaging in window */
		Wgt[i]+=Wgt[w];
	    Wgt[i]/=Wlen;
	}   /* for i */
    }	/* for c */
}
/* END of smooth_wgt */

/* ---- PDB file construction ---- */

/* start_struct: creates a PDB entry structure (cf. "pdbprot.h")
 * from the Target vectors interpreted as coordinates
 * for the first chain in Pdbtarg. If it is NULL,
 * it means that smoothing on C-alphas was
 * done and no sequences are saved. (Residue IDs will be "XXX".)
 * Return value: ptr to the entry or NULL if Targetsize<=0.
 */
Pdbentry_ *start_struct(double **Target, int Targetsize, 
	const Pdbentry_ *Pdbtarg)
{
    Pdbentry_ *Entry=NULL;
    Chain_ *Chain=NULL;
    int i;
    time_t Now;
    
    if (Targetsize<=0) return(NULL);
    
    /* unused administrative fields */
    Entry=(Pdbentry_ *) malloc(sizeof(Pdbentry_));
    strcpy(Entry->Header, "ALIGNED STRUCTURES");
    Now=time(NULL);  /* get current date */
    strftime(Entry->Date, 10, "%d-%b-%y", localtime(&Now));
    strcpy(Entry->Pdbcode, "0ROT");
    Entry->Compound=(char *) calloc(19, sizeof(char));
    strcpy(Entry->Compound, "POLYPEPTIDE CHAINS");
    Entry->Source=(char *) calloc(13, sizeof(char));
    strcpy(Entry->Source, "SIMULATIONS");
    strcpy(Entry->Expdta, "RIGID-BODY ROTATION");
    Entry->Resol=-1.0;	/* no resolution */
    Entry->Hbonds=NULL; Entry->Ssbs=NULL;
    Entry->Hbno=Entry->Ssbno=0;
    
    /* create first chain that stores the Target coords */
    Chain=(Chain_ *) malloc(sizeof(Chain_));
    Chain->Secs=NULL; Chain->Hbonds=NULL; Chain->Ssbs=NULL;
    Chain->Secsno=Chain->Hbno=Chain->Ssbno=0;
    Chain->Chid='0';	/* so that RasMol could access this chain */
    Chain->Type=(Pdbtarg==NULL)? 'A': Pdbtarg->Chains->Type;	/* C-alpha or full protein trace */
    Chain->Atomno=Targetsize;
    Chain->Aano=(Pdbtarg==NULL)? Targetsize: Pdbtarg->Chains->Aano;
    Chain->Atoms=(Atom_ *) calloc(Targetsize, sizeof(Atom_ ));
    
    /* Create the sequence array and copy it if Pdbtarg is known */
    if (Pdbtarg!=NULL)
    {
	Chain->Seq=(char *) calloc(strlen(Pdbtarg->Chains->Seq)+1, sizeof(char));
	strcpy(Chain->Seq, Pdbtarg->Chains->Seq);
    }
    else
	Chain->Seq=(char *) calloc(Targetsize+1, sizeof(char));
    
    /* copy the coordinates */
    for (i=0; i<Targetsize; i++)
    {
	if (Pdbtarg!=NULL)	/* take from target */
	    Chain->Atoms[i]=Pdbtarg->Chains->Atoms[i];
	else
	{   /* smoothed C-alpha */
	    Chain->Atoms[i].Resno=i+1;
	    strcpy(Chain->Atoms[i].Id, "CA");   /* all C-alpha */
	    Chain->Atoms[i].Alt=Chain->Atoms[i].Rid=' ';
	    Chain->Atoms[i].Aa=Chain->Seq[i]='X';   /* unknown because smoothed */
	}
	Chain->Atoms[i].Atno=i+1;  /* numbering starts from 1 */
	Chain->Atoms[i].X=Target[i][0];
	Chain->Atoms[i].Y=Target[i][1];
	Chain->Atoms[i].Z=Target[i][2];
	Chain->Atoms[i].Occu=1.0;
	Chain->Atoms[i].Bfact=0.0;
    }
    
    /* finish up */
    Entry->Chains=Chain;
    Entry->Chainno=1;
    return(Entry);
}
/* END of start_struct */

/* add_struct: adds a new chain to Entry which contains the 
 * vectors in Struct,  interpreted as coordinates belonging
 * to the first chain in Pdbtarg if it is not NULL.
 * The new chain is Targetsize long, with Chainid as a
 * character ID. It is assumed that there's a 1:1 correspondence
 * between the vectors in Struct and the C-alpha coordinates
 * in Target (see above) stored in Entry->Chain[0].
 * Structseq==NULL means that smoothing has been done and no
 * sequence information will be saved.
 * The B-factor field will contain the distance (in angstroms)
 * between the corresponding atoms in Target and Struct.
 * This is for the Quanta "4th-parameter" colouring option.
 */
void add_struct(Pdbentry_ *Entry, double **Struct, int Targetsize, 
	const Pdbentry_ *Pdbtarg, char Chainid)
{
    Chain_ *Chain;
    static int Lastatno=0;
    int i;
    
    if (Lastatno==0) Lastatno=Targetsize;
    
    /* create the chain that stores the Struct coords */
    (Entry->Chainno)++;
    Entry->Chains=(Chain_ *) realloc(Entry->Chains, (Entry->Chainno)*sizeof(Chain_));
    Chain=Entry->Chains+(Entry->Chainno)-1;
    Chain->Secs=NULL; Chain->Hbonds=NULL; Chain->Ssbs=NULL;
    Chain->Secsno=Chain->Hbno=Chain->Ssbno=0;
    Chain->Chid=Chainid;    /* set the chain ID */
    Chain->Type=(Pdbtarg==NULL)? 'A': Pdbtarg->Chains->Type;	/* C-alpha or full protein trace */
    Chain->Atomno=Targetsize;
    Chain->Aano=(Pdbtarg==NULL)? Targetsize: Pdbtarg->Chains->Aano;
    Chain->Atoms=(Atom_ *) calloc(Targetsize, sizeof(Atom_ ));
    
    /* Create the sequence array and copy it if Pdbtarg is known */
    if (Pdbtarg!=NULL)
    {
	Chain->Seq=(char *) calloc(strlen(Pdbtarg->Chains->Seq)+1, sizeof(char));
	strcpy(Chain->Seq, Pdbtarg->Chains->Seq);
    }
    else
	Chain->Seq=(char *) calloc(Targetsize+1, sizeof(char));
    
    
    /* copy the coordinates */
    for (i=0; i<Targetsize; i++)
    {
	if (Pdbtarg!=NULL)	/* take from target */
	    Chain->Atoms[i]=Pdbtarg->Chains->Atoms[i];
	else
	{   /* smoothed C-alpha */
	    Chain->Atoms[i].Resno=i+1;
	    strcpy(Chain->Atoms[i].Id, "CA");   /* all C-alpha */
	    Chain->Atoms[i].Alt=Chain->Atoms[i].Rid=' ';
	    Chain->Atoms[i].Aa=Chain->Seq[i]='X';
	}
	Chain->Atoms[i].Atno=Lastatno+i+1;  /* contig numbering */
	Chain->Atoms[i].X=Struct[i][0];
	Chain->Atoms[i].Y=Struct[i][1];
	Chain->Atoms[i].Z=Struct[i][2];
	Chain->Atoms[i].Occu=1.0;
	Chain->Atoms[i].Bfact=atom_dist(Entry->Chains->Atoms+i, Chain->Atoms+i);
    }
    Lastatno+=Targetsize;
}
/* END of add_struct */

/* target_sd(): calculate the standard deviation of distances 
 * from the target for each atom and put these into the corresponding
 * B-factor entries of the target chain. We assume that all chains
 * have already been added to Entry before the call and the structure atom
 * B-factor fields contain the distance from the target atom. The B-factor
 * values are set to 0.0 if there were less than 2 structures in
 * addition to the target.
 */
void target_sd(Pdbentry_ *Entry)
{
    int i, k;
    float Sum, D;
    Atom_ *Targatoms;
    
    
    if (Entry->Chainno<2) return;   /* too few... */
    if (Entry->Chainno<3)   /* only one struct */
    {
	for (k=0; k<Entry->Chains->Atomno; k++)
	    Entry->Chains->Atoms[k].Bfact=0.0;
	return;
    }
    
    /* general case: loop over structs and atoms */
    Targatoms=Entry->Chains->Atoms;
    for (k=0; k<Entry->Chains->Atomno; k++)
    {
	Sum=0.0;
	for (i=1; i<Entry->Chainno; i++)
	{
	    D=Entry->Chains[i].Atoms[k].Bfact;
	    Sum+=D*D;
	}
	Targatoms[k].Bfact=sqrtf(Sum/(Entry->Chainno-2));
    }
}
/* END of target_sd() */

/* ==== END OF FUNCTIONS rotpdb.c ==== */
