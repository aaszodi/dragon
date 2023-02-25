// ==== PROJECT DRAGON: FUNCTIONS Output.c++ ====

/* Lists the simulation results to a file in PDB format. 
 * Converts the result first to the C structure Pdbentry_
 * and then uses the C module "pdbprot" for the actual output.
 */

// SGI C++ 7.1, IRIX 6.2, 23. Oct. 1997. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdio.h>
#include <iostream.h>
#include <iomanip.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <errno.h>
#include <time.h>

// ---- UTILITY HEADERS ----

#include "List1.h"
#include "pdbprot.h" 

// ---- MODULE HEADERS ----

#include "Output.h"
#include "Pieces.h"
#include "Fakebeta.h"
#include "version.h"

// ---- PROTOTYPES ----

static int prepare_basename(String_& Basename);
static int mkdir_p(const char *Path);

static Pdbentry_ *make_pdbentry(const Points_& Xyz, const Polymer_& Model, 
	const Pieces_& Pieces);
static void make_secs(Chain_ *Chain, const Pieces_& Pieces);
static void make_atoms(Chain_ *Chain, const Points_& Xyz, const Polymer_& Model);

// ==== FUNCTIONS ====

/* make_outname(): constructs an output filename of the form
 * "Basename_X.Ext" where X is the run number Rcyc. 
 * If Basename contains a directory path, then this path
 * will be created if necessary (permissions permitting).
 * If path creation fails, then the dirpath is thrown away
 * from Basename (will be used in the current directory).
 */
void make_outname(String_& Basename, int Rcyc, const String_& Ext)
{
    static String_ Numstr;
    
    prepare_basename(Basename);
    Basename+="_"; 
    Numstr.long_str(Rcyc);
    Basename+=Numstr;
    Basename+=".";
    Basename+=Ext;
}
// END of make_outname()

/* prepare_basename(): if Basename consists of a directory
 * path and a filename, then the path is created if necessary.
 * If it cannot be created, then the path will be deleted
 * from Basename.
 */
static int prepare_basename(String_& Basename)
{
    int Lastslash=Basename.strrchr('/');
    if (Lastslash<0) return(0);	    // no dirpath, OK
    
    Basename[Lastslash]='\0';	// dirpath only
    int Retval=mkdir_p(Basename);   // create dirpath
    
    if (Retval)	    // something went wrong, throw away dirpath
    {
	String_ Fileonly(((const char *)Basename)+Lastslash+1);
	Basename=Fileonly;
    }
    else Basename[Lastslash]='/';   // join dirpath and filename again
    return(Retval);
}

/* mkdir_p(): creates the path Path with rwxrwxrwx permissions,
 * spiced by the current umask of the process. Does the same
 * job as the 'mkdir -p' UNIX command.
 * Return value: 0 if OK, or the UNIX error number (errno)
 * if something went wrong.
 * NOTE: C strings are used, but C++ alloc and output.
 */
static int mkdir_p(const char *Path)
{
    static const mode_t DIR_MODE=S_IRWXU|S_IRWXG|S_IRWXO;   /* mode 777 */
    char *P=NULL, *Cur=NULL, *C=NULL;
    struct stat Statbuf;
    int Plen;
    
    /* don't do anything if Path is empty */
    if (Path==NULL || !(Plen=strlen(Path))) return(0);
    
    /* make a working copy of Path, append a slash
     * at the end if there is none
     */
    if (Path[Plen-1]!='/') Plen++;
    P=new char [Plen+1];
    strcpy(P, Path);
    P[Plen-1]='/';
    
    for(Cur=(*P=='/')? P+1:P;	/* skip leading slash */
	NULL!=(C=strchr(Cur, '/')); Cur=C+1)
    {
	*C='\0';    /* truncating P up to C */
	errno=0;    /* reset error code */
	
	/* The idea is: if P (up to C) can be stat()-ed,
	 * then obviously the path is OK and we can continue.
	 * Otherwise, if the error was that P did not exist, 
	 * then create it.
	 */
	if (stat(P, &Statbuf)<0)    /* error during access */
	{
	    if (errno==ENOENT)	/* does not exist --> create */
	    {
		errno=0;
		if (mkdir(P, DIR_MODE)<0)   /* create path component now */
		{
		    cerr<<"\n? mkdir_p(): Problem creating dir \""<<P<<"\": "
			<<strerror(errno)<<endl;
		    delete [] P; return(errno);
		}
	    }
	    else    /* real trouble, bail out */
	    {
		cerr<<"\n? mkdir_p(): Problem creating dir \""<<P<<"\": "
		    <<strerror(errno)<<endl;
		delete [] P; return(errno);
	    }
	}
	
	/* replace slash and trundle on */
	*C='/';
    }	    /* for */
    free(P);
    return(0);
}
/* END of mkdir_p() */

/* pdb_result(): saves the result of the simulation in a file Pdbf
 * provided the coordinates in Xyz are 3-dimensional and there are
 * no dimension mismatches. Returns 1 on success, 0 on error.
 */
int pdb_result(const char *Pdbf, const Points_& Xyz, 
	const Polymer_& Model, const Pieces_& Pieces, 
	const Scores_& Bestsco)
{
    // paranoia
    if (Xyz.active_len() != Model.len()+2)
    {
	cerr<<"\n? pdb_result(): No. of coordinates ("<<Xyz.active_len()
	    <<") does not match model chain length ("<<Model.len()
	    <<"), cannot write PDB file\n"<<flush;
	return(0);
    }
    if (Xyz.dim()!=3)
    {
	cerr<<"\n? pdb_result(): Coordinates are not 3D,  cannot write PDB file\n"<<flush;
	return(0);
    }
    
    // make the PDB entry (a C struct, cf. "pdbprot.h")
    Pdbentry_ *Entry=make_pdbentry(Xyz, Model, Pieces);
    
    // generate remarks
    const unsigned int REMARK_NO=5, REMARK_LEN=61;
    unsigned int ri;
    char **Remarks=new char* [REMARK_NO];
    for (ri=0; ri<REMARK_NO; ri++)
	Remarks[ri]=new char [REMARK_LEN];
    
    sprintf(Remarks[0], "BOND SCORE: %.3e", Bestsco[Scores_::BOND].score());
    sprintf(Remarks[1], "BUMP SCORE: %.3e", Bestsco[Scores_::NONBD].score());
    sprintf(Remarks[2], "EXTERNAL RESTRAINT SCORE: %.3e", Bestsco[Scores_::RESTR].score());
    sprintf(Remarks[3], "SECONDARY STRUCTURE SCORE: %.3e", Bestsco[Scores_::SECSTR].score());
    sprintf(Remarks[4], "ACCESSIBILITY SCORE: %.3e", Bestsco[Scores_::ACCESS].score());
    
    // write to disk (this is a C utility from "pdbprot") 
    put_pdb(Pdbf, Entry, Remarks, REMARK_NO);
    
    // clean up, using the C utility fn from "pdbprot" for the malloc/free symmetry
    for (ri=0; ri<REMARK_NO; ri++)
	delete [] Remarks[ri];
    delete [] Remarks;
    free_pdb(Entry);
    
    return(1);	// OK
}
// END of pdb_result()

/* make_pdbentry(): constructs a Pdbentry_ struct and returns a ptr to it.
 * The coordinates are in Xyz, the chain description in Model, the secondary
 * structure description is supplied by Pieces. Note that Pdbentry_ and the
 * sub-structures within are C structs and malloc/calloc/free should be used
 * to create and destroy them. Xyz is supposed to have been masked to 
 * "full access" before the call. 
 */
static Pdbentry_ *make_pdbentry(const Points_& Xyz, const Polymer_& Model, 
	const Pieces_& Pieces)
{
    // create the entry, C style
    Pdbentry_ *Entry=(Pdbentry_ *) malloc(sizeof(Pdbentry_));
    
    // init some of the easier fields, zero others
    strcpy(Entry->Header, "PROTEIN MODEL");
    
    time_t Now=time(NULL);  // get current date
    strftime(Entry->Date, 10, "%d-%b-%y", localtime(&Now));
    
    strcpy(Entry->Pdbcode, "0DRG");
    Entry->Compound=(char *) calloc(61, sizeof(char));
    strcpy(Entry->Compound, "MODEL C-ALPHA:FAKE C-BETA CHAIN");
    Entry->Source=(char *) calloc(61, sizeof(char));
    strncpy(Entry->Source, version_string(), 60); // version string may be truncated
    strcpy(Entry->Expdta, "THEORETICAL MODEL");
    Entry->Resol=-1.0;
    Entry->Hbonds=NULL; Entry->Hbno=0;
    Entry->Ssbs=NULL; Entry->Ssbno=0;
    
    // create the chain
    Entry->Chainno=1;
    Entry->Chains=(Chain_*) malloc(sizeof(Chain_));
    
    // no H-bond and S-S bond list
    Entry->Chains->Hbonds=NULL; Entry->Chains->Hbno=0;
    Entry->Chains->Ssbs=NULL; Entry->Chains->Ssbno=0;
    
    // generate the amino acid sequence
    int i, Rno=Model.len();	// ignore N/C-terminal pseudo-alphas
    Entry->Chains->Aano=Rno;
    Entry->Chains->Seq=(char *) calloc(Rno+1, sizeof(char));
    for (i=0; i<Rno; i++)
	Entry->Chains->Seq[i]=Model.aa(i);
    Entry->Chains->Seq[Rno]='\0';   // para...
    
    // set chain ID and chain type
    Entry->Chains->Chid=' ';
    Entry->Chains->Type='P';	// "protein"
    
    // convert the C++ secondary structure info into C
    make_secs(Entry->Chains, Pieces);
    
    // store the C++ coordinates
    make_atoms(Entry->Chains, Xyz, Model);
    
    // return the full C structure
    return(Entry);
}
// END of make_pdbentry()

/* make_secs(): stores the secstr information contained in Pieces in the
 * C structure pointed to by Chain. Chain->Seq must contain the sequence
 * of the model. 
 */
static void make_secs(Chain_ *Chain, const Pieces_& Pieces)
{
    Clist1_<Sstr_> Slist=Pieces.secs();
    unsigned int Slen=Slist.len();
    if (!Slen)	// no secondary structure, exit
    {
	Chain->Secs=NULL; Chain->Secsno=0;
	return;
    }
    
    /* There is one Secstr_ element for each helix and strand in the
     * Secs array. Alloc it to the number of clusters in Pieces which is
     * guaranteed to be >= than the strand+helix number. Should be revised
     * if clusters don't coincide with secstr features.
     */
     Secstr_ *S=(Secstr_*) calloc(Pieces.clu_no(), sizeof(Secstr_));
    
    // process all secstr in Pieces
    unsigned int i, j, k, Strno, Hc=0, Sc=0;
    int T, O;
    Strand_ St;
    
    for (Slist.begin(), i=k=0; i<Slen; i++, Slist++)
    {
	/* Instead of devising an elegant virtual mechanism for converting
	 * the Helix_ and Beta_ objects, a quick pedestrian approach is chosen.
	 * Helices will be detected by the is_helix() method returning 1.
	 * Note: the secstr limits are already in [1..Rno] range.
	 */
	Strno=(*Slist)->strand_no();
	if ((*Slist)->is_helix())	// assumed to be helical
	{
	    S[k].Sectype=HELIX;
	    S[k].No=++Hc;
	    sprintf(S[k].Id, "H%d", Hc);
	    
	    S[k].Beg=((Helix_ *)&(*Slist))->beg();	// access as helix
	    S[k].End=((Helix_ *)&(*Slist))->end();
	    
	    S[k].Chid=S[k].Begrid=S[k].Endrid=' ';
	    S[k].Begaa=Chain->Seq[S[k].Beg-1];
	    S[k].Endaa=Chain->Seq[S[k].End-1];
	    
	    // helix types: right-handed 3/10, alpha, pi
	    switch(((Helix_ *)&(*Slist))->helix_type())
	    {
		case Helix_::HX310: S[k].Type=5; break;	    // 3/10
		case Helix_::HXPI: S[k].Type=3; break;	    // pi
		case Helix_::ALPHA: 
		default: S[k].Type=1; break;	    // alpha
	    }
	    k++;
	}
	else	// assumed to be a beta-sheet
	{
	    // cycle through the strands
	    for (j=0; j<Strno; j++, k++)
	    {
		St= (*(Beta_ *)&(*Slist))[j];   // access as beta: get j-th strand
		S[k].Sectype=SHEET;
		S[k].No=j+1;	    // this is the strand no. within the sheet
		sprintf(S[k].Id, "S%d", Sc+1);	// this is the sheet count
		S[k].Beg=St.beg();
		S[k].End=St.end();
		S[k].Chid=S[k].Begrid=S[k].Endrid=' ';
		S[k].Begaa=Chain->Seq[St.beg()-1];
		S[k].Endaa=Chain->Seq[St.end()-1];
		S[k].Type=St.sense();	// 0 for first, +1 for par, -1 for anti
		
		// beta-specific info
		S[k].Strandno=Strno;
		strcpy(S[k].Thisat, " CA ");
		strcpy(S[k].Otherat, " CA ");
		S[k].Thisrid=S[k].Otherid=S[k].Otherchid=' ';
		
		// phasing: only for 2nd, 3rd, ...
		if (!j) continue;
		
		/* get the first residue which has a "previous" neighbour 
		 * and use them as the "This" and "Other" residues
		 */
		for (T=St.beg(); T<=St.end() && (O=(*Slist)->hbond_prev(T))<0; T++);
		S[k].This=T; S[k].Other=O;
		S[k].Thisaa=Chain->Seq[T];
		S[k].Otheraa=Chain->Seq[O];
	    }
	    Sc++;
	}	// helix? sheet?
    }	    // for (Slist)
    
    // re-size the S array and set Chain
    Chain->Secs=(Secstr_ *) realloc(S, k*sizeof(Secstr_));
    Chain->Secsno=k;
}
// END of make_secs()

/* make_atoms(): constructs the coordinate array in Chain from the C-alpha coordinates
 * in Xyz (should be fully unmasked, preferably 3D throughout). Fake C-beta
 * coordinates are added to the C-alphas in the Chain->Atoms array. The
 * sequence must already be present in Chain->Seq.
 */
static void make_atoms(Chain_ *Chain, const Points_& Xyz, const Polymer_& Model)
{
    // make the fake C-betas
    unsigned int Rno=Xyz.len()-2;   // 
    Points_ Beta(Rno+2, 3);
    Fakebeta_::beta_xyz(Xyz, Model, Beta);
    
    // alloc the atom array
    unsigned int i, k, Ano=2*Rno+2, Glyno=0;
    Atom_ *A=(Atom_ *) calloc(Ano, sizeof(Atom_));
    
    // 0th atom is the terminal N
    A[0].Atno=1;
    strcpy(A[0].Id, " N  ");
    A[0].Alt=A[0].Rid=' ';
    A[0].Aa=Chain->Seq[0];
    A[0].Resno=1;
    A[0].X=Xyz[0][0]; A[0].Y=Xyz[0][1]; A[0].Z=Xyz[0][2];
    A[0].Occu=Model.cons(0); A[0].Bfact=Model.phob(0);
    
    // scan the chain, store coords
    for (i=k=1; i<=Rno; i++, k++)
    {
	// C-alpha
	A[k].Atno=k+1;
	strcpy(A[k].Id, " CA ");
	A[k].Alt=A[k].Rid=' ';
	A[k].Aa=Chain->Seq[i-1];
	A[k].Resno=i;
	A[k].X=Xyz[i][0]; A[k].Y=Xyz[i][1]; A[k].Z=Xyz[i][2];
	A[k].Occu=Model.cons(i-1); A[k].Bfact=Model.phob(i-1);
	
	// fake C-beta (skip for Glys!)
	if (A[k].Aa=='G') { Glyno++; continue; }
	
	k++;
	A[k].Atno=k+1;
	strcpy(A[k].Id, " CB ");
	A[k].Alt=A[k].Rid=' ';
	A[k].Aa=Chain->Seq[i-1];
	A[k].Resno=i;
	A[k].X=Beta[i][0]; A[k].Y=Beta[i][1]; A[k].Z=Beta[i][2];
	A[k].Occu=Model.cons(i-1); A[k].Bfact=Model.phob(i-1);
    }

    // last atom is the terminal C, i==Rno+1 automatically
    A[k].Atno=k+1;
    strcpy(A[k].Id, " C  ");
    A[k].Alt=A[k].Rid=' ';
    A[k].Aa=Chain->Seq[i-2];
    A[k].Resno=Rno;
    A[k].X=Xyz[i][0]; A[k].Y=Xyz[i][1]; A[k].Z=Xyz[i][2];
    A[k].Occu=Model.cons(Rno-1); A[k].Bfact=Model.phob(Rno-1);
    
    // correct for missing C-betas on Gly-s
    Ano-=Glyno;
    A=(Atom_ *) realloc(A, Ano*sizeof(Atom_));
    
    Chain->Atoms=A;
    Chain->Atomno=Ano;
}
// END of make_atoms()

// ==== END OF FUNCTIONS Output.c++ ====
