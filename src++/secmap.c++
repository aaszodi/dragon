// ==== PROGRAM secmap.c++ ====

/* Generates a mapping of secondary structures onto
 * a target sequence using a multiple alignment.
 * The secstr assignments come from DSSP files.
 */

// SGI C++ 7.1, 27-Jun-1998. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <iostream.h>
#include <iomanip.h>

// ---- MODULES ----

#include "Align.h"	/* for alignment description: from DRAGON */
#include "Secmap.h"	/* classes for secstr mapping */
#include "dsspread.h"	/* DSSP file I/O in C */

// ---- PROTOTYPES ----

static int dssp_secmap(const char *Dsspnm, int Knownno, const Align_& Align, 
	int Master, Secmap_ *Maps, int Maplen);

// ==== MAIN ====

/* The program takes the following arguments:-
 * the alignment file name, the no. of the target sequence within
 * the alignment and a list of DSSP files corresponding to known
 * structures within the alignment. The alignment name, the
 * target seq number and at least one DSSP filename are mandatory.
 * Example: secmap alignment.aln 2 1lfs.dssp 2lfs.dssp
 * Output is to stdout.
 */
int main(int argc, char *argv[])
{
    if (argc<4)
    {
	cerr<<"\n! Usage: "<<argv[0]
	    <<" alignment_file target_no DSSP_file [DSSP_file ...]\n";
	exit(EXIT_FAILURE);
    }
    
    // get the alignment
    Align_ Align;
    
    if (!Align.read_file(argv[1]))
    {
	cerr<<"\n! "<<argv[0]<<": Cannot read alignment file \""<<argv[1]<<"\"\n";
	exit(EXIT_FAILURE);
    }
    cout<<"# Alignment file: "<<argv[1]<<endl;
    
    // get the master idx
    int Master=-1;
    Master=atoi(argv[2]);	// convert text to number
    if (Master<=0 || Master>Align.seq_no())
    {
	cerr<<"\n! "<<argv[0]<<": Master="<<Master<<" is out of range [1.."
	    <<Align.seq_no()<<"]\n";
	exit(EXIT_FAILURE);
    }
    cout<<"# Target sequence number = "<<Master<<endl;
    Master--;	// move into [0..Seqno-1] range
    
    // prepare storage for the mapping
    int i, k, Maplen=Align.len();	// total alignment length with gaps
    Secmap_ *Maps=new Secmap_ [Maplen];
    const char *Posstr;
    char Targetaa;
    
    // store the target sequence with the gaps
    for (i=0, k=1; i<Maplen; i++, k++)
    {
	Posstr=Align.pos(i);
	Targetaa=Posstr[Master];
	if (Targetaa!='-')
	    Maps[i].set_aa(Targetaa, k, argc-3);
	else
	{
	    Maps[i].set_aa('-', 0, argc-3);
	    k--;
	}
    }
    
    // scan the DSSP files
    int Knownno, Dsfx, Dspos;
    for (Knownno=0, Dsfx=3; Dsfx<argc; Knownno++, Dsfx++)
    {
	if (!(Dspos=dssp_secmap(argv[Dsfx], Knownno, Align, Master, Maps, Maplen)))
	{
	    {
		cerr<<"\n? "<<argv[0]<<": Cannot process DSSP file \""<<argv[Dsfx]<<"\", skipped\n";
		Knownno--; continue;
	    }
	}
	cout<<"# Scaffold: \""<<argv[Dsfx]<<"\", sequence number = "<<Dspos<<endl;
    }
    
    // output to stdout
    for (i=0; i<Maplen; i++)
	cout<<Maps[i];
    
    delete [] Maps;
    exit(EXIT_SUCCESS);
}

// ==== GLOBAL FUNCTIONS ====

/* dssp_secmap(): makes sure that the DSSP in Dsspnm belongs to a
 * sequence in Align, then extracts all secstr info which
 * can be mapped onto the target sequence identified by Master, 
 * and stores them in Maps (Maplen long, with gaps etc).
 * Knownno provides the index of the currently processed DSSP file.
 * Returns 0 if the DSSP file cannot be
 * read or the sequence is not in the alignment, the no. of
 * the scaffold sequence in the alignment otherwise.
 */
static int dssp_secmap(const char *Dsspnm, int Knownno, const Align_& Align, 
	int Master, Secmap_ *Maps, int Maplen)
{
    // get the DSSP info
    unsigned int Dsspsize=0, Chainno=0;
    Dssprec_ *Dssp=dssp_read(Dsspnm, &Dsspsize, &Chainno);
    if (Dssp==NULL || !Dsspsize) return(0);
    
    // check sequence
    int i, k, Dssno;
    char Aa;
    char *Alnseq=NULL, *Dsspseq=new char [Dsspsize+1];
    for (i=k=0; i<Dsspsize; i++)
    {
	Aa=Dssp[i].Res;
	if (Aa=='!')
	{
	    cerr<<"\n? dssp_secmap(\""<<Dsspnm
		<<"\", ...): Chain break at pos="<<(i+1)<<endl;
	    continue;
	}
	Dsspseq[k++]=Aa;	// copy residue codes
    }
    Dsspseq[k]='\0'; // paranoid termination
    for (Dssno=0; Dssno<Align.seq_no(); Dssno++)
    {
	Align.seq(Dssno, Alnseq);	// extract Dssno:th sequence
	if (!strncmp(Alnseq, Dsspseq, k))
	    break;  // Dssno:th sequence is the same, OK
    }
    delete [] Alnseq; delete [] Dsspseq;
    if (Dssno>=Align.seq_no())
    {
	cerr<<"\n? dssp_secmap(): sequence from \""<<Dsspnm<<"\" is not in alignment\n";
	free(Dssp);
	return(0);
    }
    
    // walk along the alignment
    int p, st, sd;
    Smap_ Smap;
    for (i=p=0; p<Maplen; p++, i++)
    {
	st=Align.seq_pos(Master, p);	 // the target's resno if >=0 or gap
	sd=Align.seq_pos(Dssno, p);	// the DSSP seq's resno
	if (sd>=0)
	    for (i=0 ; i<Dsspsize && Dssp[i].Resno!=sd+1; i++);	// skip DSSP chain breaks
	
	// The i:th DSSP record corresponds now to the p-th alignment
	// position which contains the st:th target residue and the sd:th
	// scaffold residue.
	if (sd<0)
	{
	    Smap.set_nonbeta(Smap_::GAP);   // DSSP gapped at position p
	}
	else
	{
	    int s1, p1, s2, p2;
	    
	    switch(Dssp[i].Secstruct)
	    {
		case 'G':	// 3/10 helix
		    Smap.set_nonbeta(Smap_::HELIX_310);
		break;
		case 'H':	// alpha-helix
		    Smap.set_nonbeta(Smap_::HELIX_AL);
		break;
		case 'I':	// pi-helix
		    Smap.set_nonbeta(Smap_::HELIX_PI);
		break;
		case 'E':	// beta-sheet
		    /* get the partner residues mapped onto
		     * the target sequence. s1, s2 will be <0 if
		     * non-mappable
		     */
		    if (Dssp[i].Beta1)	// there was a beta partner
		    {
			p1=Align.align_pos(Dssno, Dssp[i].Beta1-1);
			s1=Align.seq_pos(Master, p1);   // target resno = 1st beta partner
		    }
		    else s1=-1;
		    if (Dssp[i].Beta2)
		    {
			p2=Align.align_pos(Dssno, Dssp[i].Beta2-1);
			s2=Align.seq_pos(Master, p2);   // target resno = 2nd beta partner
		    }
		    else s2=-1;
		    if (s1<0 && s2<0)   // no beta partners, assign OTHER
			Smap.set_nonbeta(Smap_::OTHER);
		    else Smap.set_beta(s1+1, s2+1,
			    isupper(Dssp[i].Bridge1), isupper(Dssp[i].Bridge2),
			    Dssp[i].Sheet);
		break;
		default:	// any other assignment
		    Smap.set_nonbeta(Smap_::OTHER);
		break;
	    }	    // switch
	}
	
	// store as Knownno:th mapping for the p:th position
	Maps[p].set_struct(Knownno, Smap);
    }
    
    free(Dssp);	/* was allocated in a C routine */
    return(Dssno+1);	// DSSP sequence in the alignment
}
// END of dssp_secmap()

// ==== END OF PROGRAM secmap.c++ ====
