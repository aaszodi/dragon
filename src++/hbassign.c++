// ==== PROGRAM hbassign.c++ ====

/* Generates H-bond assignments for a target structure
 * using a MULTAL-style alignment and DSSP files of
 * template structures. Use for DRAGON homology models
 * to provide H-bond restraints for refinement with QUANTA.
 */

// SGI C++ 4.0, IRIX 5.3, 26. June 1996. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <strstream.h>
#include <fstream.h>

// ---- MODULES ----

#include "Align.h"	/* for alignment description: from DRAGON */
#include "List1.h"	/* 2-way lists */
#include "dsspread.h"	/* DSSP file I/O in C */

// ---- DEFINITIONS ----

#define LINELEN 132	/* control file line length */
#define HBEN_MAX (-2.0)	/* H-bond energy threshold (kcal/mol) */

// ==== CLASSES ====

/* Class Hb_: stores the main-chain H-bond data for
 * the target molecule.
 */
class Hb_
{
    // data
    private:
    
    unsigned int Nh, Co;    // residue nos [1..x] for donor and acceptor
    int Filenum;	    // no. of different DSSP files containing this pair
    double En;	    // best bond energy in kcal/mol
    
    // methods
    public:
    
    // default ctor
    Hb_(unsigned int Nhno=0, unsigned int Cono=0, double Energy=0.0): 
	Nh(Nhno), Co(Cono), Filenum(Nhno && Cono), En(Energy) {}
    
    int filenum() const { return(abs(Filenum)); }
    void prime_filenum() { if (Filenum>0) Filenum*= -1; }
    double energy() const { return(En); }
    
    unsigned int add_bond(unsigned int Nhno, unsigned int Cono, double Energy);
    friend ostream& operator<<(ostream& Out, const Hb_& Hb);
    
};
// END OF CLASS Hb_

// ---- PROTOTYPES ----

static int dssp_hblist(const char *Dsspnm, const Align_& Align, 
	int Master, List1_<Hb_>& Hblist);
static void map_bond(List1_<Hb_>& Hblist, int Nhno, int Cono, double En);

// ==== MAIN ====

/* The program takes the name of an ASCII file which specifies
 * the alignment file name, the no. of the target sequence within
 * the alignment and a list of DSSP files corresponding to known
 * structures within the alignment. No empty lines; one param per line.
 * Example format:-
 * 
 * align.seq
 * 2
 * 1ABC.dssp
 * data/2DEF.dssp
 * ....
 * Output is to stdout in QUANTA/CHARMm *.dcn file format.
 */
int main(int argc, char *argv[])
{
    if (argc<2)
    {
	cerr<<"\n! Usage: "<<argv[0]<<" control_file\n";
	cerr<<"where the control file has the format:\n";
	cerr<<"\t<MULTAL_file>\n\t<target_no>\n\t<DSSP_file>\n\t...\n";
	exit(EXIT_FAILURE);
    }
    
    // open the control file
    ifstream Ctrlf(argv[1]);
    if (!Ctrlf)
    {
	cerr<<"\n! "<<argv[0]<<": Cannot open control file \""<<argv[1]<<"\"\n";
	exit(EXIT_FAILURE);
    }
    
    // get the alignment
    char Alnfnm[LINELEN], Line[LINELEN];
    Align_ Align;
    
    Alnfnm[0]='\0';
    Ctrlf.getline(Alnfnm, LINELEN, '\n');
    if (!Align.read_file(Alnfnm))
    {
	cerr<<"\n! "<<argv[0]<<": Cannot read alignment file \""<<Alnfnm<<"\"\n";
	exit(EXIT_FAILURE);
    }
    
    // get the master idx
    int Master=0;
    Line[0]='\0';
    Ctrlf.getline(Line, LINELEN, '\n');
    istrstream(Line)>>Master;	// convert text to number
    if (Master<=0 || Master>Align.seq_no())
    {
	cerr<<"\n! "<<argv[0]<<": Master="<<Master<<" is out of range [1.."
	    <<Align.seq_no()<<"]\n";
	exit(EXIT_FAILURE);
    }
    Master--;	// move into [0..Seqno-1] range
    
    // construct the list of H-bonds, write comments to stdout
    List1_<Hb_> Hblist;
    int Knownno;  // no. of known structures (valid DSSP stuff)
    
    cout<<"*CHARMm distance constraints faked by "<<argv[0]<<endl;
    cout<<"*Control file: \""<<argv[1]<<"\"\n";
    cout<<"*DRAGON-IV alignment file: \""<<Alnfnm<<"\"\n";
    cout<<"*Masterno="<<Master<<endl;
    
    for (Knownno=0; Ctrlf.getline(Line, LINELEN, '\n'); Knownno++)
    {
	if (Line[0]=='\0' || Line[0]=='\n')
	{ Knownno--; continue; }
	
	if (!dssp_hblist(Line, Align, Master, Hblist))
	{
	    {
		cerr<<"\n? "<<argv[0]<<": Cannot process DSSP file \""<<Line<<"\", skipped\n";
		Knownno--; continue;
	    }
	}
	cout<<"*DSSP: \""<<Line<<"\"\n";
    }
    Ctrlf.close();
    
    // output to stdout
    cout<<"NOE\nRESET\n";
    for (Hblist.begin(); Hblist!=NULL; Hblist++)
	if (Hblist->filenum()>=Knownno && Hblist->energy()<HBEN_MAX)	// found in all files
	    cout<<(*Hblist);
    cout<<"SCALE     1.0000\nEND\n";    // epilogue
    exit(EXIT_SUCCESS);
}

// ==== Hb_ METHODS ====

/* add_bond(): add the NH->OC bond to the calling object if the
 * indices (Nhno, Cono) match and stores Energy if it is lower
 * than the best seen so far. The number of successful additions
 * is counted as well. Does nothing and returns 0 if the indices
 * don't match; otherwise the no. of successful additions is returned. 
 */
unsigned int Hb_::add_bond(unsigned int Nhno, unsigned int Cono, double Energy)
{
    if (!Filenum)	// virgin, add anyway
    {
	Nh=Nhno; Co=Cono;
	Filenum=1; En=Energy;
	return(1);
    }
    
    if (Nh!=Nhno || Co!=Cono) return(0);    // mismatch
    if (En>Energy) En=Energy;	// save energy if better than prev
    if (Filenum<0)   // there were no additions from current DSSP: see prime_filenum()
    {
	Filenum*=-1; // turn back to >0: no. of files containing this bond so far
	++Filenum;  // and now add 1
    }
    return(Filenum);
}
// END of add_bond()

ostream& operator<<(ostream& Out, const Hb_& Hb)
{
    char Selestr[50];
    const char *Endstr="                                      END -";

    strcpy(Selestr, Endstr);
    ostrstream Sel(Selestr, 50);
    Sel<<"ASSIGN SELE ATOM 0XXX "<<Hb.Co<<" O";
    Out<<Selestr<<endl;
    strcpy(Selestr, Endstr);
    Sel.seekp(0);
    Sel<<"       SELE ATOM 0XXX "<<Hb.Nh<<" HN";
    Out<<Selestr<<"\n  KMIN   25.00 RMIN   1.900 KMAX   25.00 RMAX    2.10\n";
    return(Out);
}

// ==== GLOBAL FUNCTIONS ====

/* dssp_hblist(): makes sure that the DSSP in Dsspnm belongs to a
 * sequence in Align, then scans all main-chain H-bonds which
 * can be mapped onto the target sequence identified by Master, 
 * and stores them in Hblist. Returns 0 if the DSSP file cannot be
 * read or the sequence is not in the alignment, the no. of
 * added H-bonds otherwise.
 */
static int dssp_hblist(const char *Dsspnm, const Align_& Align, 
	int Master, List1_<Hb_>& Hblist)
{
    // get the DSSP info
    unsigned int Dsspsize=0, Chainno=0;
    Dssprec_ *Dssp=dssp_read(Dsspnm, &Dsspsize, &Chainno);
    if (Dssp==NULL || !Dsspsize) return(0);
    
    // check sequence
    unsigned int i, k, Sno;
    char Aa;
    char *Alnseq=NULL, *Dsspseq=new char [Dsspsize+1];
    int *Res=new int [Dsspsize];
    for (i=k=0; i<Dsspsize; i++)
    {
	Aa=Dssp[i].Res;
	if (Aa=='!')
	{
	    cerr<<"\n? dssp_hblist(\""<<Dsspnm
		<<"\", ...): Chain break at pos="<<(i+1)<<endl;
	    continue;
	}
	Res[i]=k;   // i:th position is the k:th amino acid
	Dsspseq[k++]=Aa;	// copy residue codes
    }
    Dsspseq[k]='\0'; // paranoid termination
    for (Sno=0; Sno<Align.seq_no(); Sno++)
    {
	Align.seq(Sno, Alnseq);	// extract s:th sequence
	if (!strncmp(Alnseq, Dsspseq, k))
	    break;  // Sno:th sequence is the same, OK
    }
    delete [] Alnseq; delete [] Dsspseq;
    if (Sno>=Align.seq_no())
    {
	cerr<<"dssp_hblist(): sequence from \""<<Dsspnm<<"\" is not in alignment\n";
	free(Dssp); delete [] Res;
	return(0);
    }
    
    // prime the list: first addition will incr Filenum member in Hb_-s
    for (Hblist.begin(); Hblist!=NULL; Hblist++)
	Hblist->prime_filenum();
    
    // walk along the DSSP array
    int p, s, s2, Respos;
    int Hadd=0;
    for (i=0; i<Dsspsize; i++)
    {
	// skip chain breaks
	if (Dssp[i].Res=='!') continue;
	
	// is this residue main-chain H-bonded at all?
	if (!Dssp[i].Nho[0].Offs && !Dssp[i].Nho[1].Offs &&
	    !Dssp[i].Ohn[0].Offs && !Dssp[i].Ohn[1].Offs)
		continue;   // no H-bonds
	
	Respos=Res[i];	    // there might be gaps
	p=Align.align_pos(Sno, Respos);  // this position contains the i:th DSSP residue
	s=Align.seq_pos(Master, p); // the target's resno in the p:th alignment position
	if (s<0) continue;  // gapped
	
	// get the partners' mapped positions on the Master, add bonds
	// NH-->O first
	if (Dssp[i].Nho[0].Offs)
	{
	    p=Align.align_pos(Sno, Res[Respos+Dssp[i].Nho[0].Offs]);
	    s2=Align.seq_pos(Master, p);
	    if (s2>=0) { map_bond(Hblist, s, s2, Dssp[i].Nho[0].En); Hadd++; }
	}
	if (Dssp[i].Nho[1].Offs)
	{
	    p=Align.align_pos(Sno, Res[Respos+Dssp[i].Nho[1].Offs]);
	    s2=Align.seq_pos(Master, p);
	    if (s2>=0) { map_bond(Hblist, s, s2, Dssp[i].Nho[1].En); Hadd++; }
	}
	// then O<--HN
	if (Dssp[i].Ohn[0].Offs)
	{
	    p=Align.align_pos(Sno, Res[Respos+Dssp[i].Ohn[0].Offs]);
	    s2=Align.seq_pos(Master, p);
	    if (s2>=0) { map_bond(Hblist, s2, s, Dssp[i].Ohn[0].En); Hadd++; }
	}
	if (Dssp[i].Ohn[1].Offs)
	{
	    p=Align.align_pos(Sno, Res[Respos+Dssp[i].Ohn[1].Offs]);
	    s2=Align.seq_pos(Master, p);
	    if (s2>=0) { map_bond(Hblist, s2, s, Dssp[i].Ohn[1].En); Hadd++; }
	}
    }
    
    delete [] Res;
    free(Dssp);	/* was allocated in a C routine */
    return(Hadd);	// no. of newly created bonds
}
// END of dssp_hblist()

/* map_bond(): adds a new H-bond between Nhno and Cono with energy En
 * to the list Hblist. Note that the residue nos. are [0..S-1], so
 * 1 will be added to both. 
 */
static void map_bond(List1_<Hb_>& Hblist, int Nhno, int Cono, double En)
{
    Nhno++; Cono++;
    for (Hblist.begin(); Hblist!=NULL; Hblist++)
    {
	if (Hblist->add_bond(Nhno, Cono, En))
{
// cerr<<"found NH="<<Nhno<<", CO="<<Cono<<": added E="<<En<<endl;
	    return;	// successful addition
}
    }
// cerr<<"started new NH="<<Nhno<<", CO="<<Cono<<": E="<<En<<endl;
    Hblist+=Hb_(Nhno, Cono, En);    // append to end if not found before
}
// END of map_bond()

// ==== END OF PROGRAM hbassign.c++ ====
