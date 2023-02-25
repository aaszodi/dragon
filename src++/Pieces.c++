// ==== PROJECT DRAGON: METHODS Pieces.c++ ====

/* Class for keeping track of secondary structures and
 * general segments for hierarchic projection.
 */

// SGI C++ 7.1, IRIX 6.2, 2. Apr. 1997. Andris Aszodi

// ---- STANDARD HEADERS ----

#include "fstream.h"

// ---- CLASS HEADER ----

#include "Pieces.h"

// ==== Pieces_ METHODS ====

// ---- Constructors ----

/* Inits the object to accept secondary structures within a Resno-long
 * chain. Note that Resno must always be specified explicitly to avoid
 * confusion and therefore there's no default value for it. Also, the
 * default ctor is disallowed (see "forbidden methods").
 */
Pieces_::Pieces_(unsigned int Resno):
	Rno(Resno), Secsmask(Resno+2), Clus(1), Ctype(1), Changed(false), 
	Secs(), Coils(Linsegm_(0, Resno+1)), Ptclu(NULL)
{
    // NOTE: 0th and Rno+1th points are the N/C termini!
    Clus[0].len(Rno+2); Clus[0].set_values(true);	// one big coil cluster
    Ctype[0]=COIL; make_ptidx();
}

// ---- Size ----

/* res_no(): sets the new size to R. Destroys the secstr list and 
 * represents the chain as one long coil (similar to the ctor).
 * The new secstr assignment can then be built by reading a new
 * specification.
 */
void Pieces_::res_no(unsigned int R)
{
    Secs.clear(); Coils.clear();    // explicit destruction of lists
    Rno=R;
    Clus.len(1); Ctype.len(1);
    Clus[0].len(Rno+2); Clus[0].set_values(true);	// one big cluster
    Coils+=Linsegm_(0, Rno+1);	// one coil spanning all points, including N/C termini
    Ctype[0]=COIL;
    Secsmask.len(Rno+2); Secsmask.set_values(false);  // no secstr
    make_ptidx();
    Changed=false;
}
// END of res_no()

// ---- Internal consistency maintenance ----

/* make_coils(): given a chain size, a valid list of secondary structures
 * in Secs and a valid secondary structure membership bitmask, 
 * the coil list and the cluster array are updated here. Should be called whenever
 * the secstr list or the chain size is modified. Private
 */
void Pieces_::make_coils()
{
    if (!Changed) return;   // do nothing if update has been done already
    
    // get secstr segment no
    unsigned int i, j, Ssegno=0;
    int Secmaskno;
    
    /* Scan all items in the secstr list for their masks.
     * Note that overlapping beta-sheets are allowed. If an overlap
     * is detected between a sheet and a previous cluster mask in the cluster
     * array, then the mask of the current sheet is ORed to
     * that cluster mask. This assumes that helix:helix and 
     * helix:sheet overlaps have already been filtered out by
     * the input methods.
     */
    Bits_ Smask;
    Clus.len(Secs.len());   // reset cluster array to store secstr masks only
    Ctype.len(Secs.len());  // reset cluster type array to the same size
    
    for (Secmaskno=0, Secs.begin(); Secs!=NULL; Secmaskno++, Secs++)
    {
	Smask=(*Secs)->mask(Rno+2);	// get the mask w/ the right size
	// check if sheet overlaps w/ another cluster mask already there
	if ((*Secs)->is_beta())
	{
	    for (j=0; j<Secmaskno; j++)
	    {
		if (Smask.on_no()+Clus[j].on_no()!=(Smask | Clus[j]).on_no())
		    break;	// overlap found
	    }
	    if (j==Secmaskno)	// no overlap
	    {
		Clus[Secmaskno]=Smask;
		Ctype[Secmaskno]=SHEET;
	    }
	    else
	    {
		Clus[j]|=Smask;    // overlap, OR them together
		Secmaskno--;	// reset counter by 1
	    }
	}
	else 	// store helices in cluster descriptor array
	{
	    Clus[Secmaskno]=Smask;
	    Ctype[Secmaskno]=HELIX;
	}
	
	Ssegno+=(*Secs)->strand_no();	// get individual segment numbers
	
	/* NOTE: call the secstr ideal structure update method
	 * here manually. Cf. notes in "Secstr.h"
	 */
	(*Secs)->make_idstruct();
    }
    
    /* The coils will be the contiguous 00..0 regions in the 
     * secstr mask. Walk over it and collect them
     */
    int B=-1, E;
    unsigned int Cno=Secmaskno;  // index to where the coils should be put
    Linsegm_ Coil;	// temporary
    
    Coils.clear();	// delete previous list of coils
    Clus.len(Secmaskno+Ssegno+1);	// resize cluster array (Ssegno overestimates if overlaps)
    Ctype.len(Secmaskno+Ssegno+1);	// and cluster types as well
    for (i=0; i<Rno+2; i++)	// N/C termini are always coil parts
    {
	if (Secsmask.get_bit(i))
	{
	    if (B>=0)	// just finished non-secstr region
	    {
		Coil=Linsegm_(B, E);
		Coils+=Coil;	// add new coil to the end of the list
		Clus[Cno]=Coil.mask(Rno+2); // store its mask
		Ctype[Cno++]=COIL;	    // and its type
		B=-1;
	    }
	    continue;	// otherwise just skip
	}
	else	// we are in a non-secstr region now
	{
	    if (B<0) B=i;   // restart
	    E=i;    // always mark the last false bit seen
	}
    }
    
    // finish up the last open non-secstr segment
    if (B>=0)
    {
	Coil=Linsegm_(B, E);
	Coils+=Coil;	// add last coil 
	Clus[Cno]=Coil.mask(Rno+2); // store its mask
	Ctype[Cno++]=COIL;	    // and its type
    }
    Clus.len(Cno);  // readjust cluster array size
    make_ptidx();
    
    Changed=false;  // update is complete now
}
// END of make_coils()

/* make_ptidx(): constructs the internal index array Ptclu.
 * Ptclu[i] is the index of the cluster that contains the i:th point.
 * Private
 */
void Pieces_::make_ptidx()
{
    if (Ptclu!=NULL) delete [] Ptclu;
    Ptclu=new int [Rno+2];
    
    register int i, ci;
    
    for (i=0; i<Rno+2; i++)
    {
	for (ci=0; ci<clu_no() && !Clus[ci].get_bit(i); ci++);
	if (ci>=clu_no())
	    cerr<<"\n! Pieces_::make_ptidx(): Point "<<i<<" is not found in any of the clusters\n";
	Ptclu[i]=ci;
    }
}

// ---- Input/output ----

/* read_secstr(): reads the secondary structure specification from
 * file Secf. If Secf=="" or NULL then the secondary structure
 * layout will be reset to all-coil.
 * Returns 1 if OK, -1 on reset, 0 on error.
 */
int Pieces_::read_secstr(const char *Secf)
{
    if (Secf==NULL || !strlen(Secf))	// reset to all-coil
    {
	res_no(Rno);	// this does the reset
	return(-1);
    }
    
    ifstream Inf(Secf);
    
    if (!Inf)
    {
	cerr<<"\n? Pieces_::read_secstr(\""<<Secf<<"\"): Cannot open\n";
	return(0);
    }
    
    Inf>>(*this);
    Inf.close();
    return(Inf.good() || Inf.eof());
}
// END of read_secstr()

/* >>: reads the secondary structure specification from a stream In
 * into the Pieces_ object P which is supposed to be pre-initialised.
 * Alpha-helices and beta sheets can be specified according to the
 * following syntax (also cf. the Helix_ and Sheet_ input operator comments):-
 * "<helix> <beg> <end> \n" (where <helix>=["HELIX"|"ALPHA"|"HX310"|"HXPI"] )
 * "SHEET\n"		    (or for beta-sheets)
 * "STRAND <beg> <end>\n"
 * "STRAND <beg> <end> [PAR|ANTI] <this> <other>\n"
 * .....
 * "END\n"
 * Helices and sheets may be freely mixed and comment lines (lines
 * beginning with '#') may be present between the descriptions.
 * Note that comments are not allowed between "SHEET"/"END" pairs.
 * The routine builds a temporary secstr list and updates P only
 * if the input was successful. Checks are made to ensure that
 * all items fit into the chain (based on the internal Rno value)
 * and that there is no overlap between helices and other helices/sheets.
 * Sheet/sheet overlap is allowed but a warning is printed.
 */
istream& operator>>(istream& In, Pieces_& P)
{
    if (!In)
    {
	cerr<<"\n? >>Pieces_: Cannot read from stream\n";
	return(In);
    }
    
    static const unsigned int BUFLEN=132;   // alloc input line buffer
    char Buf[BUFLEN+1];
    streampos Linbeg=0, Linend=In.tellg();  // current pos
    Helix_ Htemp;   // temporary structures and list
    Beta_ Btemp;
    List1_<Sstr_> Templist;
    Bits_ Hsmask(P.Rno+2), Secsmask(P.Rno+2, false);	// new secondary structure mask
    
    // read the stream line-by-line
    while(Linbeg=Linend, In.seekg(Linend), 
	    memset(Buf, 0, BUFLEN+1), In.getline(Buf, BUFLEN).good())
    {
	Linend=In.tellg();  // beginning of next line
	if (!strlen(Buf) || Buf[0]=='#')    // skip empty lines and comments
	    continue;
	
	if (NULL!=strstr(Buf, "SHEET"))	// beta-sheet should come
	{
	    In.seekg(Linbeg);  // go back to beginning of the line
	    In>>Btemp;	    // attempt to read the sheet
	    Linend=In.tellg();	// sheets span several input lines
	    if (!In)
	    {
		cerr<<"\n? >>Pieces_: Cannot parse into sheet:\n"<<Buf<<endl;
		In.clear(~ios::failbit & In.rdstate());	// reset "fail" bit
		continue;
	    }
	    	    
	    // check size
	    Hsmask=Btemp.mask();    // the length will correspond to last "ON"
	    if (Hsmask.len()>P.Rno+1)
	    {
		cerr<<"\n? >>Pieces_: Sheet does not fit, ignored\n";
		continue;
	    }
	    
	    // check overlap (allow sheet:sheet overlaps)
	    Hsmask.len(P.Rno+2);	// adjust mask length
	    if (Hsmask.on_no()+Secsmask.on_no()!=(Hsmask | Secsmask).on_no())
	    {
		// Overlapped with somebody. Check if sheet or helix
		Bits_ Tmask(P.Rno+2);
		bool Helixoverlap=false;
		for (Templist.begin(); Templist!=NULL && !Helixoverlap; Templist++)
		{
		    if ((*Templist)->is_beta()) continue;   // don't test against sheets
		    Tmask=(*Templist)->mask(P.Rno+2);
		    if (Tmask.on_no()+Hsmask.on_no()!=(Tmask | Hsmask).on_no())
		    {
			cerr<<"\n? >>Pieces_: Sheet overlaps w/ helix, ignored\n";
			Helixoverlap=true;
			break;
		    }
		}
		if (Helixoverlap) continue; // don't accept, get next secstr
		else
		    cerr<<"\nWARNING: >>Pieces_: Sheet overlaps w/ other sheet(s), bifurcation assumed\n";
	    }
	    
	    // looks OK, add to temp list and update secstr mask
	    Templist+=Btemp;
	    Secsmask|=Hsmask;
	    continue;
	}
	else	// helical input
	{
	    In.seekg(Linbeg);  // go back to beginning of the line
	    if (!(In>>Htemp))	    // read in the helix
	    {
		cerr<<"\n? >>Pieces_: Cannot parse into helix:\n"<<Buf<<endl;
		In.clear(~ios::failbit & In.rdstate());	// reset "fail" bit
		continue;
	    }
	    
	    // check size
	    if (Htemp.end()>P.Rno)
	    {
		cerr<<"\n? >>Pieces_: Helix does not fit, ignored\n"<<Buf<<endl;
		continue;
	    }
	    
	    // check overlap (the mask spans the N/C termini as well)
	    Hsmask=Htemp.mask(P.Rno+2);
	    if (Hsmask.on_no()+Secsmask.on_no()!=(Hsmask | Secsmask).on_no())
	    {
		cerr<<"\n? >>Pieces_: Helix overlaps w/ other secstr, ignored\n"<<Buf<<endl;
		continue;
	    }
	    
	    // looks OK, add to temp list and update secstr mask
     	    Templist+=Htemp;
	    Secsmask|=Hsmask;
	    continue;
	}
    }	    // while
    
    // update P if there was at least 1 valid secstr
    if (Templist.len()>=1)
    {
	P.Secs=Templist; P.Secsmask=Secsmask;
	P.Changed=true; P.make_coils();
    }
    return(In);
}
// END of >>

/* <<: lists the secondary structure elements to Out. */
ostream& operator<<(ostream& Out, const Pieces_& P)
{
    Clist1_<Sstr_> Sl=P.secs();
    for (Sl.begin(); Sl!=NULL; Sl++)
	Out<<(*(*Sl));
    return(Out);
}
// END of <<

// ==== END OF METHODS Pieces.c++ ====
