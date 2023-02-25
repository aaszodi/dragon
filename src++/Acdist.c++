// ==== PROJECT DRAGON: METHODS Acdist.c++ ====

/* For the storage and retrieval of side-chain atom
 * distances from the C-alpha atom or the sidechain centroid.
 */

// SGI C++, IRIX 6.2, 14. Aug. 1996. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <strstream.h>
#include <fstream.h>

// ---- MODULE HEADER ----

#include "Acdist.h"

// ==== Acdist_ METHODS ====

// ---- Static initialisation ----

const char* Acdist_::AAcodes="ACDEFGHIKLMNPQRSTVWY";

// ---- Access ----

/* ca_dist(), scc_dist(): return the distance of Atom from
 * the C-alpha or the sidechain centroid for amino acid Aa, respectively.
 * -1.0 returned and a warning printed if Aa or Atom was not found.
 */
float Acdist_::ca_dist(char Aa, const String_& Atom) const
{
    int Idx=get_idx(Aa);
    if (Idx<0)
    {
	cerr<<"\n? Acdist_::ca_dist(): Invalid amino acid \'"<<Aa<<"\'\n";
	return(-1.0);
    }
    return(Acss[Idx].ca_dist(Atom));
}

float Acdist_::scc_dist(char Aa, const String_& Atom) const
{
    int Idx=get_idx(Aa);
    if (Idx<0)
    {
	cerr<<"\n? Acdist_::scc_dist(): Invalid amino acid \'"<<Aa<<"\'\n";
	return(-1.0);
    }
    return(Acss[Idx].scc_dist(Atom));
}

// ---- Input ----

/* read_file(): reads sidechain atom distance data from 
 * the file called Fname. Does nothing if Fname is "", NULL (the default)
 * or cannot be opened. Updates only the distances which
 * are explicitly mentioned in the file; for a complete reset, 
 * call reset(). The format of a line in the distance file
 * is:-
 * "AAcode Atomname CAdist CTRdist\n"
 * where AAcode is a 1-letter amino acid code (char), Atomname is
 * an all-uppercase string (must be a PDB-type sidechain atom), 
 * CAdist is the distance of the atom from the C-alpha atom in
 * angstroms (float),  CTRdist is the distance of the atom from
 * the sidechain centroid (float). The items are separated by
 * whitespaces. Lines beginning with '#' are comments. Invalid
 * data elicit warnings and will be skipped.
 * Return value: 0 on error, 1 if OK.
 */
int Acdist_::read_file(const char *Fname)
{
    if (Fname==NULL || !strlen(Fname)) return(0);
    
    ifstream Inf(Fname);
    if (!Inf)
    {
	cerr<<"\n? Acdist_::read_file(\""<<Fname<<"\"): Cannot open\n";
	return(0);
    }
    
    Inf>>(*this);
    Inf.close();
    return(Inf.good() || Inf.eof());
}
// END of read_file()

/* >>: reads from the stream Inf. See comments to read_file() */
istream& operator>>(istream& Inf, Acdist_& Acdist)
{
    if (!Inf)
    {
	cerr<<"\n? >>Acdist_: Cannot read from stream\n";
	return(Inf);
    }
    
    static const unsigned int LINELEN=132;
    char Line[LINELEN+1];
    istrstream Instr(Line, LINELEN);
    char Aac;
    String_ Atname(10);	    // PDB atom names are 4 chars: safety margin
    float Ad, Cd;
    int Lineno, Idx;
    
    // scan all lines
    for (Lineno=1; Inf.good(); Lineno++)
    {
        Line[0]='\0'; Inf.getline(Line, LINELEN, '\n');
	if (!strlen(Line) || Line[0]=='#') continue; // skip comments and empty
	
	// read one line 
	Instr.seekg(ios::beg); Instr.clear();
	Instr>>Aac>>Atname>>Ad>>Cd;
	
	if (0>(Idx=Acdist.get_idx(Aac)))
	{
	    cerr<<"\n? >>Acdist_: Illegal amino acid code \'"
		<<Aac<<"\' in line "<<Lineno<<", skipped\n";
	    continue;
	}
	
	if (!Acdist.set_acd(Idx, Atname, Ad, Cd))	    // error msgs come from Acs_
	    cerr<<"...in line "<<Lineno<<", skipped\n";
    }
    
    return(Inf);
}
// END of >>

// ==== Acs_ METHODS ====

// ---- Access ----

/* ca_dist(), scc_dist(): return the distance of Atom from
 * the C-alpha or the sidechain centroid, respectively.
 * -1.0 returned and a warning printed if Atom was not found.
 */
float Acs_::ca_dist(const String_& Atom) const
{
    if (Atom==String_("SCC")) return(scc_dist("CA"));
    
    const Acd_ *Acdptr=get_acdptr(Atom);
    return(Acdptr==NULL? -1.0: Acdptr->Adist);
}

float Acs_::scc_dist(const String_& Atom) const
{
    if (Atom==String_("SCC")) return(0.0);

    const Acd_ *Acdptr=get_acdptr(Atom);
    return(Acdptr==NULL? -1.0: Acdptr->Cdist);
}

/* set_acd(): stores the distances of Atom from the CA and the
 * side chain centroid (Ad, Cd). Prints a warning and does
 * nothing if Atom was invalid or if Ad or Cd were negative.
 * Return value: 0 on error, non-0 otherwise.
 */
int Acs_::set_acd(const String_& Atom, float Ad, float Cd)
{
    if (Ad<0.0 || Cd<0.0)
    {
	cerr<<"\n? Acs_::set_acd(): Negative distance(s)\n";
	return(0);
    }
    
    Acd_ *Acdptr=(Acd_ *)get_acdptr(Atom);
    if (Acdptr==NULL) return(0);
    
    Acdptr->Adist=Ad;
    Acdptr->Cdist=Cd;
    return(1);
}
// END of set_acd()

// ---- Setup ----

// Modified 12-Jul-1997 when porting to GCC.

/* set_dists(): sets up the sidechain atom distances (from CA and 
 * centroid) for amino acid Aac. If Aac is lowercase, then it
 * is converted to uppercase. If Aac is illegal, then it is
 * assumed to be 'X' ("anything" or "unknown"). 'X' will store
 * nothing at all, though. In all valid cases, an internal array
 * is set up to hold the names of the sidechain atoms (according to
 * the PDB convention) and their distances. By default, these
 * distances will be taken from the most abundant rotamers (>=10%) in
 * the Ponder/Richards library as defined in Quanta 4.1.
 */
void Acs_::set_dists(char Aac)
{
    Aa=toupper(Aac);
    if (Aa<'A' || Aa>'Z') Aa='X';  // unknown
    switch(Aa)
      {
      case 'A': ala(); break;
      case 'C': cys(); break;
      case 'D': asp(); break;
      case 'E': glu(); break;
      case 'F': phe(); break;
      case 'G': gly(); break;
      case 'H': his(); break;
      case 'I': ile(); break;
      case 'K': lys(); break;
      case 'L': leu(); break;
      case 'M': met(); break;
      case 'N': asn(); break;
      case 'P': pro(); break;
      case 'Q': gln(); break;
      case 'R': arg(); break;
      case 'S': ser(); break;
      case 'T': thr(); break;
      case 'V': val(); break;
      case 'W': trp(); break;
      case 'Y': tyr(); break;
      case 'B':  // undefined or unknown amino acids
      case 'J':  // fall through...
      case 'O':
      case 'U':
      case 'X':
      case 'Z':
      default: unk();
      }
}
// END of set_dists()

// ---- Rotamer setup functions ----

/* ala(),...,unk(): These little functions set up the interatomic distance data
 * for each average rotamer. set_dists() uses them and holds
 * pointers to each in an array. 
 * unk() handles "unknown" (non-standard) amino acids. All private
 */
void Acs_::ala()
{
    Acds.len(7);
    Acds[0]=Acd_("CA", 0.00, 1.61);
    Acds[1]=Acd_("CB", 1.53, 0.08);
    Acds[2]=Acd_("HA", 1.09, 2.22);
    Acds[3]=Acd_("1HB", 2.17, 1.07);
    Acds[4]=Acd_("2HB", 2.19, 1.07);
    Acds[5]=Acd_("3HB", 2.17, 1.07);
    Acds[6]=Acd_("H", 2.13, 3.25);
}
void Acs_::cys()
{
    Acds.len(8);
    Acds[0]=Acd_("CA", 0.00, 2.36);
    Acds[1]=Acd_("CB", 1.53, 1.27);
    Acds[2]=Acd_("SG", 2.80, 0.55);
    Acds[3]=Acd_("HA", 1.09, 2.73);
    Acds[4]=Acd_("1HB", 2.15, 1.88);
    Acds[5]=Acd_("2HB", 2.17, 1.92);
    Acds[6]=Acd_("HG", 3.49, 1.59);
    Acds[7]=Acd_("H", 2.13, 3.85);
}
void Acs_::asp()
{
    Acds.len(9);
    Acds[0]=Acd_("CA", 0.00, 2.56);
    Acds[1]=Acd_("CB", 1.53, 1.50);
    Acds[2]=Acd_("CG", 2.59, 0.06);
    Acds[3]=Acd_("OD1", 2.83, 1.27);
    Acds[4]=Acd_("OD2", 3.70, 1.24);
    Acds[5]=Acd_("HA", 1.09, 2.86);
    Acds[6]=Acd_("1HB", 2.14, 2.12);
    Acds[7]=Acd_("2HB", 2.16, 2.11);
    Acds[8]=Acd_("H", 2.13, 4.08);
}
void Acs_::glu()
{
    Acds.len(12);
    Acds[0]=Acd_("CA", 0.00, 3.40);
    Acds[1]=Acd_("CB", 1.53, 2.10);
    Acds[2]=Acd_("CG", 2.69, 1.12);
    Acds[3]=Acd_("CD", 3.92, 0.54);
    Acds[4]=Acd_("OE1", 4.21, 1.33);
    Acds[5]=Acd_("OE2", 4.93, 1.72);
    Acds[6]=Acd_("HA", 1.09, 3.62);
    Acds[7]=Acd_("1HB", 2.13, 2.45);
    Acds[8]=Acd_("2HB", 2.15, 2.35);
    Acds[9]=Acd_("1HG", 2.98, 1.86);
    Acds[10]=Acd_("2HG", 3.05, 1.88);
    Acds[11]=Acd_("H", 2.13, 5.00);
}
void Acs_::phe()
{
    Acds.len(17);
    Acds[0]=Acd_("CA", 0.00, 3.71);
    Acds[1]=Acd_("CB", 1.53, 2.57);
    Acds[2]=Acd_("CG", 2.68, 1.05);
    Acds[3]=Acd_("CD1", 3.52, 1.25);
    Acds[4]=Acd_("CD2", 3.76, 1.28);
    Acds[5]=Acd_("CE1", 4.91, 1.68);
    Acds[6]=Acd_("CE2", 5.08, 1.70);
    Acds[7]=Acd_("CZ", 5.55, 1.86);
    Acds[8]=Acd_("HA", 1.09, 4.04);
    Acds[9]=Acd_("1HB", 2.11, 2.99);
    Acds[10]=Acd_("2HB", 2.16, 3.03);
    Acds[11]=Acd_("HD1", 3.45, 2.33);
    Acds[12]=Acd_("HD2", 3.89, 2.35);
    Acds[13]=Acd_("HE1", 5.68, 2.76);
    Acds[14]=Acd_("HE2", 5.94, 2.78);
    Acds[15]=Acd_("HZ", 6.64, 2.95);
    Acds[16]=Acd_("H", 2.13, 5.03);
}
void Acs_::gly()
{
    Acds.len(4);
    Acds[0]=Acd_("CA", 0.00, 1.09);
    Acds[1]=Acd_("HA", 1.09, 1.79);
    Acds[2]=Acd_("2HA", 1.09, 0.00);
    Acds[3]=Acd_("H", 2.13, 2.84);
}
void Acs_::his()
{
    Acds.len(14);
    Acds[0]=Acd_("1HB", 2.14, 2.80);
    Acds[1]=Acd_("2HB", 2.15, 2.81);
    Acds[2]=Acd_("CA", 0.00, 3.16);
    Acds[3]=Acd_("CB", 1.53, 2.24);
    Acds[4]=Acd_("CG", 2.59, 0.69);
    Acds[5]=Acd_("CD2", 3.54, 1.12);
    Acds[6]=Acd_("CE1", 4.37, 1.44);
    Acds[7]=Acd_("ND1", 3.34, 1.08);
    Acds[8]=Acd_("NE2", 4.54, 1.61);
    Acds[9]=Acd_("HA", 1.09, 3.38);
    Acds[10]=Acd_("HD1", 3.59, 2.02);
    Acds[11]=Acd_("HD2", 3.91, 2.18);
    Acds[12]=Acd_("HE1", 5.29, 2.52);
    Acds[13]=Acd_("H", 2.13, 4.52);
}
void Acs_::ile()
{
    Acds.len(16);
    Acds[0]=Acd_("CA", 0.00, 2.47);
    Acds[1]=Acd_("CB", 1.53, 1.08);
    Acds[2]=Acd_("CG1", 2.66, 1.00);
    Acds[3]=Acd_("CG2", 2.63, 1.82);
    Acds[4]=Acd_("CD1", 3.97, 1.79);
    Acds[5]=Acd_("HA", 1.09, 2.88);
    Acds[6]=Acd_("HB", 2.06, 1.71);
    Acds[7]=Acd_("1HG1", 2.96, 1.75);
    Acds[8]=Acd_("1HG2", 3.58, 2.01);
    Acds[9]=Acd_("2HG1", 2.95, 1.91);
    Acds[10]=Acd_("2HG2", 2.89, 2.76);
    Acds[11]=Acd_("3HG2", 2.92, 2.34);
    Acds[12]=Acd_("1HD1", 4.75, 2.73);
    Acds[13]=Acd_("2HD1", 4.27, 2.21);
    Acds[14]=Acd_("3HD1", 4.27, 2.09);
    Acds[15]=Acd_("H", 2.13, 4.03);
}
void Acs_::lys()
{
    Acds.len(19);
    Acds[0]=Acd_("CA", 0.00, 4.06);
    Acds[1]=Acd_("CB", 1.53, 2.82);
    Acds[2]=Acd_("CG", 2.71, 1.53);
    Acds[3]=Acd_("CD", 4.01, 0.31);
    Acds[4]=Acd_("CE", 5.32, 1.38);
    Acds[5]=Acd_("NZ", 6.49, 2.48);
    Acds[6]=Acd_("HA", 1.09, 4.19);
    Acds[7]=Acd_("1HB", 2.11, 3.21);
    Acds[8]=Acd_("2HB", 2.13, 3.10);
    Acds[9]=Acd_("1HG", 2.92, 2.06);
    Acds[10]=Acd_("2HG", 3.05, 2.07);
    Acds[11]=Acd_("1HD", 4.29, 1.31);
    Acds[12]=Acd_("2HD", 4.15, 1.29);
    Acds[13]=Acd_("1HE", 5.47, 1.97);
    Acds[14]=Acd_("2HE", 5.57, 1.98);
    Acds[15]=Acd_("1HZ", 7.33, 3.32);
    Acds[16]=Acd_("2HZ", 6.67, 2.78);
    Acds[17]=Acd_("3HZ", 6.58, 2.77);
    Acds[18]=Acd_("H", 2.13, 5.57);
}
void Acs_::leu()
{
    Acds.len(16);
    Acds[0]=Acd_("CA", 0.00, 2.82);
    Acds[1]=Acd_("CB", 1.53, 1.62);
    Acds[2]=Acd_("CG", 2.76, 0.32);
    Acds[3]=Acd_("CD1", 3.96, 1.53);
    Acds[4]=Acd_("CD2", 3.52, 1.53);
    Acds[5]=Acd_("HA", 1.09, 3.01);
    Acds[6]=Acd_("1HB", 2.10, 2.23);
    Acds[7]=Acd_("2HB", 2.11, 2.12);
    Acds[8]=Acd_("HG", 2.96, 1.41);
    Acds[9]=Acd_("1HD1", 4.84, 2.24);
    Acds[10]=Acd_("2HD1", 4.10, 2.28);
    Acds[11]=Acd_("3HD1", 4.29, 2.04);
    Acds[12]=Acd_("1HD2", 4.49, 2.23);
    Acds[13]=Acd_("2HD2", 3.84, 2.04);
    Acds[14]=Acd_("3HD2", 3.36, 2.26);
    Acds[15]=Acd_("H", 2.13, 4.44);
}
void Acs_::met()
{
    Acds.len(14);
    Acds[0]=Acd_("CA", 0.00, 3.26);
    Acds[1]=Acd_("CB", 1.53, 2.06);
    Acds[2]=Acd_("CG", 2.69, 1.24);
    Acds[3]=Acd_("SD", 3.88, 0.85);
    Acds[4]=Acd_("CE", 4.68, 1.79);
    Acds[5]=Acd_("HA", 1.09, 3.36);
    Acds[6]=Acd_("1HB", 2.12, 2.52);
    Acds[7]=Acd_("2HB", 2.15, 2.25);
    Acds[8]=Acd_("1HG", 2.92, 2.15);
    Acds[9]=Acd_("2HG", 3.27, 1.81);
    Acds[10]=Acd_("1HE", 5.55, 2.65);
    Acds[11]=Acd_("2HE", 5.13, 2.32);
    Acds[12]=Acd_("3HE", 4.30, 2.09);
    Acds[13]=Acd_("H", 2.13, 4.90);
}
void Acs_::asn()
{
    Acds.len(11);
    Acds[0]=Acd_("CA", 0.00, 2.57);
    Acds[1]=Acd_("CB", 1.53, 1.53);
    Acds[2]=Acd_("CG", 2.57, 0.05);
    Acds[3]=Acd_("ND2", 3.54, 1.30);
    Acds[4]=Acd_("OD1", 3.05, 1.28);
    Acds[5]=Acd_("HA", 1.09, 2.86);
    Acds[6]=Acd_("1HB", 2.15, 2.12);
    Acds[7]=Acd_("2HB", 2.18, 2.13);
    Acds[8]=Acd_("2HD2", 3.77, 1.98);
    Acds[9]=Acd_("1HD2", 4.29, 2.00);
    Acds[10]=Acd_("H", 2.13, 4.12);
}
void Acs_::pro()
{
    Acds.len(11);
    Acds[0]=Acd_("CA", 0.00, 1.97);
    Acds[1]=Acd_("CB", 1.58, 1.28);
    Acds[2]=Acd_("CG", 2.49, 0.59);
    Acds[3]=Acd_("CD", 2.46, 1.28);
    Acds[4]=Acd_("HA", 1.09, 2.69);
    Acds[5]=Acd_("1HB", 2.24, 2.14);
    Acds[6]=Acd_("2HB", 2.22, 1.94);
    Acds[7]=Acd_("1HG", 3.03, 1.54);
    Acds[8]=Acd_("2HG", 3.39, 1.50);
    Acds[9]=Acd_("1HD", 3.35, 2.13);
    Acds[10]=Acd_("2HD", 3.03, 1.94);
}
void Acs_::gln()
{
    Acds.len(14);
    Acds[0]=Acd_("CA", 0.00, 3.39);
    Acds[1]=Acd_("CB", 1.53, 2.08);
    Acds[2]=Acd_("CG", 2.70, 1.16);
    Acds[3]=Acd_("CD", 3.89, 0.53);
    Acds[4]=Acd_("OE1", 4.80, 1.67);
    Acds[5]=Acd_("NE2", 4.32, 1.44);
    Acds[6]=Acd_("HA", 1.09, 3.62);
    Acds[7]=Acd_("1HB", 2.12, 2.47);
    Acds[8]=Acd_("2HB", 2.14, 2.33);
    Acds[9]=Acd_("1HG", 2.96, 1.91);
    Acds[10]=Acd_("2HG", 3.11, 1.90);
    Acds[11]=Acd_("1HE2", 5.23, 2.26);
    Acds[12]=Acd_("2HE2", 3.99, 1.90);
    Acds[13]=Acd_("H", 2.13, 4.91);
}
void Acs_::arg()
{
    Acds.len(21);
    Acds[0]=Acd_("CA", 0.00, 4.95);
    Acds[1]=Acd_("CB", 1.53, 3.67);
    Acds[2]=Acd_("CG", 2.67, 2.29);
    Acds[3]=Acd_("CD", 4.10, 1.03);
    Acds[4]=Acd_("NE", 5.21, 0.53);
    Acds[5]=Acd_("CZ", 6.52, 1.59);
    Acds[6]=Acd_("NH1", 7.35, 2.61);
    Acds[7]=Acd_("NH2", 7.20, 2.42);
    Acds[8]=Acd_("HA", 1.09, 5.19);
    Acds[9]=Acd_("1HB", 2.13, 3.90);
    Acds[10]=Acd_("2HB", 2.13, 3.93);
    Acds[11]=Acd_("1HG", 2.90, 2.64);
    Acds[12]=Acd_("2HG", 2.89, 2.61);
    Acds[13]=Acd_("1HD", 4.39, 1.64);
    Acds[14]=Acd_("2HD", 4.41, 1.69);
    Acds[15]=Acd_("HE", 5.09, 1.41);
    Acds[16]=Acd_("1HH1", 8.32, 3.51);
    Acds[17]=Acd_("2HH1", 7.14, 2.84);
    Acds[18]=Acd_("1HH2", 8.19, 3.37);
    Acds[19]=Acd_("2HH2", 6.84, 2.50);
    Acds[20]=Acd_("H", 2.13, 6.32);
}
void Acs_::ser()
{
    Acds.len(8);
    Acds[0]=Acd_("CA", 0.00, 2.01);
    Acds[1]=Acd_("CB", 1.53, 0.78);
    Acds[2]=Acd_("OG", 2.46, 0.66);
    Acds[3]=Acd_("HA", 1.09, 2.60);
    Acds[4]=Acd_("1HB", 2.14, 1.50);
    Acds[5]=Acd_("2HB", 2.17, 1.52);
    Acds[6]=Acd_("HG", 3.30, 1.30);
    Acds[7]=Acd_("H", 2.13, 3.49);
}
void Acs_::thr()
{
    Acds.len(11);
    Acds[0]=Acd_("CA", 0.00, 2.02);
    Acds[1]=Acd_("CB", 1.53, 0.62);
    Acds[2]=Acd_("OG1", 2.44, 1.19);
    Acds[3]=Acd_("CG2", 2.67, 1.32);
    Acds[4]=Acd_("HA", 1.09, 2.55);
    Acds[5]=Acd_("HB", 2.10, 1.49);
    Acds[6]=Acd_("HG1", 2.40, 1.78);
    Acds[7]=Acd_("1HG2", 3.57, 1.78);
    Acds[8]=Acd_("2HG2", 2.95, 2.22);
    Acds[9]=Acd_("3HG2", 2.94, 1.86);
    Acds[10]=Acd_("H", 2.13, 3.55);
}
void Acs_::val()
{
    Acds.len(13);
    Acds[0]=Acd_("CA", 0.00, 2.07);
    Acds[1]=Acd_("CB", 1.53, 0.65);
    Acds[2]=Acd_("CG1", 2.65, 1.35);
    Acds[3]=Acd_("CG2", 2.63, 1.35);
    Acds[4]=Acd_("HA", 1.09, 2.54);
    Acds[5]=Acd_("HB", 2.10, 1.48);
    Acds[6]=Acd_("1HG1", 3.59, 1.82);
    Acds[7]=Acd_("2HG1", 2.93, 2.25);
    Acds[8]=Acd_("3HG1", 2.95, 1.94);
    Acds[9]=Acd_("1HG2", 3.58, 1.82);
    Acds[10]=Acd_("2HG2", 2.91, 1.94);
    Acds[11]=Acd_("3HG2", 2.91, 2.24);
    Acds[12]=Acd_("H", 2.13, 3.69);
}
void Acs_::trp()
{
    Acds.len(21);
    Acds[0]=Acd_("CA", 0.00, 3.89);
    Acds[1]=Acd_("CB", 1.53, 3.06);
    Acds[2]=Acd_("CG", 2.59, 1.65);
    Acds[3]=Acd_("CD1", 3.56, 2.12);
    Acds[4]=Acd_("CD2", 3.49, 0.49);
    Acds[5]=Acd_("NE1", 4.66, 1.93);
    Acds[6]=Acd_("CE2", 4.60, 0.92);
    Acds[7]=Acd_("CE3", 3.96, 1.69);
    Acds[8]=Acd_("CZ2", 5.83, 2.03);
    Acds[9]=Acd_("CZ3", 5.33, 2.45);
    Acds[10]=Acd_("CH2", 6.13, 2.58);
    Acds[11]=Acd_("HA", 1.09, 4.01);
    Acds[12]=Acd_("1HB", 2.14, 3.59);
    Acds[13]=Acd_("2HB", 2.16, 3.52);
    Acds[14]=Acd_("HD1", 3.86, 3.22);
    Acds[15]=Acd_("HE1", 5.52, 2.81);
    Acds[16]=Acd_("HE3", 3.64, 2.59);
    Acds[17]=Acd_("HH2", 7.18, 3.66);
    Acds[18]=Acd_("HZ2", 6.72, 2.98);
    Acds[19]=Acd_("HZ3", 5.96, 3.51);
    Acds[20]=Acd_("H", 2.13, 5.33);
}
void Acs_::tyr()
{
    Acds.len(18);
    Acds[0]=Acd_("CA", 0.00, 4.09);
    Acds[1]=Acd_("CB", 1.53, 3.09);
    Acds[2]=Acd_("CG", 2.70, 1.52);
    Acds[3]=Acd_("CD1", 3.53, 1.43);
    Acds[4]=Acd_("CD2", 3.69, 1.45);
    Acds[5]=Acd_("CE1", 4.84, 1.40);
    Acds[6]=Acd_("CE2", 4.96, 1.42);
    Acds[7]=Acd_("CZ", 5.42, 1.38);
    Acds[8]=Acd_("OH", 6.72, 2.76);
    Acds[9]=Acd_("HA", 1.09, 4.28);
    Acds[10]=Acd_("1HB", 2.10, 3.56);
    Acds[11]=Acd_("2HB", 2.16, 3.58);
    Acds[12]=Acd_("HD1", 3.55, 2.53);
    Acds[13]=Acd_("HE1", 5.63, 2.49);
    Acds[14]=Acd_("HD2", 3.85, 2.54);
    Acds[15]=Acd_("HE2", 5.78, 2.51);
    Acds[16]=Acd_("HH", 7.27, 3.24);
    Acds[17]=Acd_("H", 2.13, 5.42);
}
void Acs_::unk()
{
    cerr<<"\n? Acs_::set_dists(X): Unknown amino acid, prepare for trouble...\n";
    Acds.len(0);
}

// ==== END OF METHODS Acdist.c++ ====
