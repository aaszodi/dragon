// ==== PROJECT DRAGON: METHODS Restr.c++ ====

/* Stores, smooths and retrieves distance restraints
 * in a uniform way.
 */

// SGI C++ 7.1, IRIX 6.2, 21. May 1998. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <fstream.h>
#include <strstream.h>
#include <string.h>
#include <ctype.h>

// ---- MODULE HEADERS ----

#include "Restr.h"
#include "portrandom.h"

// ---- DEFINITIONS ----

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif

// ==== Restr_ METHODS ====

// ---- Input,output ----

/* >>: reads from the stream In assuming the following format:-
 * "Pos1 Pos2 Lowlim Uplim Strict Atom1 Atom2", where the first five things
 * are the Restr_ data members. Atom1 and Atom2 are atom names corresponding
 * to the PDB conventions,  such as "CA" or "1HB2". Lowercase strings are
 * converted to uppercase but their validity is not checked here. 
 * The string "SCC" may also be specified: it means the side-chain
 * centroid. The maximal string length is 4 characters.
 * Residue numbers are [1..Rno] but the upper limit
 * is not checked here. 
 * 
 * NOTE: In earlier versions (3.x and <=4.7) distance restraints
 * could be specified only between C-alpha and pseudo-C-beta (SCC) atoms.
 * The format was:-
 * "Pos1 Pos2 Lowlim Uplim Strict ABstr", where ABstr corresponded to
 * the 2-char regular expression "[aAbB][aAbB]" with 'a' or 'A' indicating
 * a C-alpha atom, 'b' or 'B' a dummy "C-beta" (SCC). To maintain
 * backward compatibility, this format is also recognised.
 */
istream& operator>>(istream& In, Restr_& R)
{
    static const int ATNAMELEN=4;   // max. no. of chars in atom name
    
    int P1, P2;
    double L, U, S;
    String_ A1(ATNAMELEN), A2(ATNAMELEN);
    
    if (!In)
    {
	cerr<<"\n? ...>>Restr_: Cannot read from stream\n";
	return(In);
    }
    
    // read positions
    In>>P1>>P2;
    if (!In)
    {
	cerr<<"\n? ...>>Restr_: Read error at residue numbers\n";
	return(In);
    }
    
    if (P1<=0 || P2<=0)
    {
	cerr<<"\n? ...>>Restr_: Residue numbers must be >0\n";
	In.clear(ios::failbit|In.rdstate());	// set error flag
	return(In);
    }
    
    if (P1==P2)
    {
	cerr<<"\n? ...>>Restr_: Restraint between 1 residue not allowed\n";
	In.clear(ios::failbit|In.rdstate());	// set error flag
	return(In);
    }
    
    // read distance limits
    In>>L>>U;
    if (!In)
    {
	cerr<<"\n? ...>>Restr_: Read error at distance limits\n";
	return(In);
    }
    
    if (L<0.0 || U<0.0)
    {
	cerr<<"\n? ...>>Restr_: Negative distance limit(s) not allowed\n";
	In.clear(ios::failbit|In.rdstate());	// set error flag
	return(In);
    }
    if (L>U) { double T=L; L=U; U=T; }	// swap silently
    
    // read strictness
    In>>S;
    if (!In)
    {
	cerr<<"\n? ...>>Restr_: Read error at strictness\n";
	return(In);
    }
    
    // read atom specs
    
    /* NOTE: the >>(char*) does not handle EOLN conditions correctly
     * and sometimes reads garbage. Doesn't seem to be my fault.
     * Purify detects UMR conditions deep in "istream". Theoretically, 
     * >> should put a '\0' into (char*) if extraction didn't work:
     * this does not happen in this implementation.
     */
    In>>A1;
    if (!In || !strlen(A1))
    {
	cerr<<"\n? ...>>Restr_: Read error at first atom specification string\n";
	In.clear(ios::failbit|In.rdstate());	// set error flag
	return(In);
    }
    A1.toupper();   // uppercase

    /* Check if A1 is the old format */
    if (strlen(A1)==2 && 
	strchr("AB", A1[0])!=NULL && strchr("AB", A1[1])!=NULL)
    {
	A2=(A1[1]=='A')? "CA": "SCC";	// fake second atom name
	A1=(A1[0]=='A')? "CA": "SCC";	// and then the first
    }
    else
    {
	In>>A2;
	if (!In || !strlen(A2))
	{
	    cerr<<"\n? ...>>Restr_: Read error at second atom specification string\n";
	    In.clear(ios::failbit|In.rdstate());	// set error flag
	    return(In);
	}
	A2.toupper();   // uppercase
    }
    
    // seems more or less OK, copy values
    R.pos(1, P1); R.pos(2, P2);	// keep [1..] range
    R.low(L); R.up(U);
    R.atom(1, A1); R.atom(2, A2);
    R.strict(S);	    // was missing until 24-Nov-95...
    return(In);
}
// END of >>Restr_

/* Output: lists to Out corresponding to the input format. */
ostream& operator<<(ostream& Out, const Restr_& R)
{
    Out<<R.pos(1)<<" "<<R.pos(2)<<" "
	<<R.low()<<" "<<R.up()<<" "<<R.strict()<<" "
	<<R.atom(1)<<" "<<R.atom(2)<<endl;
    return(Out);
}
// END of <<Restr_

// ==== METHODS ====

// ---- Static init ----

const float Restraints_::CA_BONDANGLE=2.33F;   // max. CA:CA:CA virtual bond angle (radians)
const float Restraints_::CA_BUMP=2.46F;	// C-alpha bump radius
const float Restraints_::CA_1_MIN=3.75F;	// CA:CA virtual bond length
const float Restraints_::CA_1_MAX=3.85F;
const float Restraints_::CA_2_MIN=6.00F;	// second neighbour distance
const float Restraints_::CA_2_MAX=7.00F;

const float Restraints_::NTCA_BUMP=3.95F; // 1.49+CA_BUMP
const float Restraints_::NT_BONDLEN=1.47F;  // H3N-C bond length
const float Restraints_::NT_2_MIN=5.0F;    // min. Nt:second CA
const float Restraints_::NT_2_MAX=5.8F;    // max. Nt:second CA
const float Restraints_::CTCA_BUMP=4.46F; // 2.0+CA_BUMP
const float Restraints_::CT_BONDLEN=1.54F;  // C-COO bond length
const float Restraints_::CT_2_MIN=4.9F;    // min. Ct:second CA
const float Restraints_::CT_2_MAX=5.9F;    // max. Ct:second CA
const float Restraints_::NTCT_BUMP=3.49F; // NT_BUMP+CT_BUMP

const float Restraints_::STR1=2.0F;	// CA:CA neighbour strictness
const float Restraints_::STR2=1.5F;	// CA:CA second neighbour strictness
const float Restraints_::STRA=1.0F;	// general CA:CA vdW bump strictness
const float Restraints_::STRB=0.7F;	// strictness for bumps involving sidechains

// ---- Access ----

/* check_index(): returns 0 and generates a warning if either of
 * its arguments are larger than the current size. Otherwise returns
 * 1 but swaps its arguments if i<j so that the access methods will
 * always know that i>=j after index checking. This is needed because
 * the upper limits are stored in the upper triangles.
 * Private
 */
int Restraints_::check_index(unsigned int& i, unsigned int& j) const
{
    if (i>=Size)
    {
	cerr<<"\n? Restraints_::check_index("<<i<<", ...) out of range [0.."<<Size-1<<"]\n";
	return(0);
    }
    if (j>=Size)
    {
	cerr<<"\n? Restraints_::check_index(..., "<<j<<") out of range [0.."<<Size-1<<"]\n";
	return(0);
    }
    if (i<j) { register unsigned int t=i; i=j; j=t; }
    return(1);
}
// END of check_index()

/* set_size(): adjusts the size for an Rno-long model chain (that
 * is actually made up of Rno+2 points). Zeroes all internal data
 * if the size has changed.
 * Returns old size.
 */
unsigned int Restraints_::set_size(unsigned int Rno)
{
    Rno+=2;
    if (Size==Rno) return(Size);    // do nothing
    
    unsigned int Oldsize=Size; Size=Rno;
    Lowup.set_size(Size); Lowup.set_values();
    Lowup2.set_size(Size); Lowup2.set_values();
    Strict.set_size(Size); Strict.set_values();
    Maxsepar.len(Size); Maxsepar.set_values(0.0);
    return(Oldsize);
}
// END of set_size()

/* low(),low2(),up(),up2(): retrieve or set the unsquared or squared
 * lower or upper limits between the ith and jth C-alphas.
 * Calling low(i, j, L) is equivalent to a call
 * to low2(i, j, L*L). In general, i, j is in the range [1..Rno]
 * and j may be larger or smaller than i. i==j returns 0.0.
 * Negative distance limits are silently converted to their abs values.
 * On index range error, no "set" operations are done and the "get"
 * methods return 0.0 with a warning.
 */
double Restraints_::low(unsigned int i, unsigned int j) const
{
    return(check_index(i, j)? Lowup[i][j]: 0.0);
}
void Restraints_::low(unsigned int i, unsigned int j, double Low)
{
    if (!check_index(i, j)) return;
    Lowup[i][j]=fabs(Low); Lowup2[i][j]=Low*Low;
}
double Restraints_::low2(unsigned int i, unsigned int j) const
{
    return(check_index(i, j)? Lowup2[i][j]: 0.0);
}
void Restraints_::low2(unsigned int i, unsigned int j, double Low2)
{
    if (!check_index(i, j)) return;
    Low2=fabs(Low2); Lowup2[i][j]=Low2; Lowup[i][j]=sqrt(Low2);
}
double Restraints_::up(unsigned int i, unsigned int j) const
{
    return(check_index(i, j)? Lowup[j][i]: 0.0);
}
void Restraints_::up(unsigned int i, unsigned int j, double Up)
{
    if (!check_index(i, j)) return;
    Lowup[j][i]=fabs(Up); Lowup2[j][i]=Up*Up;
}
double Restraints_::up2(unsigned int i, unsigned int j) const
{
    return(check_index(i, j)? Lowup2[j][i]: 0.0);
}
void Restraints_::up2(unsigned int i, unsigned int j, double Up2)
{
    if (!check_index(i, j)) return;
    Up2=fabs(Up2); Lowup2[j][i]=Up2; Lowup[j][i]=sqrt(Up2);
}
// END of low()...up2()

/* strict(): gets or sets the strictness value for the i:j C-alpha restraint.
 * Negative strictnesses are silently converted to absolute values.
 * Returns the correct values for the non-specific restraints
 * (bonds and bumps).
 */
double Restraints_::strict(unsigned int i, unsigned int j) const
{
    if (!check_index(i, j)) return(0.0);	// nonsense
    if (i-j==1) return(STR1);   // CA:CA bond
    if (i-j==2) return(STR2);	// geminal CA distance
    return(Strict[i][j]>0.0? Strict[i][j]: STRA);   //spec or unspec restr
}
// END of strict()

// ---- Restraint setup ----

/* merge_restr(): inspects the restraint between "i" and "j"
 * and updates it with the UNsquared values Low, Up if these
 * define a narrower range. If the original restraint was uninitialised
 * (ie. [0, 0]) then it is overwritten. The strictness Str overwrites
 * the original strictness (it is possible
 * to narrow down a stricter restraint by a less strict one).
 * Return value: The number of limits accepted (0, 1, 2).
 * Private
 */
inline
int Restraints_::merge_restr(unsigned int i, unsigned int j, 
	double Low, double Up, double Str)
{
    if (Low>Up) return(0);
    
    register double Olow=low(i, j), Oup=up(i, j);
    int Modify=0;
    
    if (Low<Oup && (Olow==0.0 || Olow<=Low)) { Modify++; low(i, j, Low); }
    if (Up>Olow && (Oup==0.0 || Up<=Oup)) { Modify++; up(i, j, Up); }
    if (Modify) strict(i, j, Str);
    return(Modify);
}
// END of merge_restr()

/* setup_restr(): clears the restraint matrices and sets them up
 * according to the list of external restraints (which should already
 * be prepared in the calling object), secondary structure (from Pieces)
 * and intra-monomer atom distances (from Polymer). Call only once
 * before the simulations.
 */
void Restraints_::setup_restr(const Pieces_& Pieces, const Polymer_& Polymer)
{
    // zero everything
    Lowup.set_values(); Lowup2.set_values();
    Strict.set_values(); Maxsepar.set_values(0.0);
    flory_constr();
    
    // set bond/bump, secstr, external
    setup_bondbump();
    setup_secstrestr(Pieces);
    setup_extrestr(Polymer);

    // smooth restraints
    smooth_restr(0);
}
// END of setup_restr()

/* setup_bondbump(): adds the bond lengths, C-alpha van der Waals
 * bump limits and the Flory maximal separation upper limit to each
 * C-alpha pair. Assumes that the strictnesses were 0.0 and doesn't set them.
 * Private
 */
void Restraints_::setup_bondbump()
{
    register unsigned int d, i, j;
    
    // set the C-alpha virtual bonds (overwrites everybody)
    for (i=2; i<Size-1; i++)
    {
	j=i-1;
	low(i, j, CA_1_MIN); up(i, j, CA_1_MAX);
    }

    // set the 2nd neighbour limits (overwrites everybody)
    for (i=3; i<Size-1; i++)
    {
	j=i-2;
	low(i, j, CA_2_MIN); up(i, j, CA_2_MAX);
    }
	
    // set the more distant neighbours
    for (d=3; d<Size-1; d++)
	for (i=d; i<Size-1; i++)
	{
	    j=i-d;
	    low(i, j, 2.0*CA_BUMP);
	    up(i, j, max_separ(d));
	}
	
    // set the limits between the N-terminus and the others
    low(1, 0, NT_BONDLEN); up(1, 0, NT_BONDLEN);
    low(2, 0, NT_2_MIN); up(2, 0, NT_2_MAX);
    for (i=3; i<Size-1; i++)
    {
	low(i, 0, NTCA_BUMP); up(i, 0, max_separ(i));
    }
    
    // set the limits between the C-terminus and the others
    for (i=1; i<Size-3; i++)
    {
	low(Size-1, i, CTCA_BUMP); up(Size-1, i, max_separ(Size-i-1));
    }

    low(Size-1, Size-3, CT_2_MIN); up(Size-1, Size-3, CT_2_MAX);
    low(Size-1, Size-2, CT_BONDLEN); up(Size-1, Size-2, CT_BONDLEN);
    
    // limits between the N- and C-terminus
    low(0, Size-1, NTCT_BUMP); up(0, Size-1, max_separ(Size-1));
}
// END of setup_bondbump()

/* setup_extrestr(): adds the external restraints defined in the
 * internal restraint list Restrs to the restraint matrices. 
 * CA:CA restraints will be deleted from Restrs and will be adjusted
 * together with the secstr and bonds/bumps etc. Non-CA:CA
 * restraints will be converted to CA:CA for smoothing purposes
 * but with a 0.0 weight.
 * Private
 */
void Restraints_::setup_extrestr(const Polymer_& Polymer)
{
    static const String_ CA("CA");  // string comparison
    static const float CA_MINDIST=2.0*CA_BUMP;
    
    register unsigned i, j;
    double L, U, D1, D2;
    String_ A1, A2;
    
    for (Restrs.begin(); Restrs!=NULL; Restrs++)
    {
	i=Restrs->pos(1); j=Restrs->pos(2);
	if (abs(i-j)<2)	// no CA:CA "bond" restraints
	{
	    Restrs.del();
	    continue;
	}
	
	L=Restrs->low(); U=Restrs->up();
	A1=Restrs->atom(1); A2=Restrs->atom(2);
	if (A1==CA && A2==CA)	// CA:CA restraint, put into matrices
	{
	    merge_restr(i, j, L, U, Restrs->strict());
	    Restrs.del();   // remove from list
	    continue;
	}
	
	// modify the restraint with the CA:Ax atom distances
	D1=(A1==CA)? 0.0: Polymer.ca_dist(i-1, A1);
	D2=(A2==CA)? 0.0: Polymer.ca_dist(j-1, A2);
	L-=(D1+D2); U+=(D1+D2);
	if (L<CA_MINDIST) L=CA_MINDIST;	// keep hard vdW

	merge_restr(i, j, L, U, 0.5*Restrs->strict());
    }
}
// END of setup_extrestr()

/* setup_secstrestr(): adds the distance restraints defined by the 
 * secondary structure. These overwrite everything they see.
 * The list of secondary structures comes from Pieces.
 * Private
 */
void Restraints_::setup_secstrestr(const Pieces_& Pieces)
{
    // Factors by which the distances are allowed to deviate
    const float LODEVFACT=0.99F, HIDEVFACT=1.01F;
    
    Clist1_<Sstr_> Slist=Pieces.secs();    // secondary struct list iterator
    if (!Slist.len()) return;	    // there were none
    
    Trimat_ Idist(Size), Strimat(Size);	// temp storage
    Idist.set_values(); Strimat.set_values();
    
    // copy the ideal UNsquared distances and the strictnesses
    for (Slist.begin(); Slist!=NULL; Slist++)
	(*Slist)->ideal_dist(Idist, Strimat);

    // update the restraints
    register unsigned int i, j;
    register double D, S;
    
    for (i=2; i<Size-1; i++)
	for (j=1; j<i; j++)
	{
	    if ((S=Strimat[i][j])<=0.0) continue;   // no modification
	    D=Idist[i][j];
	    low(i, j, D*LODEVFACT); up(i, j, D*HIDEVFACT);
	    Strict[i][j]=S;
	}
}
// END of setup_secstrestr()

/* smooth_restr(): tries to narrow the ranges of those restraints which
 * have a strictness of 0.0 ("soft"). Assumes that all "hard" restraints
 * are already set up. Performs Pass passes (or as many as necessary
 * if Pass=0, default 1) and does not check the
 * incompatibility between the "hard" restraints.
 * Return value: the number of inequality violations found.
 * Private
 */
int Restraints_::smooth_restr(unsigned int Pass)
{
    // NOTE: this is a real Piglet Algorithm. Use with caution.
    
    // NOTE: this EPSILON was needed by GCC when compiling with -O2
    static const double EPSILON=FLT_EPSILON;
    
    Trimat_ Newbound(Size);    // temp storage for new bounds
    register unsigned int i, j, k;
    register double Lnew, Unew, Ltemp, Btemp;
    int Cyc, Adjno, Violno;
    
    // fix upper limits first (following Kuntz...)
    cout<<"SMUP: "<<flush;
    for (Cyc=0; !Pass || Cyc<Pass; Cyc++)
    {
	Newbound.set_values(0.0); Adjno=0;
    	cout<<'.'<<flush;

	// scan all pairs, collect new unsquared upper limits in Newbound
	for (i=0; i<Size; i++)
	{
	    for (j=0; j<i; j++)
	    {
		if (hard(i, j)) continue;	// was "hard" restraint, no change

		Btemp=up(i, j);
		for (k=0; k<Size; k++)
		{
		    if (k==i || k==j) continue;
		    Unew=up(i, k)+up(j, k);
		    
		    // upper limit can be narrowed
		    if (Btemp>Unew+EPSILON)
		    {
			Btemp=Unew; ++Adjno;
		    }
		}	// for k
		if (Btemp<up(i, j)-EPSILON) Newbound[i][j]=Btemp;
		
	    }	    // for j
	}	// for i
	if (!Adjno) break;  // no further modifications
	
	// copy modified upper limits back
	for (i=0; i<Size; i++)
	    for (j=0; j<i; j++)
		if (Newbound[i][j]>0.0) up(i, j, Newbound[i][j]);
    }	    // for Cyc
    cout<<Cyc<<endl;
    
    // fix lower limits
    cout<<"SMLOW: "<<flush;
    for (Cyc=0; !Pass || Cyc<Pass; Cyc++)
    {
	Newbound.set_values(0.0); Adjno=Violno=0;
	cout<<'.'<<flush;
	// scan all pairs, collect new lower limits in Newbound
	for (i=0; i<Size; i++)
	{
	    for (j=0; j<i; j++)
	    {
		if (hard(i, j)) continue;	// was "hard" restraint
		
		Btemp=low(i, j);
		for (k=0; k<Size; k++)
		{
		    if (k==i || k==j) continue;
		    Lnew=low(i, k)-up(j, k);
		    Ltemp=low(j, k)-up(i, k);
		    if (Ltemp>Lnew) Lnew=Ltemp;
		    if (up(i, j)<Lnew)
		    {
			Violno++; continue; // inequality violation :-(
		    }
		    
		    // lower limit can be narrowed
		    if (Btemp<Lnew-EPSILON) { Btemp=Lnew; ++Adjno; }
		}	// for k
		if (Btemp>low(i, j)+EPSILON) Newbound[i][j]=Btemp;
		
	    }	    // for j
	}	// for i
	if (!Adjno) break;  // no further modifications
	
	// copy modified lower limits back
	for (i=0; i<Size; i++)
	    for (j=0; j<i; j++)
		if (Newbound[i][j]>0.0) low(i, j, Newbound[i][j]);
	
	if (Violno) break;  // something is fishy, stop here
    }	    // for Cyc
    cout<<Cyc<<", triangle violations="<<Violno<<endl;
    return(Violno);
}
// END of smooth_restr()

// ---- Distance matrix initialisation ----

/* init_distmat(): produces a random (squared) distance matrix with entries from a
 * Gaussian distribution. The average and S.D is calculated from the
 * expected radius of the molecule, assuming a spherical shape.
 * The random number generator is initialised with
 * Randseed. Only random numbers falling between the distance bounds
 * are accepted. "Softly" restrained positions are modified by the
 * hydrophobic distance estimates.
 */
void Restraints_::init_distmat(Trimat_& Dist, Polymer_& Polymer, 
	long Randseed) const
{
    register unsigned int i, j, Ptno=Dist.rno();
    
    if (Ptno!=Size)
    {
	cerr<<"\n? Restraints_::init_distmat(): Size mismatch (adjusted)\n";
	Dist.set_size(Ptno=Size);
    }
    
    register double Drand, Low, Up, Destim;
    double Rexp, Avgdist, Dev, Strict;
    
    Rexp=exp_rad(Ptno-2);
    Avgdist=36.0*Rexp/35.0; Dev=sqrt(1.2)*Rexp;

    init_portrand(Randseed);     /* init RNG */
    for (i=0; i<Ptno; i++)
    {
	for (j=0; j<=i; j++)
	{
	    // zero the main diagonal for safety's sake
	    if (i==j)
	    {
		Dist[i][j]=0.0;
		continue;
	    }
	    
	    // fill "fixed" distances with the desired value
	    Low=low(i, j); Up=up(i, j);
	    if (Low==Up)
	    {
		Dist[i][j]=low2(i, j);	// squared
		continue;
	    }
	    
	    /* Try to generate Gaussian random values within the limits
	     * with average and SD corresponding to the distance distribution
	     * within a protein of the same size. If this doesn't work, 
	     * use a simple uniform distribution.
	     */
	    Drand=portrandom_gauss();   // zero mean, unit variance
	    Drand=(Drand*Dev)+Avgdist;	// prescribed mean & variance
	    if (Drand<Low || Drand>Up)
		Drand=(Up-Low)*port_random()+Low;  // uniform

	    // adjust with hydrophobic estimate if "soft", skip N-, C-termini
	    if (i>0 && i<Ptno-1 && j>0 && !hard(i, j))
	    {
		Destim=Polymer.estim_dist(i-1, j-1);  // hydrophobic estimate
		if (Destim<0.0) continue;   // error?
		if (Destim>Up) Destim=0.95*Up;	// bracket if necessary
		if (Destim<Low) Destim=1.05*Low;
		Strict=strict(i, j)*Polymer.cons(i-1)*Polymer.cons(j-1); // strictness from conservation
		Drand=(1.0-Strict)*Drand+Strict*Destim;  // modified random dist
	    }
	    Dist[i][j]=Drand*Drand;
	}
    }
}
// END of init_distmat()

// ---- Auxiliaries ----

/* flory_constr: calculates the average end-to-end distances (UNsquared)
 * for a freely rotating chain (ref: Flory(1969)) and uses
 * these as upper residue dist limits. Stores the results in
 * the Maxsepar member array. 
 * Private
 */
void Restraints_::flory_constr()
{
    unsigned int Rno=Size-2;
    if (!Rno) return;
  
    /* Rexp is the expected radius: the quantity REXP_MAX*Rexp
     * will serve as an upper bound for the d>=2 entries.
     * REXP_MAX>=2.0 is recommended (max. dist is at least as
     * long as the expected diagonal).
     */
    double Rexp=exp_rad(Rno);
    static const double REXP_MAX=2.5, CA_BONDLEN_2=3.8*3.8;
    
    register double Rx, Alpha, Alpowd, C0, C1;
    register unsigned int d;

    // init 
    Rx=Rexp*REXP_MAX; Rx*=Rx;	// square max. distance
    Alpha=cos(M_PI-CA_BONDANGLE);
    C0=(1.0+Alpha)/(1.0-Alpha);	// constant term
    C1=2.0*Alpha/(1.0-Alpha)/(1.0-Alpha);   // coefficient
    
    Maxsepar[0]=0.0;		// self-distance is always 0
    Maxsepar[1]=CA_BONDLEN_2;	// next neighbour: bond length squared
    
    // cycle through all other diagonal values
    Alpowd=Alpha;   // Alpha^d calc'd cunningly in the cycle
    for (d=2; d<Size; d++)
    {
	Alpowd*=Alpha;	/* ^d */
	Maxsepar[d]=(C0-C1*(1.0-Alpowd)/d)*d*CA_BONDLEN_2;
	if (Maxsepar[d]>Rx) break;	// upper bound reached
    }
    for (; d<Size; d++) Maxsepar[d]=Rx;	// flat upper bound for rest
    
    // very primitively, take the square root of everybody... 21-May-96
    for (d=0; d<Size; d++) Maxsepar[d]=sqrt(Maxsepar[d]);
}
// END of flory_constr()

// ---- Input/output ----

/* Output: lists to Out the restraints. */
ostream& operator<<(ostream& Out, const Restraints_& D)
{
    if (!(D.Restrs)) Out<<"# <no external restraints>\n";
    else
    {
	Clist1_<Restr_> Rl=D.ext_restr(); // const access
	for (Rl.begin(); Rl!=NULL; Rl++) Out<<(*Rl);
    }
    return(Out);
}
// END of <<

/* read_restrs(): reads restraints from a file Fname. Assumes that
 * the file contains comments (lines beginning w/ '#') and
 * lines corresponding to the restraint format (see above in the
 * Restr_ class input for details), one line per restraint each.
 * If Fname==NULL or "" then the restraint list is cleared.
 * Checks the restraints for range conformance (residue nos in the
 * range [1..Rno]) (Rno comes from the length of the Maxsepar array internally.)
 * If Rno==0, then nothing is done (no restraints would pass).
 * Interatom restraints will be converted into CA and SCC
 * (sidechain centroid)-based restraints using the polymer 
 * description in Polymer and the atom:(CA|SCC) distances
 * in Acdist.
 * The old restraint list will be replaced with the new.
 * Return value: 0 on error, 1 if OK.
 */
int Restraints_::read_restrs(const char *Fname, const Polymer_& Polymer)
{
    if (!Polymer.len())
    {
	cerr<<"\n? Restraints_::read_restrs(): Please read sequence before restraints\n";
	return(0);
    }
    unsigned int Rno=Size-2;
    if (Polymer.len()!=Rno)
    {
	cerr<<"\n? Restraints_::read_restrs(): Polymer length mismatch ("<<Polymer.len()<<"!="<<Rno<<")\n";
	return(0);
    }
    if (Fname==NULL || !strlen(Fname))	// clear restraint list
    {
	Restrs.clear(); return(0);
    }
    
    ifstream Inf(Fname);    // open for reading
    if (!Inf)
    {
	cerr<<"\n? Restraints_::read_restrs(\""<<Fname<<"\"): Cannot open\n";
	return(0);
    }
    
    Inf>>(*this);   // errors are just skipped
    Inf.close(); 
    
    // convert interatomic distances into CA and SCC-based distances
    convert_restraints(Polymer);
    
    return(Inf.good() || Inf.eof());	// OK
}
// END of read_restrs()

/* Input: reads restraints (and comment lines beginning with '#')
 * from the input stream Inf. See >>Restr_ for the restraint line
 * format. The old restraint list is cleared. The new list after
 * input will contain a mixture of interatomic and CA:SCC restraints.
 * To convert to an all-CA:SCC list, invoke convert_restraints()
 * after >> (this is done automagically in read_restrs() above).
 * Note that the density is not set here!
 */
istream& operator>>(istream& Inf, Restraints_& Restraints)
{
    if (!Inf)
    {
	cerr<<"\n? >>Restraints_: Cannot read from stream\n";
	return(Inf);
    }
    
    unsigned int Rno=Restraints.Size-2;
    if (!Rno)
    {
	cerr<<"\n? >>Restraints_(): Please read sequence before restraints\n";
	return(Inf);
    }
    
    static const unsigned int LINELEN=120;
    char Line[LINELEN+1];
    istrstream Instr(Line, LINELEN+1);
    Restr_ Rs;
    unsigned int Lineno;
    
    Restraints.Restrs.clear(); // forget all previous restraints
    
    // scan all lines (zero Line before read to avoid UMRs in >>String_ )
    for (Lineno=0; Inf.good(); Lineno++)
    {
        memset(Line, 0, LINELEN+1);
	Inf.getline(Line, LINELEN, '\n');
	if (Line[0]=='#' || !strlen(Line))
	    continue; // skip comments and empty lines
	
	// read one line
	Instr.seekg(ios::beg); Instr.clear(); Instr>>Rs;
	if (!Instr) // error set by Restr_ input
	{
	    cerr<<"\n? >>Restraints_: Cannot read line "<<(Lineno+1)<<", skipped\n";
	    continue;
	}
	
	// check for "overhangs"
	if (!Rs.pos(1) || Rs.pos(1)>Rno || 
	    !Rs.pos(2) || Rs.pos(2)>Rno)
	{
	    cerr<<"\n? >>Restraints_: Residue number(s) "<<Rs.pos(1)<<", "<<Rs.pos(2)
		<<" in line "<<(Lineno+1)<<": out of range [1.."
		<<Rno<<"], skipped\n";
	    continue;
	}
	
	// OK, add to list 
	Restraints.Restrs+=Rs;
    }
    return(Inf);
}
// END of operator>> 

/* add_restrs(): adds the CA|SCC restraints specified in the list
 * Rs to the internal restraint list. Filters out those restraints
 * which are between side-chain atoms. Expects NON_SQUARED entries
 * in the list. Useful for restraints coming from the homology
 * modelling module. Call after read_restrs() since that method will
 * clear the restraint list.
 * Return value: the number of restraints processed.
 */
int Restraints_::add_restrs(const List1_<Restr_>& Rs)
{
    Clist1_<Restr_> Ra(Rs); // const iterator
    List1_<Restr_> Rnew;    // put "new" restraints here
    
    int Rsno=0;
    
    for (Ra.begin(); Ra!=NULL; Ra++)
    {
	// pass CA/SCC restraints straight on
	if ((Ra->atom(1)==String_("CA") || Ra->atom(1)==String_("SCC")) &&
	    (Ra->atom(2)==String_("CA") || Ra->atom(2)==String_("SCC")))
	{
	    if (!absorb_restraint(Restrs, *Ra)) 
		Rnew+=(*Ra);    // save on "new list" (could not be absorbed)
	}
	else
	{
	    cerr<<(*Ra)<<"\n? Restraints_::add_restrs(): This restraint is not between CA|SCC atoms\n";
	    continue;
	}
	Rsno++;
    }
    Restrs+=Rnew;   // append new list to original
    return(Rsno);
}
// END of add_restrs()

// ---- Restraint conversion ----

/* add_restraint(): adds restraint R to the internal restraint list
 * (which will contain CA/SCC restraints only). If there was a
 * restraint between the atoms of R, then it is "absorbed" into
 * Restrs (see absorb_restrain()). If R is the first restraint between its atoms, 
 * then it is simply appended to Restrs. Private
 */
inline
void Restraints_::add_restraint(const Restr_& R)
{
    if (!absorb_restraint(Restrs, R))
	Restrs+=R;	// first restraint between these atoms, append to list
    Restrs.end();
}
// END of add_restraint()

/* absorb_restraint(): check if the restraint R is between
 * the same CA/SCC atoms as one of the restraints in the list Rlist.
 * If this is the case, then the range of that restraint is 
 * modified: the restraint range is narrowed if R is stricter and either its
 * Lowlim is higher or Uplim is lower than that of the previous restraint
 * found in the list (R is "absorbed" into the list).
 * Otherwise no action is taken.
 * Return value: 1 if R could be "absorbed", 0 if not. Static private
 */
int Restraints_::absorb_restraint(List1_<Restr_>& Rlist, const Restr_& R)
{
    for (Rlist.begin(); Rlist!=NULL; Rlist++)
    {
	// is *Rlist between the same atoms as R?
	if (Rlist->pos(1)==R.pos(1) && Rlist->pos(2)==R.pos(2) &&	// X:Y == X:Y 
	    Rlist->atom(1)==R.atom(1) && Rlist->atom(2)==R.atom(2)
	    ||
	    Rlist->pos(1)==R.pos(2) && Rlist->pos(2)==R.pos(1) &&	// X:Y == Y:X 
	    Rlist->atom(1)==R.atom(2) && Rlist->atom(2)==R.atom(1))
	{
	    // try to narrow restraint range
	    if (Rlist->low()<=R.low() && Rlist->strict()<=R.strict())
	    {
		Rlist->low(R.low());
		Rlist->strict(R.strict());
	    }
	    if (Rlist->up()>=R.up() && Rlist->strict()<=R.strict())
	    {
		Rlist->up(R.up());
		Rlist->strict(R.strict());
	    }
	    
	    return(1);	// R is "absorbed" into *Rlist
	}
    }
    return(0);
}
// END of absorb_restraint()

/* convert_restraints(): converts the restraints in the internal restraint list
 * which may have been specified between side-chain atoms, into
 * (looser) restraints between C-alpha atoms (CA) and/or side
 * chain centroids (SCC). Currently DRAGON uses a reduced
 * chain representation with CAs and SCCs only.
 */
void Restraints_::convert_restraints(const Polymer_& Polymer)
{
    register float Cad1, Sccd1, Cad2, Sccd2;
    Restr_ R;
    
    /* What follows is an awful hack. The original restraint list
     * in the calling object is deep-copied to Ra and cleared;
     * Ra's items are then processed back into the internal list.
     */
    List1_<Restr_> Ra(Restrs);	// copy list into Ra
    Restrs.clear();	// clear previous restraints
    
    // scan all restraints
    for (Ra.begin(); Ra!=NULL; Ra++)
    {
	// obtain CA and SCC distances, skip if atom spec was wrong
	if (!get_cascc(Polymer, Ra->pos(1), Ra->atom(1), Cad1, Sccd1))
	    continue;
	if (!get_cascc(Polymer, Ra->pos(2), Ra->atom(2), Cad2, Sccd2))
	    continue;

	// pass CA/SCC restraints straight on
	if ((Ra->atom(1)==String_("CA") || Ra->atom(1)==String_("SCC")) &&
	    (Ra->atom(2)==String_("CA") || Ra->atom(2)==String_("SCC")))
	{
	    add_restraint(*Ra);
	    continue;
	}
	
	/* Generate 4 CA|SCC restraints from *Ra. The basic idea is that
	 * if the limits are given between atoms A1 and A2, the
	 * limits between X and Y (X, Y=[CA|SCC]) can be deduced
	 * by considering the following arrangements:
	 * 
	 * |<........minimal dist(A1, A2)........>|
	 * A1---X<....minimal dist(X, Y)...>Y----A2
	 * 
	 * and
	 * 
	 * |<........maximal dist(X, Y)..........>|
	 * X---A1<..maximal dist(A1, A2)...>A2----Y
	 * 
	 * taking the bumps into account if necessary.
	 * The resulting restraints are quite generous and it
	 * is possible to narrow them further by doing some
	 * extra geometry. THIS IS NOT IMPLEMENTED YET! 27-Nov-95
	 */
	// these members are the same for all
	R.pos(1, Ra->pos(1)); R.pos(2, Ra->pos(2));
	R.strict(Ra->strict());
	
	// make a CA:CA restraint
	R.atom(1, "CA"); R.atom(2, "CA");
	R.low(Ra->low()-(Cad1+Cad2));
	R.up(Ra->up()+Cad1+Cad2);
	add_restraint(R);
	
	// make a CA:SCC restraint
	R.atom(2, "SCC");
	R.low(Ra->low()-(Cad1+Sccd2));
	R.up(Ra->up()+Cad1+Sccd2);
	add_restraint(R);
	
	// make a SCC:CA restraint
	R.atom(1, "SCC"); R.atom(2, "CA");
	R.low(Ra->low()-(Sccd1+Cad2));
	R.up(Ra->up()+Sccd1+Cad2);
	add_restraint(R);
	
	// make a SCC:SCC restraint
	R.atom(2, "SCC");
	R.low(Ra->low()-(Sccd1+Sccd2));
	R.up(Ra->up()+Sccd1+Sccd2);
	add_restraint(R);
    }
}
// END of convert_restraints()

/* get_cascc(): returns the distance of Atom from CA and SCC
 * in Cad and Sccd, respectively, in the Pos-th position in 
 * the model sequence. The distances, AA identity etc. are
 * supplied by Polymer. Static private
 */
int Restraints_::get_cascc(const Polymer_& Polymer, unsigned int Pos, 
	const String_& Atom, float& Cad, float& Sccd)
{
    Cad=Polymer.ca_dist(Pos-1, Atom);
    if (Cad<0.0)
    {
	cerr<<"\n? Restraints_::get_cascc(): Nonexistant atom \""
	    <<Atom<<"\" specified for residue "
	    <<Polymer.aa(Pos)<<"-"<<Pos<<endl;
	return(0);
    }
    Sccd=Polymer.scc_dist(Pos-1, Atom);
    return(1);
}
// END of get_cascc()

// ==== END OF METHODS Restr.c++ ====
