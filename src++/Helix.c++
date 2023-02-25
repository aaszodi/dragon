// ==== PROJECT DRAGON: METHODS Helix.c++ ====

/* Alpha-helix topology and geometry class. Cf. "Segment.h"
 * and "Sstrbase.h" for base class information.
 */

// SGI C++, IRIX 6.2, 14. Aug. 1996. (C) Andras Aszodi

// ---- STANDARD HEADERS ----

#include <math.h>
#include <string.h>

// ---- UTILITY HEADERS ----

#include "Vector.h"
#include "Sqmat.h"
#include "Hirot.h"

// ---- MODULE HEADER ----

#include "Helix.h"

// ---- DEFINITIONS ----

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

// ==== Helix_ METHODS ====

// ---- Static data member initialisation ----

const int Helix_::HX310_DIAG=2;  // 3/10-helix H-bond (i,i+2)
const int Helix_::ALPHA_DIAG=3;  // alpha-helix H-bond (i,i+3)
const int Helix_::HXPI_DIAG=4;  // pi-helix H-bond (i,i+4)

// The helical params are from the Schulz/Schirmer book (1979)
const double Helix_::RADIUS_310=1.9,
    Helix_::PITCH_310=2.0,
    Helix_::TURN_310=2.09;
const double Helix_::RADIUS_ALPHA=2.3,
    Helix_::PITCH_ALPHA=1.5,
    Helix_::TURN_ALPHA=1.75; 
const double Helix_::RADIUS_PI=2.8,
    Helix_::PITCH_PI=1.1,
    Helix_::TURN_PI=1.46;

Array_<double> Helix_::Dist310(0);	// zero length
Array_<double> Helix_::Distalpha(0);
Array_<double> Helix_::Distpi(0);

// ---- Constructors ----

/* set_diag(): sets the diagonal phase of the helix according to the 
 * type stored in the Htype member. Protected 
 */
inline
void Helix_::set_diag()
{
    switch (Htype)
    {
	case HX310: Diag=HX310_DIAG; break;
	case HXPI: Diag=HXPI_DIAG; break;
	case ALPHA:
	default: Diag=ALPHA_DIAG; break;
    }
}
// END of set_diag()

/* Inits the helix to begin at Start and end at Stop to have type Ht
 * (the default is ALPHA).
 * Stop>=Start but the Linsegm_ base class ctor takes care of this.
 * Will be set to store 2 tetrahedral point sets.
 */
Helix_::Helix_(unsigned int Start, unsigned int Stop, Helixtype_ Ht):
	Linsegm_(Start, Stop), Sstrbase_(2), Htype(Ht)
{
    set_diag();
    if (!Beg)
    {
	// cerr<<"\n? Helix_("<<Start<<",..."): Reserved residue number (0)\n";
	Beg=1; End=Beg+Diag;
    }
    if ((End-Beg)<Diag)
    {
	// cerr<<"\n? Helix_("<<Start<<", "<<Stop<<"): Too short\n";
	End=Beg+Diag;
    }
    Id.len(End-Beg+1); Id.dim(3);
}

/* Inits the helix with a Linsegm_ object with type Ht (default ALPHA). */
Helix_::Helix_(const Linsegm_& Ls, Helixtype_ Ht): 
	Linsegm_(Ls), Sstrbase_(2), Htype(Ht)
{
    set_diag();
    if (!Beg)
    {
	// cerr<<"\n? Helix_(Linsegm_): Reserved residue number (0)\n";
	Beg=1; End=Beg+Diag;
    }
    if ((End-Beg)<Diag)
    {
	cerr<<"\n? Helix_(Linsegm_): Too short\n";
	End=Beg+Diag;
    }
    Id.len(End-Beg+1); Id.dim(3);
}

// ---- Virtual constructor ----

void Helix_::operator()(Sstrbase_*& Bptr) const { Bptr=new Helix_(*this); }

// ---- H-bond topology ----

/* hbond_prev(), hbond_next(): return the no. of the previous
 * or the next residue H-bonded to Res or -1 if there's no partner
 * (at helix ends) or -2 if Res is not a member of the helix.
 * A warning is also printed in this case.
 */
int Helix_::hbond_prev(unsigned int Res) const
{
    if (!member(Res))
    {
	cerr<<"? Helix_::hbond_prev(): Residue "<<Res<<" isn't a member\n";
	return(-2);
    }
    return((Res>=Beg+Diag)? int(Res)-Diag: -1);
}

int Helix_::hbond_next(unsigned int Res) const
{
    if (!member(Res))
    {
	cerr<<"? Helix_::hbond_next(): Residue "<<Res<<" isn't a member\n";
	return(-2);
    }
    return((Res<=End-Diag)? int(Res)+Diag: -1);
}
// END of hbond_prev(), hbond_next()

// ---- Tetrahedral points ----

/* make_ths(): for detangling, 2 tetrahedra will be fit on each helix
 * with point indices (B, B+2, E-3, E-1) and (B+1, B+3, E-2, E) where
 * B, E are the beginning and end indices of the helix, respectively.
 * These indices are stored in the Thedra array inherited from Sstrbase_.
 * 4-residue alpha helices have only 1 tetrahedron: (B, B+1, B+2, B+3), 
 * 5- and 6-residue helices have 2 tetrahedra but the index layout is
 * special.
 * Protected virtual (called by make_idstruct() only).
 */
void Helix_::make_ths()
{
    if (len()<4)    // too short
    {
	Thedra.len(0); return;
    }
    
    switch(len())
    {
	case 4:	    // 1 tetrahedron only
	Thedra.len(1);
	Thedra[0].P1=beg(); Thedra[0].P2=beg()+1; 
	Thedra[0].P3=beg()+2; Thedra[0].P4=beg()+3;
	break;
	
	case 5:	    // 2 tetrahedra, but overlapping on middle 3 points
	Thedra.len(2);
	Thedra[0].P1=beg();
	Thedra[0].P2=Thedra[1].P1=beg()+1;
	Thedra[0].P3=Thedra[1].P2=beg()+2;
	Thedra[0].P4=Thedra[1].P3=beg()+3;
	Thedra[1].P4=end();
	break;
	
	case 6:	    // 2 tetrahedra, overlapping on middle 2 points
	Thedra.len(2);
	Thedra[0].P1=beg(); Thedra[1].P1=beg()+1;
	Thedra[0].P2=Thedra[1].P2=beg()+2;
	Thedra[0].P3=Thedra[1].P3=beg()+3;
	Thedra[0].P4=end(); Thedra[1].P4=end()-1;
	break;
	
	default:    // general case with length>=7
	Thedra.len(2);
	Thedra[0].P1=beg(); Thedra[1].P1=beg()+1;
	Thedra[0].P2=beg()+2; Thedra[1].P2=beg()+3;
	Thedra[0].P3=end()-3; Thedra[1].P3=end()-2;
	Thedra[0].P4=end()-1; Thedra[1].P4=end();
	break;
    }
}
// END of make_ths()

// ---- Ideal geometry ----

/* Updates the ideal distance array which is one of three static
 * members, one for each type of helix. If the calling object contains
 * a helix longer than the corresponding dist array, then the missing
 * values are calculated and put onto the end. Protected
 */
inline
void Helix_::update_iddist(Array_<double>& Dist) const
{
    register unsigned int Oldlen, i;
    if ((Oldlen=Dist.len())<len())
    {
	Dist.len(len());    // grow the array
	for (i=Oldlen; i<len(); i++)
	    Dist[i]=diff_len(Id[0], Id[i]);	// add new UNsquared values
    }
}
// END of update_iddist()

/* make_idstruct(): generates the 3D ideal right-handed helix
 * in Id if the sentinel Changed (inherited from Segmbase_) is true.
 * Returns length or 0 if something fails. 
 */
unsigned int Helix_::make_idstruct()
{
    if (!Changed) return(Id.len());   // no action needed
    
    // construct the tetrahedron index array
    make_ths();
    
    Id.len_dim(len(), 3);
    
    // make a helix with appropriate radius, pitch and turn
    unsigned int Retval;
    switch(Htype)
    {
	case HX310:
	    Retval=make_helix(Id, RADIUS_310, PITCH_310, TURN_310);
	break;
	case HXPI:
	    Retval=make_helix(Id, RADIUS_PI, PITCH_PI, TURN_PI);
	break;
	case ALPHA:
	default:
	    Retval=make_helix(Id, RADIUS_ALPHA, PITCH_ALPHA, TURN_ALPHA);
	break;
    }
    

    // center locally
    if (Retval)
    {
	Vector_ Ctr=Id.centroid();
	Id-=Ctr; Changed=false;    // also reset sentinel
    }
    
    // check if the ideal distances (one array for each type) need updating
    switch(Htype)
    {
	case HX310:
	    update_iddist(Dist310);
	break;
	case HXPI:
	    update_iddist(Distpi);
	break;
	case ALPHA:
	default:
	    update_iddist(Distalpha);
	break;
    }
    
    return(Retval);
}
// END of make_idstruct()

/* copy_iddist(): copies the ideal distances (held in Dist,
 * which is one of the static dist arrays for each type of helix)
 * into a distance matrix Dmat and updates Strimat with strictness Strict
 * (which is a member of Sstrbase_). Protected
 */
inline
void Helix_::copy_iddist(Trimat_& Dmat, Trimat_& Strimat, 
	const Array_<double>& Dist) const
{
    register unsigned int d, i, j;
    
    for (d=0; d<len(); d++)
    {
	for (i=beg()+d; i<=end(); i++)
	{
	    j=i-d;
	    if (Strimat[i][j]<=Strict)	// previous was weaker
	    {
		Dmat[i][j]=Dist[d];    // apply Dist[d]
		Strimat[i][j]=Strict;   // with strictness Strict
	    }
	}
    }
}
// END of copy_iddist()

/* ideal_dist(): puts the ideal alpha-helical UNsquared distances into
 * a distance matrix Dmat in the right position. Does nothing if
 * the helix does not fit in the matrix. Prints a warning if Changed==true, 
 * since this indicates that the size was changed without updating the
 * ideal structure and therefore what is returned may be incorrect.
 * (The situation will be dealt with elegantly when the "mutable"
 * keyword finds its way into the compiler.) The adjustment is
 * applied with a variable strictness (Strict is in Sstrbase_ ).
 * Also sets the corresponding strictness entries in Strimat.
 */
void Helix_::ideal_dist(Trimat_& Dmat, Trimat_& Strimat) const
{
    // check update status
    if (Changed)
    {
	cerr<<"\n? Helix_::ideal_dist(): make_idstruct() should have been called\n";
	return;
    }
    
    // check sizes
    if (Dmat.rno()<=end() || Strimat.rno()<=end())
    {
	cerr<<"\n? Helix_::ideal_dist(): Matrix too small\n";
	return;
    }
    
    // copy ideal squared distances with prescribed strictness
    switch(Htype)
    {
	case HX310: copy_iddist(Dmat, Strimat, Dist310); break;
	case HXPI: copy_iddist(Dmat, Strimat, Distpi); break;
	case ALPHA:
	default: copy_iddist(Dmat, Strimat, Distalpha); break;
    }
}
// END of ideal_dist2()

/* ideal_struct(): applies the ideal helical coordinates stored
 * inside onto the point set Model. Model must be large enough to contain
 * the helix and when masked with the helix's mask, the active region
 * must be 3-dimensional. If this is the case, then the ideal structure
 * will be RMS fitted onto Model's active region, the original segment
 * replaced by the rotated/transposed ideal at the given strictness
 * in Sstrbase, and the RMS value returned.
 * -1.0 is returned on error. Model's original activation pattern is
 * always retained. Prints a warning if an update is needed.
 */
double Helix_::ideal_struct(Points_& Model) const
{
    // check update status
    if (Changed)
    {
	cerr<<"\n? Helix_::ideal_struct(): make_idstruct() should have been called\n";
	return(-1.0);
    }
    
    // mask the point array with own mask, check local dimension
    if (Model.len()<=end())
    {
	cerr<<"\n? Helix_::ideal_struct(): Does not fit in\n";
	return(-1.0);
    }
    Bits_ Oldmask=Model.mask(mask(Model.len()));
    if (Model.dim()!=3)
    {
	// cerr<<"\n? Helix_::ideal_struct(): Model is not 3D\n";
	Model.mask(Oldmask); return(-1.0);
    }
    
    // center model
    Vector_ Mctr=Model.centroid();
    Model-=Mctr;
    
    // perform a "best RMS rotation"
    Hirot_ Hr;
    Hr.best_rot(Id, Model);
    double Rms=Hr.get_rms(Id, Model);
    if (Rms<0.0) return(Rms);	// something disastrous happened in hi_rot()
    
    /* now replace the active segment in Model by the rotated Id
     * using the Strict weight
     */
    float Strict1=1.0-Strict;
    for (register unsigned int i=0; i<len(); i++)
	Model[i]=Strict1*Model[i]+Strict*Hr.rot_matrix()*Id[i];
    Model+=Mctr;    // transpose back to original centroid
    
    Model.mask(Oldmask);	// reset original mask
    return(Rms);
}
// END of ideal_struct()

// ---- Handedness ----

/* check_torsion(): walks over the helix in 3D in Model and calculates all
 * (i, i+3) torsion angles. For right-handed helices these should
 * be all positive. Good and Bad will be set to the no. of correct
 * and incorrect torsion angles.
 * Return value: 1 if Good>=Bad, -1 if Good<Bad, 0 if not in 3D.
 */
int Helix_::check_torsion(Points_& Model, 
	unsigned int& Good, unsigned int& Bad) const
{
    Bits_ Oldmask=Model.mask(mask(Model.len()));
    if (Model.dim()!=3)
    {
	Model.mask(Oldmask);
	return(0);
    }
    
    register unsigned int i;
    register double Tors;
    
    Good=Bad=0;
    for (i=0; i+3<len(); i++)
    {
	Tors=pos4_angle(Model[i], Model[i+1], Model[i+2], Model[i+3]);
	if (-M_PI>Tors)
	{
	    cerr<<"\n? Helix_::check_torsion(): collinearity\n";
	    continue;
	}
	if (Tors<0.0) ++Bad; else ++Good;
    }
    
    Model.mask(Oldmask);
    return((Good>=Bad)? 1: -1);
}
// END of check_torsion()

// ---- Input ----

/* >>: reads in a helix from the stream In. Expects the following
 * format:-
 * "<str> <beg> <end> [strict]\n" 
 * where the spaces represent one or more
 * whitespaces, <beg> and <end> are positive integers. 
 * <str> is one of the following strings: "HX310" for a 3/10 helix, 
 * "ALPHA" or "HELIX" for an alpha-helix, "HXPI" for a pi-helix.
 * [strict] is an optional floating-point number between 0.0 and 1.0, 
 * which, if specified, tells the system to what extent to apply the ideal
 * helical structure in the model.
 * Input stops after <end> or [strict] and the newline is not consumed.
 * On error, the "fail" bit is set and the Helix_ object will not
 * be modified. 
 */
istream& operator>>(istream& In, Helix_& H)
{
    if (!In) return(In);    // don't bother with unreadable streams
    
    char Hbuf[6];
    Helix_::Helixtype_ Ht;
    In>>setw(6)>>Hbuf;
    if (!In)
    {
	cerr<<"\n? >>Helix_: string (\"HX310\", \"ALPHA\", \"HELIX\" or \"HXPI\") expected\n";
	In.clear(ios::failbit|In.rdstate());
	return(In);
    }
    
    if (!strcmp(Hbuf, "HX310")) Ht=Helix_::HX310;
    else if (!strcmp(Hbuf, "ALPHA") || !strcmp(Hbuf, "HELIX")) Ht=Helix_::ALPHA;
    else if (!strcmp(Hbuf, "HXPI")) Ht=Helix_::HXPI;
    else
    {
	cerr<<"\n? >>Helix_: Invalid helix type \""<<Hbuf<<"\"\n";
	In.clear(ios::failbit|In.rdstate());
	return(In);
    }
    
    int B=0, E=0;
    In>>B>>E;
    if (!In || B<=0 || E<=0)
    {
	cerr<<"\n? >>Helix_: Invalid limits: "<<B<<", "<<E<<endl;
	In.clear(ios::failbit|In.rdstate());
	return(In);
    }
    
    // check if strictness was specified; if yes, read it
    float Str=1.0;  // default
    int Rdstate=In.rdstate();	// save previous state
    In>>Str;
    if (!In)	// was not there: restore previous state
	In.clear(Rdstate);
    
    if (Str<=0.0)
    {
	cerr<<"\n? >>Helix_: Strictness "<<Str<<"<=0.0, helix will be ignored\n";
	In.clear(ios::failbit|In.rdstate());
	return(In);
    }
    if (Str>1.0) Str=1.0;   // adjust silently
    
    // OK, update the helix
    if (B>E)
    {
	int T=B; B=E; E=T;  // swap
    }
    H.limits(B, E); // keep [1..Rno] range V4.8.1
    H.Strict=Str;   // store strictness
    H.Htype=Ht;	// store type first
    H.set_diag();	// and then adjust the diagonal phase
    return(In);
}
// END of >>

// ---- Output ----

/* write_to(): writes the calling object to the stream Out.
 * The format is: "<type> <begin> <end> [strict]\n" where the residue
 * numbers start with 1 and <type> is "HX310", "ALPHA" or "HXPI".
 * [strict] is printed only if it is not 1.0.
 */
void Helix_::write_to(ostream& Out) const
{
    switch(Htype)
    {
	case HX310: Out<<"HX310 "; break;
	case HXPI: Out<<"HXPI "; break;
	case ALPHA: 
	default: Out<<"ALPHA "; break;
    }
    Out<<beg()<<' '<<end();
    if (Strict!=1.0) Out<<' '<<Strict;
    Out<<endl;
}
// END of write_to()

// ==== END OF METHODS Helix.c++ ====
