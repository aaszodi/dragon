// ==== PROJECT DRAGON: METHODS Secstr.c++ ====

/* Classes for handling secondary structures: H-bond topology
 * and ideal geometry. Chain topology comes from the Segment module.
 */

// SGI C++ 4.0, IRIX 5.3, 15. Aug. 1995. (C) Andras Aszodi

// ---- STANDARD HEADERS ----

#include <math.h>
#include <string.h>

// ---- UTILITY HEADERS ----

#include "Vector.h"
#include "Sqmat.h"
#include "Hirot.h"

// ---- MODULE HEADER ----

#include "Secstr.h"

// ---- DEFINITIONS ----

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

// ==== Sstrbase_ METHODS ====

/* make_helix(): constructs an ideal helix in a point set Hel which
 * must be properly masked and the active region should be 3-dimensional.
 * The helix parameters are supplied by Radius, Pitch and Turn and
 * the structure will be grown so that the N->C direction corresponds
 * to the positive direction on the X axis. Phasing (default 1)
 * determines whether the first point is on +Y (>0) or -Y (when <=0).
 * Return value: the length of the helix built or 0 on error. Protected static
 */
unsigned int Sstrbase_::make_helix(Points_& Hel, double Radius, 
	double Pitch, double Turn, int Phasing)
{
    if (Hel.dim()!=3)
    {
	cerr<<"\n? Sstrbase_::make_helix(): Sorry,  not 3D\n";
	return(0);
    }
    
    unsigned int L=Hel.active_len();
    if (!L) return(0);	// totally inactive
    
    // flip structure around X-axis (makes sense for beta-strands only)
    if (Phasing<=0) Radius*=(-1);
    
    for (register unsigned int i=0; i<L; i++)
    {
	Hel[i][0]=i*Pitch;  // X, Y, Z coords
	Hel[i][1]=Radius*cos(i*Turn);
	Hel[i][2]=Radius*sin(i*Turn);
    }
    return(L);
}
// END of make_helix()

/* pos4_angle: given the Cartesian coordinates of 4 points (P[1-4]),
 * the torsion angle defined by them (along 2-3) is returned. The
 * value should be between -Pi..+Pi with the sign indicating the
 * usual handedness convention. -2*Pi is returned if any 3 of the 
 * points are colinear. Protected static
 */
double Sstrbase_::pos4_angle(const Vector_& P1, const Vector_& P2,
	const Vector_& P3, const Vector_& P4)
{
    // get normal vectors of planes (123) and (234)
    Vector_ V2=P3-P2;   // torsion is taken along this
    Vector_ W1=cross_prod(P2-P1,V2);
    Vector_ W2=cross_prod(V2,P4-P3);

    // colinearity check: return -2Pi if fails
    register double W1len, W2len;

    if (0.0==(W1len=W1.vec_len()) ||
	0.0==(W2len=W2.vec_len()))
	    return(-2.0*M_PI);

    // get torsion angle (i.e. angle between planes (123),(234))
    register double Costheta;
    
    Costheta=(W1*W2)/(W1len*W2len);
    Costheta=acos(Costheta);

    // get sign of theta (handedness)
    return(((V2*cross_prod(W1, W2))>=0.0)? Costheta: -Costheta);
}
// END of pos4_angle()

/* <<: calls a virtual output function. Details in derived classes */
ostream& operator<<(ostream& Out, const Sstrbase_& S)
{
    S.write_to(Out);
    return(Out);
}
// END of <<

// ==== Helix_ METHODS ====

// ---- Static data member initialisation ----

const int Helix_::HELIX_ALPHA_DIAG=3;  // alpha-helix H-bond (i,i+3)
Array_<double> Helix_::Dist(0);	// use default init (zero length)

// ---- Constructors ----

/* Inits the helix to begin at Start and end at Stop.
 * Stop>=Start but the Linsegm_ base class ctor takes care of this.
 * Will be set to store 2 tetrahedral point sets.
 * There must be at least HELIX_ALPHA_DIAG residues.
 */
Helix_::Helix_(unsigned int Start, unsigned int Stop):
	Linsegm_(Start, Stop), Sstrbase_(2)
{
    if ((End-Beg)<HELIX_ALPHA_DIAG)
    {
	// cerr<<"\n? Helix_("<<Start<<", "<<Stop<<"): Too short\n";
	End=Beg+HELIX_ALPHA_DIAG;
    }
    Id.len(End-Beg+1); Id.dim(3);
}

/* Inits the helix with a Linsegm_ object. */
Helix_::Helix_(const Linsegm_& Ls): 
	Linsegm_(Ls), Sstrbase_(2)
{
    if ((End-Beg)<HELIX_ALPHA_DIAG)
    {
	cerr<<"\n? Helix_(Linsegm_): Too short\n";
	End=Beg+HELIX_ALPHA_DIAG;
    }
    Id.len(End-Beg+1); Id.dim(3);
}

// ---- Virtual constructors ----

void Helix_::operator()(Segmbase_*& Bptr) const { Bptr=new Helix_(*this); }
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
    return((Res>=Beg+HELIX_ALPHA_DIAG)? int(Res)-HELIX_ALPHA_DIAG: -1);
}

int Helix_::hbond_next(unsigned int Res) const
{
    if (!member(Res))
    {
	cerr<<"? Helix_::hbond_next(): Residue "<<Res<<" isn't a member\n";
	return(-2);
    }
    return((Res<=End-HELIX_ALPHA_DIAG)? int(Res)+HELIX_ALPHA_DIAG: -1);
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

/* make_idstruct(): generates the 3D ideal right-handed alpha-helix
 * in Id if the sentinel Changed (inherited from Segmbase_) is true.
 * Returns length or 0 if something fails. 
 */
unsigned int Helix_::make_idstruct()
{
    if (!Changed) return(Id.len());   // no action needed
    
    // construct the tetrahedron index array
    make_ths();
    
    // alpha-helix parameters
    static const double RADIUS=2.29, PITCH=1.50, TURN=1.75;
    
    Id.len(len()); Id.mask(true); Id.dim(3);
    
    // make a helix with alpha-helical radius, pitch and turn
    unsigned int Retval=make_helix(Id, RADIUS, PITCH, TURN);

    // center locally
    if (Retval)
    {
	Vector_ Ctr=Id.centroid();
	Id-=Ctr; Changed=false;    // also reset sentinel
    }
    
    // check if the ideal distances need updating
    unsigned int Oldlen, i;
        
    if ((Oldlen=Dist.len())<len())
    {
	Dist.len(len());    // grow the array
	for (i=Oldlen; i<len(); i++)
	    Dist[i]=diff_len2(Id[0], Id[i]);	// add new values
    }
    
    return(Retval);
}
// END of make_idstruct()

/* ideal_dist2(): puts the ideal alpha-helical squared distances into
 * a distance matrix Dmat in the right position. Does nothing if
 * the helix does not fit in the matrix. Prints a warning if Changed==true, 
 * since this indicates that the size was changed without updating the
 * ideal structure and therefore what is returned may be incorrect.
 * (The situation will be dealt with elegantly when the "mutable"
 * keyword finds its way into the compiler.)
 */
void Helix_::ideal_dist2(Trimat_& Dmat) const
{
    // check update status
    if (Changed)
    {
	cerr<<"\n? Helix_::ideal_dist2(): make_idstruct() should have been called\n";
	return;
    }
    
    // check sizes
    if (Dmat.rno()<=end())
    {
	cerr<<"\n? Helix_::ideal_dist2(): Matrix too small\n";
	return;
    }
    
    // copy ideal squared distances
    register unsigned int d, i, j;
    
    for (d=0; d<len(); d++)
	for (i=beg()+d; i<=end(); i++)
	{
	    j=i-d; Dmat[i][j]=Dist[d];
	}
}
// END of ideal_dist2()

/* ideal_struct(): applies the ideal alpha-helical coordinates stored
 * inside onto the point set Model. Model must be large enough to contain
 * the helix and when masked with the helix's mask, the active region
 * must be 3-dimensional. If this is the case, then the ideal structure
 * will be RMS fitted onto Model's active region, the original segment
 * replaced by the rotated/transposed ideal, and the RMS value returned.
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
    Sqmat_ Rot(3);
    Vector_ W(len()); W.set_values(1.0);    // uniform weight
    hi_rot(Id, Model, W, Rot);
    double Rms=get_rms(Id, Model, W, Rot);
    if (Rms<0.0) return(Rms);	// something disastrous happened in hi_rot()
    
    // now replace the active segment in Model by the rotated Id
    for (register unsigned int i=0; i<len(); i++)
	Model[i]=Id[i];
    Model*=Rot;	    // rotate ideal
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
 * " HELIX <beg> <end> \n" 
 * where the spaces represent one or more
 * whitespaces, <beg> and <end> are positive integers. Input stops
 * after <end> and the newline is not consumed. On error, 
 * the "fail" bit is set and the Helix_ object will not
 * be modified. 
 */
istream& operator>>(istream& In, Helix_& H)
{
    if (!In) return(In);    // don't bother with unreadable streams
    
    char Hbuf[6];
    In>>setw(6)>>Hbuf;
    if (!In || strcmp(Hbuf, "HELIX"))
    {
	cerr<<"\n? >>Helix_: Invalid descriptor "<<Hbuf<<": HELIX expected\n";
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
    
    // OK, update the helix
    if (B>E)
    {
	int T=B; B=E; E=T;  // swap
    }
    H.limits(B-1, E-1);
    return(In);
}
// END of >>

// ---- Output ----

/* write_to(): writes the calling object to the stream Out.
 * The format is: "HELIX <begin> <end>\n" where the residue
 * numbers start with 1. Protected
 */
void Helix_::write_to(ostream& Out) const
{
    Out<<"HELIX "<<(beg()+1)<<' '<<(end()+1)<<endl;
}
// END of write_to()

// ==== Beta_ METHODS ====

// ---- Constructors ----

/* Inits with a sheet. */
Beta_::Beta_(const Sheet_& Sh): Sheet_(Sh)
{
    unsigned int L=mask().len();
    Idup.len(L); Idup.mask(true); Idup.dim(3); 
    Iddown.len(L); Iddown.mask(true); Iddown.dim(3);
    Dist.set_size(L);
}

// ---- Virtual constructors ----

void Beta_::operator()(Segmbase_*& Bptr) const { Bptr=new Beta_(*this); }
void Beta_::operator()(Sstrbase_*& Bptr) const { Bptr=new Beta_(*this); }

// ---- H-bond topology ----

/* hbond_prev(), hbond_next(): return the H-bonding partners of residue Resno.
 * Return value is -1 if Resno is at the edge of the sheet and hasn't got
 * a partner in that particular direction: -2 if Resno is not a member of
 * the sheet at all (a warning is also printed in this case).
 */
int Beta_::hbond_prev(unsigned int Resno) const
{
    int Idx=strand_res(Resno);
    if (Idx<0)
    {
	cerr<<"? Sheet_::hbond_prev(): Residue "<<Resno<<" is not in sheet\n";
	return(-2);
    }
    if (!Idx) return(-1);   // Resno in first strand, can't have prev. partner

    int Prev;
    if (Strands[Idx].sense()==Strand_::PAR)	// parallel
	Prev=Strands[Idx-1].beg()+Strands[Idx].phase()+(Resno-Strands[Idx].beg());
    else	// anti
	Prev=Strands[Idx-1].end()-Strands[Idx].phase()-(Resno-Strands[Idx].beg());
    return((Strands[Idx-1].member(Prev))? Prev: -1);	// may be on an overhang
}
// END of hbond_prev()

int Beta_::hbond_next(unsigned int Resno) const
{
    int Idx=strand_res(Resno);
    if (Idx<0)
    {
	cerr<<"? Sheet_::hbond_next(): Residue "<<Resno<<" is not in sheet\n";
	return(-2);
    }
    if (Idx==Strands.len()-1) return(-1);   // Resno in last strand, can't have next. partner
    
    int Next;
    if (Strands[Idx+1].sense()==Strand_::PAR)	// parallel
	Next=Resno-Strands[Idx].beg()-Strands[Idx+1].phase()+Strands[Idx+1].beg();
    else    // anti
	Next=Strands[Idx].end()-Resno-Strands[Idx+1].phase()+Strands[Idx+1].beg();
    return((Strands[Idx+1].member(Next))? Next: -1);	// may be on an overhang
}
// END of hbond_next()

// ---- Tetrahedral points ----

/* make_ths(): makes the array of tetrahedra indices used by the
 * detangling routines. A tetrahedron is spanned by the end-points
 * of two neigbouring strands so there are S-1 tetrahedra for a
 * S-strand sheet.
 * Protected virtual (called by make_idstruct() only).
 */
void Beta_::make_ths()
{
    if (strand_no()<=1)
    {
	Thedra.len(0); return;
    }
    
    Thedra.len(strand_no()-1);
    for (register unsigned int i=0; i<strand_no()-1; i++)
    {
	Thedra[i].P1=Strands[i].beg(); Thedra[i].P2=Strands[i].end();
	Thedra[i].P3=Strands[i+1].beg(); Thedra[i].P4=Strands[i+1].end();
    }
}

// ---- Ideal geometry ----

/* make_idstruct(): generates two ideal beta-sheets ("up" and "down"
 * phasing) if Changed (inherited from Segmbase_) is true, and puts
 * the 3D coordinates into Idup and Iddown. Returns no. of residues or
 * 0 if something went wrong.
 */
unsigned int Beta_::make_idstruct()
{
    // beta geometry (helical params, strand separation, twist angle)
    static const double RADIUS=0.96, PITCH=3.32, TURN=3.25, 
	STRSEP=4.90, TW_ANGLE=-0.349;
    
    if (!Changed) return(mask().on_no());	// do nothing
    
    // construct tetrahedral points
    make_ths();
    
    // determine the maximal width of the sheet
    register unsigned int i, Sno=strand_no();
    int Minoffs=INT_MAX, Maxoffs=-INT_MAX, Eoffs;
    int *Boffs=new int [Sno];
    
    for (i=0; i<Sno; i++)
    {
	Boffs[i]=offs_strd(i, 0);
	if (Boffs[i]<Minoffs) Minoffs=Boffs[i];
	if (Boffs[i]>Maxoffs) Maxoffs=Boffs[i];
	
	Eoffs=offs_strd(i, Strands[i].end()-Strands[i].beg());
	if (Eoffs<Minoffs) Minoffs=Eoffs;
	if (Eoffs>Maxoffs) Maxoffs=Eoffs;
	
	// swap end offsets for strands which are anti wrt the first
	if (sense(0, i)==Strand_::ANTI) Boffs[i]=Eoffs;
    }
    
    // generate an "up" and "down" long strand as prototypes of all strands
    unsigned int Width=Maxoffs-Minoffs+1;
    Points_ Protoup(Width), Protodown(Width);   // automatically 3D
    Sqmat_ Rot(3);	// for rotating around X
    double Xangcorr=(TURN-M_PI)*Width/2.0;
    
    make_helix(Protoup, RADIUS, PITCH, TURN, 1);   // ideal beta-strand
    make_helix(Protodown, RADIUS, PITCH, TURN, -1);
    
    /* rotate around the X-axis so that the middle portion will
     * be approximately orthogonal to the X:Z plane
     */
    Rot[0][0]=1.0;
    Rot[1][1]=Rot[2][2]=cos(Xangcorr);
    Rot[1][2]=Rot[2][1]=sin(Xangcorr);
    Rot[2][1]*=(-1.0);
    Protoup*=Rot; Protodown*=Rot;
    
    /* Adjust the sizes of the ideal coordinate arrays.
     * Note that the arrays must be long enough to accommodate the
     * full sheet and in general will contain lots of unused vectors.
     * This is wasteful but compact storage would be tedious to implement.
     */
    Bits_ Betamask=mask();
    unsigned int L=Betamask.len();
    
    Idup.len(L); Idup.mask(Betamask); Idup.dim(3);
    Iddown.len(L); Iddown.mask(Betamask); Iddown.dim(3);
    
    /* Now copy appropriate portions of the prototype strand into
     * the strands. The first strand gets an exact copy, the
     * next will be shifted along the Z-axis by STRSEP etc.
     */
    register unsigned int j, Actlen;
    Bits_ Strmask;	// strand mask
    register double Strshift;	// strand shift
    int Dir;	// strand direction wrt first
    
    for (i=0, Dir=1; i<Sno; i++)
    {
	Strmask=Strands[i].mask(L);
	Idup.mask(Strmask);	// work on current strand only
	Iddown.mask(Strmask);
	Actlen=Idup.active_len();   // no. of residues in current strand
	
	if (Strands[i].sense()==Strand_::ANTI) Dir*=(-1);   // swap direction
	
	if (Dir>=0)
	{	// parallel to first
	    for (j=0; j<Actlen; j++)
	    {
		Idup[j]=Protoup[j+Boffs[i]-Minoffs];	// copy coords
		Iddown[j]=Protodown[j+Boffs[i]-Minoffs];
	    }
	}
	else
	{	// antiparallel to first: reverse copy
	    for (j=0; j<Actlen; j++)
	    {
		Idup[Actlen-j-1]=Protoup[j+Boffs[i]-Minoffs];	// copy coords
		Iddown[Actlen-j-1]=Protodown[j+Boffs[i]-Minoffs];
	    }
	}
	
	// add strand separation shift
	Strshift=i*STRSEP;
	for (j=0; j<Idup.active_len(); j++)
	{
	    Idup[j][2]+=Strshift;
	    Iddown[j][2]+=Strshift;
	}
    }

    // mask to sheet and center
    Idup.mask(Betamask);
    Vector_ Ctr=Idup.centroid();
    Idup-=Ctr;
    Iddown.mask(Betamask);
    Ctr=Iddown.centroid();
    Iddown-=Ctr;
    
    // add the sheet twist by rotating the strands around the Z-axis
    Rot.set_values(); Rot[2][2]=1.0;
    
    for (i=1; i<Sno; i++)
    {
	// get rotation matrix: could be faster
	Rot[0][0]=Rot[1][1]=cos(TW_ANGLE*i);
	Rot[1][0]=Rot[0][1]=sin(TW_ANGLE*i);
	Rot[0][1]*=(-1.0);
	
	// rotate the strand
	Strmask=Strands[i].mask(L);
	Idup.mask(Strmask); Iddown.mask(Strmask); 
	Idup*=Rot; Iddown*=Rot;
    }
    
    // generate ideal squared distances from the "Idup" coordinates
    Dist.set_size(L); Dist.set_values();
    Idup.mask(true);	// fully accessible, use real indices
    for (i=0; i<L; i++)
    {
	if (!member(i)) continue;
	for (j=0; j<=i; j++)
	{
	    if (!member(j)) continue;
	    Dist[i][j]=diff_len2(Idup[i], Idup[j]);
	}
    }
    
    // mask to sheet and center again
    Idup.mask(Betamask);
    Ctr=Idup.centroid(); Idup-=Ctr;
    Iddown.mask(Betamask);
    Ctr=Iddown.centroid(); Iddown-=Ctr;
    
    // clean up
    Changed=false;  // reset sentinel
    delete [] Boffs;
    return(Betamask.on_no());
}
// END of make_idstruct()

/* ideal_dist2(): puts the ideal beta-sheet squared distances into
 * a distance matrix Dmat in the right position. Does nothing if
 * the sheet does not fit in the matrix. Prints a warning if Changed==true, 
 * since this indicates that the size was changed without updating the
 * ideal structure and therefore what is returned may be incorrect.
 * (The situation will be dealt with elegantly when the "mutable"
 * keyword finds its way into the compiler.)
 */
void Beta_::ideal_dist2(Trimat_& Dmat) const
{
    // check update status
    if (Changed)
    {
	cerr<<"\n? Beta_::ideal_dist2(): make_idstruct() should have been called\n";
	return;
    }
    
    // size check
    if (Dmat.rno()<Dist.rno())
    {
	cerr<<"\n? Beta_::ideal_dist2(): Matrix too small\n";
	return;
    }
    
    // copy the ideal coordinates
    register unsigned int i, j;
    
    for (i=0; i<Dist.rno(); i++)
    {
	if (!member(i)) continue;
	
	for (j=0; j<=i; j++)
	{
	    if (!member(j)) continue;
	    Dmat[i][j]=Dist[i][j];
	}
    }
}
// END of ideal_dist2()

/* ideal_struct(): applies the ideal sheet coordinates stored
 * inside onto the point set Model. Model must be large enough to contain
 * the sheet and when masked with the sheet's mask, the active region
 * must be 3-dimensional. If this is the case, then the ideal structure
 * will be RMS fitted onto Model's active region, the original segment
 * replaced by the rotated/transposed ideal, and the RMS value returned.
 * The phasing will be chosen so that the ideal sheet with a better RMS
 * will be fitted (no knowledge of the phasing is necessary).
 * -1.0 is returned on error. Model's original activation pattern is
 * always retained. Prints a warning if an update is needed.
 */
double Beta_::ideal_struct(Points_& Model) const
{
    // check update status
    if (Changed)
    {
	cerr<<"\n? Beta_::ideal_struct(): make_idstruct() should have been called\n";
	return(-1.0);
    }
    
    // mask the point array with own mask, check local dimension
    Bits_ Betamask=mask();
    if (Model.len()<Betamask.len())
    {
	cerr<<"\n? Beta_::ideal_struct(): Does not fit in\n";
	return(-1.0);
    }
    Betamask.len(Model.len());	// adjust mask length to model data length
    Bits_ Oldmask=Model.mask(Betamask);
    if (Model.dim()!=3)
    {
	// cerr<<"\n? Beta_::ideal_struct(): Model is not 3D\n";
	Model.mask(Oldmask); return(-1.0);
    }
    
    // center model
    Vector_ Mctr=Model.centroid();
    Model-=Mctr;
    
    // perform "best RMS rotation" for both phasings
    Sqmat_ Rotup(3), Rotdown(3);
    Vector_ W(Betamask.len()); W.set_values(1.0);    // uniform weight
    
    hi_rot(Idup, Model, W, Rotup);
    hi_rot(Iddown, Model, W, Rotdown);
    
    double Rmsup, Rmsdown;
    
    Rmsup=get_rms(Idup, Model, W, Rotup);
    Rmsdown=get_rms(Iddown, Model, W, Rotdown);
    if (Rmsup<0.0 || Rmsdown<0.0) return(-1.0);	// something disastrous
    
    // now replace the active segment in Model by the rotated Id
    register unsigned int i, L=Betamask.on_no();
    
    if (Rmsup<=Rmsdown)	    // choose better fit
    {
	for (i=0; i<L; i++)
	    Model[i]=Idup[i];
	Model*=Rotup;
    }
    else
    {
	for (i=0; i<L; i++)
	    Model[i]=Iddown[i];
	Model*=Rotdown;
	Rmsup=Rmsdown;	// will be returned
    }
    Model+=Mctr;    // transpose back to original centroid
    
    Model.mask(Oldmask);	// reset original mask
    return(Rmsup);
}
// END of ideal_struct()

// ---- Handedness ----

/* check_torsion(): walks over each strand in 3D in Model and calculates all
 * the (i+1, i, k, m) torsion angles where k is i's partner in the next
 * strand, and m is i+1's partner. The angle should be negative.
 * Good and Bad will be set to the no. of correct
 * and incorrect torsion angles.
 * Return value: 1 if Good>=Bad, -1 if Good<Bad, 0 if not in 3D.
 */
int Beta_::check_torsion(Points_& Model, 
	unsigned int& Good, unsigned int& Bad) const
{
    Bits_ Oldmask=Model.mask(true);
    if (Model.dim()!=3)
    {
	Model.mask(Oldmask);
	return(0);
    }
    
    register unsigned int i, s, B, E;
    register int k, m;
    register double Tors;
    
    Good=Bad=0;
    for (s=0; s<strand_no()-1; s++)
    {
	B=Strands[s].beg(); E=Strands[s].end();
	for (i=B; i<E; i++)
	{
	    k=hbond_next(i); m=hbond_next(i+1);
	    if (k<0 || m<0) continue;	// no partner
	    
	    Tors=pos4_angle(Model[i+1], Model[i], Model[k], Model[m]);
	    if (Tors<-M_PI) continue;
	    
	    if (Tors<0.0) ++Good; else ++Bad;
	}
    }
    
    Model.mask(Oldmask);
    return((Good>=Bad)? 1: -1);
}
// END of check_torsion()

// ---- Input ----

/* >>: reads in a sheet description from the input stream In.
 * The format is:-
 * "SHEET\n"
 * "STRAND <beg> <end>\n"
 * "STRAND <beg> <end> [PAR|ANTI] <this> <other>\n"
 * .....
 * "END\n"
 * The residue position numbers are >=1.
 * Will not update the Beta_ object Beta if errors are detected.
 * The "fail" bit is set on error.
 */
istream& operator>>(istream& In, Beta_& Beta)
{
    if (!In) return(In);    // unhealthy stream, reject
    
    static const unsigned int BUFLEN=132;
    char Buf[BUFLEN+1];
    istrstream Istr(Buf, BUFLEN+1);
    
    // get the first line which should contain the word "SHEET" only
    while (In.getline(Buf, BUFLEN, '\n') && !strlen(Buf));
    if (!In || NULL==strstr(Buf, "SHEET"))
    {
	cerr<<"\n? >>Beta_: SHEET expected\n";
	In.clear(ios::failbit|In.rdstate());
	return(In);
    }
    
    // get the description of the first strand
    while (In.getline(Buf, BUFLEN, '\n') && !strlen(Buf));
    if (!In || strncmp(Buf, "STRAND", 6))
    {
	cerr<<"\n? >>Beta_: STRAND expected in line:\n"<<Buf<<endl;
	In.clear(ios::failbit|In.rdstate());
	return(In);
    }
    
    Istr.seekg(6);  // pos after the "STRAND" 
    int B=0, E=0;
    Istr>>B>>E;
    if (!Istr || B<=0 || E<=0)
    {
	cerr<<"\n? >>Beta_: Invalid limits in first STRAND: "<<B<<", "<<E<<endl;
	In.clear(ios::failbit|In.rdstate());
	return(In);
    }
    
    // so far, so good: start a temporary Beta_ object
    Beta_ Btemp(Strand_(B-1, E-1));
    
    // read all remaining strands and the "END\n" line
    char Pa[5];
    Strand_::Sense_ Sense;
    int T, O;
    
    while (In.getline(Buf, BUFLEN, '\n') && strncmp(Buf, "END", 3))
    {
	if (!strlen(Buf)) continue;	// skip "\n" lines
	if (strncmp(Buf, "STRAND", 6))
	{
	    cerr<<"\n? >>Beta_: STRAND expected in line:\n"<<Buf<<endl;
	    In.clear(ios::failbit|In.rdstate());
	    return(In);
	}
	
	Istr.seekg(6);  // pos after the "STRAND" 
	B=E=0;
	Istr>>B>>E;
	if (!Istr || B<=0 || E<=0)
	{
	    cerr<<"\n? >>Beta_: Invalid STRAND limits: "<<B<<", "<<E<<endl;
	    In.clear(ios::failbit|In.rdstate());
	    return(In);
	}
    
	// get the par/anti descriptor string
	Istr>>setw(5)>>Pa;
	Sense=Strand_::NONE;
	
	if (!strncmp(Pa, "PAR", 3)) Sense=Strand_::PAR;
	else if (!strncmp(Pa, "ANTI", 4)) Sense=Strand_::ANTI;
	
	if (!Istr || Sense==Strand_::NONE)
	{
	    cerr<<"\n? >>Beta_: [PAR|ANTI] expected in line:\n"<<Buf<<endl;
	    In.clear(ios::failbit|In.rdstate());
	    return(In);
	}
	
	// get the this/other descriptors
	T=O=0;
	Istr>>T>>O;
	if (!Istr || T<=0 || O<=0)
	{
	    cerr<<"\n? >>Beta_: Invalid this/other phase info: "<<T<<", "<<O<<endl;
	    In.clear(ios::failbit|In.rdstate());
	    return(In);
	}
	
	// attempt to add the new strand
	if (!(Btemp.add_strand(Strand_(B-1, E-1, Sense), T-1, O-1)))
	{
	    cerr<<"\n? >>Beta_: Invalid strand\n";
	    In.clear(ios::failbit|In.rdstate());
	    return(In);
	}
    }	    // while
    
    // single stranded sheets are disallowed
    if (Btemp.strand_no()<=1)
    {
	cerr<<"\n? >>Beta_: Sheets must have at least two strands\n";
	In.clear(ios::failbit|In.rdstate());
	return(In);
    }
    
    // check if the last line was "END". If not, warn but accept sheet
    if (strncmp(Buf, "END", 3))
	cerr<<"\n? >>Beta_: END expected\n";
    
    Beta=Btemp;	// accept new sheet
    return(In);
}
// END of >>

// ---- Output ----

/* write_to(): lists the calling object to the output stream Out.
 * The format is:-
 * "SHEET\n"
 * "STRAND <beg> <end>\n"
 * "STRAND <beg> <end> [PAR|ANTI] <this> <other>\n"
 * .....
 * "END\n"
 * The residue position numbers are >=1.
 * Does nothing if the sheet is empty. Protected
 */
void Beta_::write_to(ostream& Out) const
{
    unsigned int Sno=strand_no();
    if (!Sno) return;	// empty
    
    Out<<"SHEET\n";
    
    // first strand has no sense or registration
    Out<<"STRAND "<<(Strands[0].beg()+1)<<' '<<(Strands[0].end()+1)<<endl;
    
    // the rest of the strands have senses and registrations
    Strand_ Str;
    unsigned int i, R;
    int Rprev;
    
    for (i=1; i<Sno; i++)
    {
	Str=Strands[i];
	Out<<"STRAND "<<(Str.beg()+1)<<' '<<(Str.end()+1);
	Out<<((Str.sense()==Strand_::PAR)? " PAR ":" ANTI ");
	
	// locate a residue R on this strand that has a previous partner
	for (R=Str.beg(); R<=Str.end(); R++)
	{
	    if ((Rprev=hbond_prev(R))>=0)
		break;
	}
	Out<<(R+1)<<' '<<(Rprev+1)<<endl;
    }
    
    Out<<"END\n";   // terminate
}
// END of write_to()

// ==== END OF METHODS Secstr.c++ ====
