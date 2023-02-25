// ==== PROJECT DRAGON: METHODS Beta.c++ ====

/* Beta-sheet topology and geometry class. Cf. "Segment.h" and
 * "Sstrbase.h" for base class information.
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

#include "Beta.h"

// ---- DEFINITIONS ----

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

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

// ---- Virtual constructor ----

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
    
    // generate ideal UNsquared distances from the "Idup" coordinates
    Dist.set_size(L); Dist.set_values();
    Idup.mask(true);	// fully accessible, use real indices
    for (i=0; i<L; i++)
    {
	if (!member(i)) continue;
	for (j=0; j<=i; j++)
	{
	    if (!member(j)) continue;
	    Dist[i][j]=diff_len(Idup[i], Idup[j]);
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

/* ideal_dist(): puts the ideal beta-sheet UNsquared distances into
 * a distance matrix Dmat in the right position. Does nothing if
 * the sheet does not fit in the matrix. Prints a warning if Changed==true, 
 * since this indicates that the size was changed without updating the
 * ideal structure and therefore what is returned may be incorrect.
 * (The situation will be dealt with elegantly when the "mutable"
 * keyword finds its way into the compiler.) The ideal distances
 * will be applied at strictness Strict (from Sstrbase_) into Strimat.
 */
void Beta_::ideal_dist(Trimat_& Dmat, Trimat_& Strimat) const
{
    // check update status
    if (Changed)
    {
	cerr<<"\n? Beta_::ideal_dist(): make_idstruct() should have been called\n";
	return;
    }
    
    // size check
    if (Dmat.rno()<Dist.rno() || Strimat.rno()<Dist.rno())
    {
	cerr<<"\n? Beta_::ideal_dist(): Matrix too small\n";
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
	    if (Strimat[i][j]<=Strict)
	    {
		Dmat[i][j]=Dist[i][j];
		Strimat[i][j]=Strict;
	    }
	}
    }
}
// END of ideal_dist()

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
    Hirot_ Hr;
    
    // "up" phasing
    Hr.best_rot(Idup, Model);
    double Rmsup=Hr.get_rms(Idup, Model);
    if (Rmsup<0.0) return(-1.0);	// something disastrous
    Sqmat_ Rotup=Hr.rot_matrix();	// save the "up" matrix
    
    // "down" phasing
    Hr.best_rot(Iddown, Model);	    // the "down" matrix stays inside
    double Rmsdown=Hr.get_rms(Iddown, Model);
    if (Rmsdown<0.0) return(-1.0);
    
    /* now replace the active segment in Model by the rotated Id
     * taking the Strict into account
     */
    register unsigned int i, L=Betamask.on_no();
    register float Strict1=1.0-Strict;
    
    if (Rmsup<=Rmsdown)	    // choose better fit
    {
	for (i=0; i<L; i++)
	    Model[i]=Strict1*Model[i]+Strict*Rotup*Idup[i];
    }
    else
    {
	for (i=0; i<L; i++)
	    Model[i]=Strict1*Model[i]+Strict*Hr.rot_matrix()*Iddown[i];
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
 * "SHEET [strict]\n"
 * "STRAND <beg> <end>\n"
 * "STRAND <beg> <end> [PAR|ANTI] <this> <other>\n"
 * .....
 * "END\n"
 * The residue position numbers are >=1. The optional [strict]
 * parameter (a floating-point number between 0.0 and 1.0) controls
 * the strictness at which the ideal beta-structure should be 
 * applied. If missing, this strictness is 1.0.
 * Will not update the Beta_ object Beta if errors are detected.
 * The "fail" bit is set on error.
 */
istream& operator>>(istream& In, Beta_& Beta)
{
    if (!In) return(In);    // unhealthy stream, reject
    
    static const unsigned int BUFLEN=132;
    char Buf[BUFLEN+1];
    istrstream Istr(Buf, BUFLEN+1);
    
    // get the first line
    char *Sheetpos=NULL;
    while (In.getline(Buf, BUFLEN, '\n').good() && (!strlen(Buf) || Buf[0]=='#'));
    if (!In || NULL==(Sheetpos=strstr(Buf, "SHEET")))
    {
	cerr<<"\n? >>Beta_: SHEET expected\n";
	In.clear(ios::failbit|In.rdstate());
	return(In);
    }
    
    // try to get strictness if specified
    Istr.seekg(Sheetpos-Buf+5);	// pos after "SHEET"
    int Rdstate=Istr.rdstate();	// save old state
    float Str=1.0;  // default
    Istr>>Str;	    // try to read
    if (!Istr) Istr.clear(Rdstate); // wasn't there, restore previous state
    if (Str<=0.0)
    {
	cerr<<"\n? >>Beta_: Strictness "<<Str<<"<=0.0, sheet ignored\n";
	In.clear(ios::failbit|In.rdstate());
	return(In);
    }
    if (Str>1.0) Str=1.0;   // set silently
    
    // get the description of the first strand
    while (In.getline(Buf, BUFLEN, '\n').good() && (!strlen(Buf) || Buf[0]=='#'));
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
    Beta_ Btemp(Strand_(B, E));	    // [1..Rno] range
    Btemp.Strict=Str;
    
    // read all remaining strands and the "END\n" line
    char Pa[5];
    Strand_::Sense_ Sense;
    int T, O;
    
    while (In.good())
    {
        In.getline(Buf, BUFLEN, '\n');
	if (!strncmp(Buf, "END", 3)) break;   // no more
	if (!strlen(Buf) || Buf[0]=='#') continue;	// skip "\n" lines and comments
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
	if (!(Btemp.add_strand(Strand_(B, E, Sense), T, O)))
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
	cerr<<"\n? >>Beta_: please finish sheet description with \"END\" next time\n";
    
    Beta=Btemp;	// accept new sheet
    return(In);
}
// END of >>

// ---- Output ----

/* write_to(): lists the calling object to the output stream Out.
 * The format is:-
 * "SHEET [strict]\n"
 * "STRAND <beg> <end>\n"
 * "STRAND <beg> <end> [PAR|ANTI] <this> <other>\n"
 * .....
 * "END\n"
 * The residue position numbers are >=1. The optional strictness
 * value [strict] is written only if it is different from 1.0.
 * Does nothing if the sheet is empty. Protected
 */
void Beta_::write_to(ostream& Out) const
{
    unsigned int Sno=strand_no();
    if (!Sno) return;	// empty
    
    if (Strict==1.0)
	Out<<"SHEET\n";
    else
	Out<<"SHEET "<<Strict<<endl;
    
    // first strand has no sense or registration
    Out<<"STRAND "<<(Strands[0].beg())<<' '<<(Strands[0].end())<<endl;
    
    // the rest of the strands have senses and registrations
    Strand_ Str;
    unsigned int i, R;
    int Rprev;
    
    for (i=1; i<Sno; i++)
    {
	Str=Strands[i];
	Out<<"STRAND "<<(Str.beg())<<' '<<(Str.end());
	Out<<((Str.sense()==Strand_::PAR)? " PAR ":" ANTI ");
	
	// locate a residue R on this strand that has a previous partner
	for (R=Str.beg(); R<=Str.end(); R++)
	{
	    if ((Rprev=hbond_prev(R))>=0)
		break;
	}
	Out<<R<<' '<<Rprev<<endl;
    }
    
    Out<<"END\n";   // terminate
}
// END of write_to()

// ==== END OF METHODS Beta.c++ ====
