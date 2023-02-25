// ==== PROJECT DRAGON: METHODS Access.c++ ====

/* Conic accessibility calculations. */

// SGI C++, IRIX 6.2, 25. Nov. 1996. Andris Aszodi

// ---- STANDARD HEADERS ---- 

#include <math.h>
#include <string.h>
#include <fstream.h>
#include <strstream.h>

// ---- MODULE HEADER ----

#include "Access.h"

// ---- LOCAL DEFINITIONS AND CONSTANTS ----

#ifndef M_PI_2
#define M_PI_2		1.57079632679489661923	/* Pi/2 */
#endif

#ifndef DBL_EPSILON
#define EPSILON 1.0e-15
#else
#define EPSILON DBL_EPSILON
#endif

// ==== Access_ METHODS ====

// ---- Side chain accessibility ----

/* betacone_shield: calculates the local shieldedness values for 
 * each fake C-beta atom and stores them in the private Relsh array.
 * The internal array sizes are updated to
 * correspond to the size of Dista.
 * Note that currently there are no fake C-betas on the 0th and
 * (Rno+1)th "residues" (N/C-termini): for these, the C-alpha cones are used.
 * For the k-th point, all points which are closer than NBRADIUS are
 * selected, their distances from their common local centroid
 * is calculated and the angle of the smallest cone that encompasses the
 * whole set and centred on 'k' with an axis going through
 * the centroid is determined.
 */
void Access_::betacone_shield(const Trimat_& Dista, const Polymer_& Polymer) 
{
    // reset internal sizes
    unsigned int Rno;	// the real chain length
    set_size(Rno=Dista.rno()-2);    // Dista is larger (N/C-terminal points)
    
    // create the fake C-beta dist matrices 
    static Fakebeta_ Fakebeta;	// just resize below upon next call
    Fakebeta.update(Dista, Polymer);
    
    static const double NBRADIUS=8.0, NBRADIUS2=NBRADIUS*NBRADIUS;
    
    register unsigned int i, j, k, ci, cj, Closeno;
    register double Ang, Largang, D, Trisum, Isum;
    
    /* The "canonical ordering" here is that 0 is meaningless (N-terminus),
     * 1..Rno contains the BETAs, Rno+1 is the C-terminus, and Rno+1..2*Rno+3
     * the ALPHAs (where Rno+1 and 2*Rno+3 are the pseudo-alphas corresponding
     * to the NH3+ and COO- at the ends). This way the three dist matrices can be joined
     * into one big overall distmat without any clumsy indexing tricks.
     */
    
    /* scan all sidechain ("beta") points between 1..Rno. Note that
     * the Relsh[] array uses the [0..Rno-1] range for the shieldedness values.
     */
    for (k=1; k<=Rno; k++)
    {
	/* select close points: an index < Rno means the index-th beta,
	 * index>=Rno means the (index-Rno)-th alpha. Start at i==1
	 * because i==0 is the nonexistant N-terminal beta. The dist
	 * matrices are accessed via the "safe index" mechanism.
	 */
	Trisum=0.0;
	for (Closeno=0, i=1; i<=2*Rno+3; i++)
	{
	    // skip the current point and the nonexistent beta on the C-terminal moiety
	    if (i==k || i==Rno+1) continue;
	    
	    if (i<=Rno)	// beta:beta
		D=Fakebeta.bb(i, k);
	    else    // beta:alpha
		D=Fakebeta.ab(i-Rno-2, k);
	    if (D<0.0 || D>NBRADIUS2) continue;	    // too far away from k or non-metric
	    
	    Close[Closeno]=i;	// store index
	    Dik[Closeno++]=D;	// store D(i,k)^2
	    Trisum+=D;	// start summing for local centroid
	}

	if (Closeno<=1)	    // too few points, make it very exposed
	{
	    Relsh[k-1]=-1.0; continue;
	}
	
	/* calc the distances from the local centroid using
	 * Lagrange's Theorem. The first step is to sum all
	 * interpoint distances: the i-k distances were done
	 * in the previous cycle.
	 */
	for (i=0; i<Closeno; i++)
	{
	    ci=Close[i]; 
	    for (j=0; j<i; j++)
	    {
		cj=Close[j];
		if (ci<=Rno && cj<=Rno)	// beta-beta
		    D=Fakebeta.bb(ci, cj);
		else if (ci<=Rno && cj>=Rno+2) // beta-alpha
		    D=Fakebeta.ab(cj-Rno-2, ci);
		else if (ci>=Rno+2 && cj<=Rno) // alpha-beta
		    D=Fakebeta.ab(ci-Rno-2, cj);
		else    // alpha-alpha
		{
		    ci-=Rno+2; cj-=Rno+2;
		    D=Dista(ci, cj);
		    ci+=Rno+2;	// restore index
		}
		Trisum+=D;
	    }
	}
	Trisum/=(Closeno+1)*(Closeno+1);
    
	/* now get squared distances for the i-th point
	 * from the centroid. If the local dist set is
	 * non-metric enough then this dist may be negative;
	 * cheat by taking the abs value
	 */
	for (i=0; i<Closeno; i++)
	{
	    Isum=Dik[i]; ci=Close[i];
	    for (j=0; j<Closeno; j++)
	    {
		cj=Close[j];
		if (ci<=Rno && cj<=Rno)	/* beta-beta */
		    D=Fakebeta.bb(ci, cj);
		else if (ci<=Rno && cj>=Rno+2)
		    D=Fakebeta.ab(cj-Rno-2, ci);
		else if (ci>=Rno+2 && cj<=Rno)
		    D=Fakebeta.ab(ci-Rno-2, cj);
		else    // alpha-alpha
		{
		    ci-=Rno+2; cj-=Rno+2;
		    D=Dista(ci, cj);
		    ci+=Rno+2;	/* restore */
		}
		Isum+=D;
	    }
	    Di0[i]=fabs(Isum/(Closeno+1)-Trisum);
	}
    
	/* get the dist of the k-th point from the centroid
	 * in the same way and put into D
	 */
	D=0.0;
	for (i=0; i<Closeno; i++) D+=Dik[i];
	D=fabs(D/(Closeno+1)-Trisum);
	
	/* calc angle using the Cosine Rule for each entry in
	 * Close[] and determine the maximum. The angle is
	 * subtended by Dik and D, Di0 is the hypotenuse.
	 * If Di0->0, Ang->0; if either Dik or D->0, then
	 * we just jump (23-Aug-94)
	 */
	Largang=-1000.0;
	for (i=0; i<Closeno; i++)   // scan all angles
	{
	    if (D<EPSILON || Dik[i]<EPSILON) continue;	/* 23-Aug-94 */
	    Ang=(Dik[i]+D-Di0[i])/(2.0*sqrt(D*Dik[i]));
	    Ang=(fabs(Ang)<=1.0)? acos(Ang): Largang;	/* 14-Nov-94 */
	    if (Ang>Largang) Largang=Ang;
	}
	
	/* save the rel. shield of largest angle for the k-th point:
	 * non-metricities (and bugs) may lead to an Rsh=-6.38e+2
	 * value. In the past, these were replaced by 0.0, now
	 * get_shield() is instructed to ignore these.
	 * No warnings are printed (23-Aug-94)
	 */
	Relsh[k-1]=(Largang-M_PI_2)/M_PI_2;

    }	    // for k

}
// END of betacone_shield()

// ---- Accessibility scoring and adjustment ----

/* NOTE: distance "space" accessibility adjustment is NOT ported
 * from Version 3.x.
 */

/* solvent_xyz(): calculates the accessibility of the structure Xyz
 * given in Euclidean coordinates and performs an accessibility
 * adjustment. The sequence is taken from Polymer. Non-H-bonded
 * residues are moved outwards from the core so they need special 
 * treatment: Hbond is a bit-vector [0..Rno+1] in which every residue participating
 * in a H-bond is marked true (can be obtained from a Pieces_ object).
 * Dim mismatches are checked.
 * Return value: 0 on error, 1 if OK.
 */
int Access_::solvent_xyz(const Polymer_& Polymer, const Bits_& Hbond, 
	Points_& Xyz)
{
    // Hbond should be [0..Rno+1]
    unsigned int Rno=Polymer.len();
    if (Hbond.len()!=Rno+2)
    {
	cerr<<"\n? Access_::solvent_xyz(): Dim mismatch (Hbond:Points)\n";
	return(0);
    }
    
    // generate the shieldedness in private array
    if (!betacone_xyz(Polymer, Xyz))
	return(0);
    
    register unsigned int i;
    Shstate_ Shi;

    /* exposed residues:move "in",fact<1, buried:move "out",fact>1 */
    static const double ADJFACTORS[VERY_BURIED+1]=
	{0.90, 0.95, 0.99, 1.00, 1.01, 1.05, 1.10};

    /* Normal [0..Rno-1] residue indexing for everybody but Xyz.
     * N/C terminal pseudo-alphas are not adjusted
     */
    for (i=0; i<Rno; i++)
    {
	if (Surface.get_bit(i)) Shi=VERY_BURIED;
	else if (Buried.get_bit(i)) Shi=VERY_EXPOSED;
	else Shi=(Shstate_)get_shield(Relsh[i], Polymer.aa(i));

	/* move non-H-bonded buried residues outward,
	 * even if the shieldedness status was OK. The 0.40
	 * relsh limit is somewhat arbitrary
	 */
	if (Relsh[i]>=0.40 && !Hbond.get_bit(i+1))
	{
	    Xyz[i+1]*=ADJFACTORS[VERY_BURIED];
	    continue;
	}
	
	if (Shi!=AVERAGE)   // apply adjustment
	    Xyz[i+1]*=ADJFACTORS[Shi];
	
    }	    // for i
    
    return(1);
}
// END of solvent_xyz()

/* score_dist(),score_xyz(): calculate an accessibility score
 * (optimum 0.0) for distance spaces and Euclidean objects, 
 * respectively.
 * Return a very high value on error.
 */
float Access_::score_dist(const Polymer_& Polymer, const Trimat_& Dista)
{
    betacone_shield(Dista, Polymer);
    return(get_score(Polymer));
}
// END of score_dist()

float Access_::score_xyz(const Polymer_& Polymer, const Points_& Xyz)
{
    // generate the shieldedness in private array
    if (!betacone_xyz(Polymer, Xyz))
	return(1e10);
    return(get_score(Polymer));
}
// END of score_xyz()

/* betacone_xyz(): calculates the accessibility in Euclidean space
 * and stores the results in the internal array Relsh[].
 * Returns 0 on error, 1 if OK. Private
 */
int Access_::betacone_xyz(const Polymer_& Polymer, const Points_& Xyz)
{
    static Trimat_ Dista;   // re-allocated when necessary
    
    register unsigned int Rno=Polymer.len();
    
    // Xyz holds the extra N/C terminal moiety coordinates
    if (Polymer.len()!=Xyz.active_len()-2)
    {
	cerr<<"\n? Access_::betacone_xyz(): Dim mismatch (Polymer:Points)\n";
	return(0);
    }
    
    if (!Xyz.dim())
    {
	cerr<<"\n? Access_::betacone_xyz(): Dim mismatch within structure\n";
	return(0);
    }
    
    // generate the shieldedness in private array
    Xyz.dist_mat2(Dista);
    betacone_shield(Dista, Polymer);
    return(1);
}
// END of betacone_xyz()

/* get_score(): calculates the accessibility score. Private */
float Access_::get_score(const Polymer_& Polymer) const
{
    register unsigned int i, Rno=Polymer.len();
    Shstate_ Shi;
    float Score=0.0;

    /* artificial score values */
    static const double SCORES[VERY_BURIED+1]=
	{3.00, 1.00, 0.30, 0.00, 0.30, 1.00, 3.00};

    /* Normal [0..Rno-1] residue indexing for everybody but Xyz.
     * N/C terminal pseudo-alphas are not adjusted
     */
    for (i=0; i<Rno; i++)
    {
	Shi=(Shstate_)get_shield(Relsh[i], Polymer.aa(i));

	// prescribed accessibility: punish only if on wrong side
	if (Surface.get_bit(i) && Shi>=AVERAGE || 
	    Buried.get_bit(i) && Shi<=AVERAGE)
		{ Score+=SCORES[Shi]; continue; }
		
	// other residues
	Score+=SCORES[Shi];
	
    }	    // for i
    Score/=Rno;
    return(Score);
}
// END of get_score()

// ---- Size ----

/* set_size(): resets the size of the internal arrays to Rno,
 * the number of residues. Returns old size.
 * This routine must be called when the chain size changes.
 */
unsigned int Access_::set_size(unsigned int Rno)
{
    unsigned int Oldsize=Relsh.len();
    
    if (Oldsize!=Rno)
    {
	Relsh.len(Rno); Di0.len(2*(Rno+2)); 
	Dik.len(2*(Rno+2)); Close.len(2*(Rno+2));
	Surface.len(Rno); Buried.len(Rno);
	Surface.set_values(false); Buried.set_values(false);    // clear
    }
    return(Oldsize);
}
// END of set_size()

// ---- Shieldedness categories ----

/* get_shield: returns the shieldedness status of an amino
 * acid Aa with relative shieldedness Rsh. This function
 * determines whether Relsh means that this particular Aa
 * is exposed or buried or falls into the middle part of
 * its shieldedness distribution. (Cf. the Shstate_ type above.)
 * Private static
 */
int Access_::get_shield(double Rsh, char Aa)
{
    // 20 standard amino acids + B,Z
    static const int AANUM=22;
    static const char AAS[AANUM+1]="ABCDEFGHIKLMNPQRSTVWYZ";
    
    /* "experimental" burial limits. I have run a statistics
     * calculation on a set of non-homologous proteins and
     * these are the limits in the shieldedness distribution
     * for the amino acids.
     */
    static const double EXPBURLIMS[AANUM][VERY_BURIED]=
    {
	{-0.15,  -0.08,   0.00,  0.77,  0.81,  0.84}, /* A */
	{-0.15,  -0.09,  -0.07,  0.52,  0.56,  0.72}, /* B */
	{0.21,  0.31,  0.40,  0.84,  0.86,  0.89},	  /* C */
	{-0.27,  -0.21,  -0.16,  0.42,  0.50,  0.63}, /* D */
	{-0.31,  -0.25,  -0.20,  0.32,  0.42,  0.53}, /* E */
	{0.17,  0.26,  0.34,  0.80,  0.83,  0.87},	  /* F */
	{-0.12,  -0.06,  -0.01,  0.72,  0.77,  0.81}, /* G */
	{-0.18,  -0.10,  -0.02,  0.62,  0.70,  0.76}, /* H */
	{0.13,  0.23,  0.31,  0.83,  0.85,  0.89},    /* I */
	{-0.34,  -0.27,  -0.22,  0.22,  0.29,  0.38}, /* K */
	{0.12,  0.22,  0.32,  0.83,  0.85,  0.89},	  /* L */
	{0.00,  0.12,  0.23,  0.79,  0.83,  0.86},    /* M */
	{-0.25,  -0.18,  -0.13,  0.61,  0.70,  0.76}, /* N */
	{-0.22,  -0.15,  -0.09,  0.63,  0.70,  0.76}, /* P */
	{-0.26,  -0.20,  -0.15,  0.46,  0.55,  0.69}, /* Q */
	{-0.25,  -0.18,  -0.13,  0.44,  0.54,  0.64}, /* R */
	{-0.20,  -0.14,  -0.09,  0.70,  0.75,  0.80}, /* S */
	{-0.13,  -0.07,  -0.01,  0.70,  0.75,  0.80},  /* T */
	{0.09,  0.19,  0.28,  0.81,  0.84,  0.87},    /* V */
	{0.14,  0.24,  0.29,  0.82,  0.85,  0.90},    /* W */
	{-0.01,  0.11,  0.16,  0.72,  0.77,  0.81},   /* Y */
	{0.03,  0.04,  0.05,  0.48,  0.54,  0.57}     /* Z */
    };

    const char *Pos;
    register unsigned int Idx;
    
    // Rsh may be out of range->pretend it's OK (cf. betacone_shield())
    if (fabs(Rsh)>1.0) return(AVERAGE);
    
    // get the index of the Aa in AAS[]
    Pos=strchr(AAS, Aa);
    if (Pos==NULL) return(AVERAGE); // unknown AA, do nothing
    Idx=Pos-AAS;
    
    // get shieldedness state 
    int Sh; // assumes the Shield_ values are contiguous
    
    for (Sh=VERY_EXPOSED; Sh<AVERAGE; Sh++)   // exposed side
	if (Rsh<EXPBURLIMS[Idx][Sh])
	    return(Sh);
    for (Sh=VERY_BURIED; Sh>AVERAGE; Sh--)   // buried side
	if (Rsh>EXPBURLIMS[Idx][Sh-1])
	    return(Sh);
    return(AVERAGE);	// within the extremes
}
// END of get_shield()

// ---- Input/output ----

/* read_file(): reads known accessibilities from a file called Accfnm.
 * If it is NULL (the default) or "", then the surface/buried bit-vectors
 * are cleared. If Accfnm cannot be opened, then nothing is changed.
 * The file format is as follows:-
 * Lines starting with '#' are comments, 
 * [sS] <int> <int> ....\n are the nos. of residues known to be on the surface, 
 * [bB] <int> <int> ....\n are the nos. of residues known to be buried.
 * Any number of residues may be listed after the initial s, S, b, B character
 * provided the total length of the line does not exceed 132 characters,  
 * and any number of lines are allowed. The [sSbB] character must be the
 * first character on the line and the numbers are separated by whitespaces.
 * Duplicates are ignored. If a residue
 * is specified as surface and buried at the same time, it will be regarded
 * as an ordinary residue. Out-of-range residue numbers (outside 1..Rno)
 * are ignored with a warning.
 * Return value: 0 on error, 1 if OK.
 */
unsigned int Access_::read_file(const char *Accfnm)
{
    if (Accfnm==NULL || !strlen(Accfnm))    // just clear
    {
	Surface.set_values(false);
	Buried.set_values(false);
	return(0);
    }
    
    unsigned int Rno=Relsh.len();
    if (!Rno)	// = intentional
    {
	cerr<<"\n? Access_::read_file(): Please set size before input\n";
	return(0);
    }
    
    ifstream Inf(Accfnm);
    if (!Inf)
    {
	cerr<<"\n? Access_::read_file(\""<<Accfnm<<"\"): Cannot open\n";
	return(0);
    }
    Inf>>(*this);   // forgot this! Idiot! 1-Mar-96
    Inf.close();
    return(Inf.good() || Inf.eof());
}
// END of read_file()

/* >>: Reads the list of surface/buried residues from the stream Inf.
 * See also the comments for read_file().
 */
istream& operator>>(istream& Inf, Access_& Acc)
{
    if (!Inf)
    {
	cerr<<"\n? >>Access_: Cannot read from stream\n";
	return(Inf);
    }
    
    unsigned int Rno;
    if (!(Rno=Acc.Relsh.len()))	// = intentional
    {
	cerr<<"\n? >>Access_: Please set size before input\n";
	return(Inf);
    }
    
    static const unsigned int LINELEN=132;
    char Line[LINELEN+1];
    istrstream Instr(Line+1, LINELEN);    // buffer starts after 1st char!
    unsigned int Lineno;
    int Resno;
    bool S;
    
    // process file line-by-line
    Acc.Surface.set_values(false);
    Acc.Buried.set_values(false);    // reset
    for (Lineno=1; Inf.good(); Lineno++)
    {
        Line[0]='\0'; Inf.getline(Line, LINELEN, '\n');
	// skip empty lines, comments
	if (Line[0]=='#' || !strlen(Line)) continue;
	
	if (NULL==strchr("sSbB", Line[0]))
	{
	    cerr<<"\n? >>Access_: First char must be [#sSbB] in line "<<Lineno<<endl;
	    continue;
	}
	S=bool(Line[0]=='s' || Line[0]=='S');   // surface residues follow
	
	// process residue numbers in the line
	Instr.seekg(ios::beg); Instr.clear();
	while (Instr>>Resno)
	{
	    if (Instr.fail())
	    {
		cerr<<"\n? >>Access_: Malformed line ("<<Lineno<<")\n";
		break;
	    }
	    if (Resno<1 || Resno>Rno)	// check range
	    {
		cerr<<"\n? >>Access_: Residue no. "<<Resno
		    <<" is outside range [1.."<<Rno<<"] in line "<<Lineno<<endl;
		continue;
	    }
	    
	    // store residue no. in 0..Rno-1 range now
	    if (S) Acc.Surface.set_bit(Resno-1);
	    else Acc.Buried.set_bit(Resno-1);
	}
    }
    
    // check residues marked as surface and buried at the same time
    Bits_ Both=Acc.Surface & Acc.Buried;
    if (Both.on_no())
    {
	cerr<<"\n? >>Access_: Residues nonsensically specified as \"buried on surface\":\n";
	for (int i=0; i<Rno; i++)
	    if (Both.get_bit(i)) cerr<<(i+1)<<' ';
	cerr<<endl;
	Both.operator~();	    // negate in place (Sun CC could not deal with this)
	Acc.Surface&=Both; Acc.Buried&=Both;	// switch off ambiguous
    }
    return(Inf);
}
// END of >>

/* <<: lists the residues marked as surface/buried to Out nicely. */
ostream& operator<<(ostream& Out, const Access_& Access)
{
    if (!Access.Surface.on_no() && !Access.Buried.on_no())
    {
	Out<<"# No residues with known accessibilities\n";
	return(Out);
    }
    
    unsigned int i, k, Rno=Access.Relsh.len();
    Out<<"# List of residues with known accessibilities\n";
    if (Access.Surface.on_no())
    {
	Out<<"# Residues known to be on the surface";
	for (i=k=0; i<Rno; i++)
	{
	    if (!Access.Surface.get_bit(i)) continue;
	    if (k % 10 ==0) Out<<"\nS ";
	    Out<<(i+1)<<' '; k++;
	}
	Out<<endl;
    }
    
    if (Access.Buried.on_no())
    {
	Out<<"# Residues known to be buried";
	for (i=k=0; i<Rno; i++)
	{
	    if (!Access.Buried.get_bit(i)) continue;
	    if (k % 10 ==0) Out<<"\nB ";
	    Out<<(i+1)<<' '; k++;
	}
	Out<<endl;
    }
    return(Out);
}
// END of <<

// ==== END OF METHODS Access.c++ ====
