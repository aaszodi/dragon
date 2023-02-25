// ==== PROJECT DRAGON: FUNCTIONS Steric.c++ ====

/* Steric adjustment routines. */

// SGI C++ 7.1, IRIX 6.2, 10. May 1998. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <iostream.h>
#include <iomanip.h>
#include <math.h>

// ---- UTILITY HEADERS ----

#include "Array.h"
#include "Vector.h"
#include "Bits.h"
#include "String.h"
#include "Hirot.h"

// ---- MODULE HEADER ----

#include "Steric.h"

// ---- FLOAT/DOUBLE ISSUES ----

/* NOTE: SGI provides single-precision floating point functions
 * such as sqrtf() etc. Some machines (SUNs in particular) don't
 * know about this: Define NO_MATHFLOATFUNC on the command line
 */
#ifdef NO_MATHFLOATFUNC
#define sqrtf sqrt
#endif

// ==== GLOBAL FUNCTIONS ====

/* update_bonddist(): re-calculates the first and second
 * squared neighbour distances in Model and puts them into the
 * corresponding off-diagonals of Dista.
 */
void update_bonddist(const Points_& Model, Trimat_& Dista)
{
    register unsigned int i, Rno=Dista.rno();
    
    for (i=1; i<Rno; i++)
	Dista[i][i-1]=diff_len2(Model[i], Model[i-1]);
    for (i=2; i<Rno; i++)
	Dista[i][i-2]=diff_len2(Model[i], Model[i-2]);
}
// END of update_bonddist()

// ==== Steric_ METHODS ====

// ---- Ideal distances ----

/* make_iddist(): if the Actual distance is outside the range
 * defined by Low<=Up (not tested), then it is mapped inside the range and
 * an UNsquared ideal distance is calculated. If Actual is close
 * to Up, then the ideal distance will be approx. Up-(Actual-Up), 
 * that is, reflected on Up, if Actual is infinite, it is mapped to Low.
 * Similarly, Actual values close to Low will be "reflected", 
 * while Actual==0 will be mapped to Up. If Actual is within the
 * range, then it is returned (no change). The distances are
 * all UNsquared. Static private
 */
inline
float Steric_::make_iddist(float Actual, float Low, float Up)
{
    if (Low<=Actual && Actual<=Up) return(Actual);	// within range
    if (Low==Up) return(Low);	// no range
    
    register float z, ul;
    if (Actual>Up)  // too large
    {
	z=Actual-Up; ul=Up-Low;
	z=ul*z/(ul+z);	// smooth transform: z is the "mirrored" difference
	return(Up-z);
    }
    else    // Actual<Low
    {
	z=Low-Actual;
	if (Up>=2.0F*Low) Up=1.99F*Low;	// curve turns up otherwise or 1st order
	ul=(Up-2.0F*Low)/(Low*Low);
	z+=Low+ul*z*z;
	return((z>=Up)? 0.99F*Up: z);
    }
}
// END of make_iddist()

/* limit_iddist(): in some adjustments ("midpoint kick", SCC-adjustments)
 * the C-alpha limits are not taken into account. This function makes sure
 * the unsquared Ideal distance falls between the limits defined by
 * Restraints for the (i, j)th C-alpha pair. No validity checks.
 * Return value: the modified ideal distance (now between the limits).
 * Static private
 */
inline
float Steric_::limit_iddist(float Ideal, const Restraints_& Restraints, 
	unsigned int i, unsigned int j)
{
    register float Low=Restraints.low(i, j), Up=Restraints.up(i, j);
    if (Ideal<Low) Ideal=Low;
    else if (Ideal>Up) Ideal=Up;
    return(Ideal);
}
// END of limit_iddist()

/* ideal_dist(): fills up the ideal distance matrix within the
 * calling object. Dista is the actual CA:CA distance matrix (squared), 
 * Fakebeta can be queried for sidechain centroid (SCC) distances, 
 * Restraints holds the external restraints and various bump lengths, 
 * Polymer supplies CA:SCC distances, Pieces provides the cluster
 * layout, Checkflags controls the adjustment.
 * If the SCORE flag was specified in Checkflags, then *Scores will contain
 * a score update on return, otherwise it is left unchanged.
 * If SCORE is set and Vl!=NULL then the violation list pointed to by Vl is made.
 */
void Steric_::ideal_dist(const Trimat_& Dista, const Fakebeta_& Fakebeta, 
	const Restraints_& Restraints, const Polymer_& Polymer, 
	const Pieces_& Pieces, int Checkflags, Scores_* Scores, Viollist_ *Vl)
{
    static const float BB_FAR=12.0F;	// max. CA dist for beta test
    static const float AB_FAR=9.0F;	// max. CA dist for alpha:beta test

    register float D2, Cad, Cad2, D, Id, Idb, 
	Bumpab, Bumpbb, Bumpbb2, Bmax, 
	Calow, Caup, Castrict, Maxstrict=-1.0;
    register unsigned int Rno=Polymer.len(), // no. of residues
	d, i, j, b;
    register int Cluno;
    Bits_ Lok(Rno+2);
    Viol_ Viol;
    
    // checkflag check
    if (!(Checkflags & ALL))   // no clus info, reset to ALL
    {
	cerr<<"\n? Steric_::ideal_dist(): No cluster check flags, ALL set\n";
	Checkflags|=ALL;
    }
    if (!(Checkflags & (RESTR | BOND)))   // no restraint choice, reset to RESTR
    {
	cerr<<"\n? Steric_::ideal_dist(): No restraint flags, RESTR set\n";
	Checkflags|=RESTR;
    }
    if (Vl!=NULL && !(Checkflags & SCORE))
    {
	cerr<<"\n? Steric_::ideal_dist(): Viollist ptr !=NULL, SCORE set\n";
	Checkflags|=SCORE;
    }
    if ((Checkflags & SCORE) && Scores==NULL)
    {
	cerr<<"\n? Steric_::ideal_dist(): Score ptr==NULL, SCORE cleared\n";
	Checkflags&=~SCORE;
    }
    
    Lastflags=Checkflags;    // save last adjustment type
    Checkflags&=ALL;	// use only the WITHIN/BETWEEN bits from now on
    
    /* Init the strictness and ideal distance matrices.
     * This depends on the adjustment required.
     */
    Idist.set_values(0.0);  // all
    Strimat.set_values(0.0);
    
    // init the BOND/NONBD/RESTR/SECSTR scores
    if (Lastflags & SCORE)
    {
	(*Scores)[Scores_::BOND].sum_reset();
	(*Scores)[Scores_::NONBD].sum_reset();
	(*Scores)[Scores_::RESTR].sum_reset();
	(*Scores)[Scores_::SECSTR].sum_reset();
    }
    
    /* init the Lok bit-vector (ON if 0<Lambda[i]<1),
     * the 0:th and Rno+1:th bits are always OFF (N/C terminal points)
     */
    for (i=1; i<=Rno; i++)
	Lok.set_bit(i, bool(Fakebeta.lambda(i)>0.0F && Fakebeta.lambda(i)<1.0F));
    
    // do external non-CA:CA restraints first
    if (Restraints.restr_no() && (Lastflags & REXT) && !(Lastflags & BOND))
    {
	// set up the "CA" const string for comparison
	static const String_ CA("CA");
	Clist1_<Restr_> Rlist(Restraints.ext_restr());  // iterator
	bool Liok=false, Ljok=false;
	
	/* check the external non-CA:CA dist restraints. Note that
	 * non-positive distances may occur here. These restraints
	 * will be taken into account only if there are no bump violations.
	 */
	for (Rlist.begin(); Rlist!=NULL; Rlist++)
	{
	    i=Rlist->pos(1); j=Rlist->pos(2);
	    
	    /* decide whether to update this particular pair:
	     * skip if the check is done within a cluster and (i, j) are
	     * not in it, or if the check is between different clusters
	     * and (i, j) happen to be in the same one.
	     */
	    if (Checkflags!=ALL)
	    {
		Cluno=Pieces.members(i, j);
		if (Checkflags==WITHIN && Cluno<0 || Checkflags==BETWEEN && Cluno>=0)
		    continue;
	    }
	    
	    // true if the i:th (j:th) SCC is distinct from the CA
	    Liok=bool(Fakebeta.lambda(i)>0.0F && Fakebeta.lambda(i)<1.0F);
	    Ljok=bool(Fakebeta.lambda(j)>0.0F && Fakebeta.lambda(j)<1.0F);
	    Cad2=Dista(i, j);
	    if (Rlist->atom(1)==CA)	// first atom is CA
	    {
		if (Rlist->atom(2)==CA)  // CA:CA restraint, will be dealt with later
		    continue;		    // and shouldn't be here anyway
		else    // CA:SCC restraint
		    D2=Ljok? Fakebeta.ab(i, j): Cad2;
	    }
	    else	// first atom must be "SCC", i.e. fake C-beta
	    {
		if (Rlist->atom(2)==CA)   // SCC:CA restraint
		    D2=Liok? Fakebeta.ab(j, i): Cad2;
		else	// SCC:SCC restraint
		    D2=(Liok && Ljok)? Fakebeta.bb(i, j): Cad2;
	    }
    
	    /* do not bother with wildly non-metric data: 
	     * experience has shown that D2==0.0 can wreak havoc, 
	     * so these distances are not adjusted
	     */
	    if (D2<=0.0) continue;
	    
	    // adjust if weight is not lower than previous
	    if (Rlist->strict()>=Strimat(i, j))
	    {
		Cad=sqrtf(Cad2);
		Strimat(i, j)=Rlist->strict();
		if (Maxstrict<Rlist->strict())
		    Maxstrict=Rlist->strict();
		if (D2<Rlist->low2() || D2>Rlist->up2())    // violation
		{
		    D=sqrtf(D2);
		    Id=make_iddist(D, Rlist->low(), Rlist->up());
		    Id*=Cad/D;	// scale to CA:CA
		    Idist(i, j)=limit_iddist(Id, Restraints, i, j); // limit to CA:CA
		    
		    if (Lastflags & SCORE)  // do the scoring
		    {
			(*Scores)[Scores_::RESTR]+=
			    Viol.rel_viol(D, Rlist->low(), Rlist->up(), Rlist->strict());
			if (Vl!=NULL)
			{
			    Viol.atom(1, Rlist->atom(1), i, Viol_::RESTR);
			    Viol.atom(2, Rlist->atom(2), j);
			    Vl->add_viol(Viol);
			}
		    }
		}
		else Idist(i, j)=Cad;	// keep actual with restr's weight
	    }
	}	// for Rlist
    }	    // if (there are external restraints)

    // scan all distances (or just 1st,2nd for BOND checks)
    Pieces_::Clutype_ Clutyp;
    int Dmax=(Lastflags & BOND)? 3: Rno+2;
    for (d=1; d<Dmax; d++) 
    {
	for (i=d; i<Rno+2; i++)
	{
	    j=i-d;
	    Cluno=Pieces.members(i, j);	// -1 if in separate clusters
	    Clutyp=Pieces.clu_type(Cluno);  // UNKNOWN if Cluno==-1
	    
	    /* decide whether to update this particular pair:
	     * skip if cluster membership is not what is required
	     * by Checkflags
	     */
	    if (Checkflags!=ALL)
	    {
		if (Checkflags==WITHIN && Cluno<0 || Checkflags==BETWEEN && Cluno>=0)
		    continue;
	    }
	    
	    // Check specific restraints?
	    if (!(Lastflags & BOND) && Restraints.specific(i, j))
	    {
		// same secondary structure but RINT is not set, skip
		if (!(Lastflags & RINT) &&
		    (Clutyp==Pieces_::HELIX || Clutyp==Pieces_::SHEET))
			continue;
		
		// external restraint (different clus or COIL) but REXT is not set, skip
		if (!(Lastflags & REXT) &&
		    (Clutyp==Pieces_::UNKNOWN || Clutyp==Pieces_::COIL))
			continue;
	    }
	    
	    // do not bother if this restraint is not strict enough
	    Castrict=Restraints.strict(i, j);
	    if (Castrict<Strimat[i][j]) continue;
	    
	    /* check minimal and maximal alpha-alpha separation:
	     * if violated, do not check for beta bumps.
	     * Squared values are used for the comparison.
	     */
	    Cad2=Dista[i][j]; Cad=sqrtf(Cad2);
	    Calow=Restraints.low(i, j);
	    Caup=Restraints.up(i, j);
		
	    // CA violation detected
	    if (Cad<Calow || Cad>Caup)
	    {
		// get ideal distance (unsquared)
		Idist[i][j]=make_iddist(Cad, Calow, Caup);

		// increase the strictness of CA:CA virtual bonds
		if (d<3)
		{
		    Viol.rel_viol(Cad, Calow, Caup, Castrict);
		    float Err=1.0F+Viol.rel_error();
		    Err*=Err; Err*=Err;    // 4th power...
		    Castrict*=Err;  // increase strictness
		}
		Strimat[i][j]=Castrict;	// set strictness
		if (Maxstrict<Castrict)
		    Maxstrict=Castrict;
		if (Lastflags & SCORE)  // do the scoring
		{
		    // get the score and violation typing
		    Scores_::Scotype_ Scotyp;
		    Viol_::Violtype_ Violtyp;
		    if (d<3)	// these are 1:2 and 1:3 virtual bonds
		    {
			Scotyp=Scores_::BOND;
			Violtyp=Viol_::BOND;
		    }
		    else if (Restraints.specific(i, j))	// external or secstr restraints
		    {
			// figure out if (i,j) are in the same secstr element
			if (Cluno==-1)	// not even the same segment, must be external
			{
			    Scotyp=Scores_::RESTR;
			    Violtyp=Viol_::RESTR;
			}
			else	// same segment: external or secstr restraint?
			{
			    switch(Clutyp)
			    {
				case Pieces_::HELIX:
				Scotyp=Scores_::SECSTR; Violtyp=Viol_::HELIX;
				break;
				case Pieces_::SHEET:
				Scotyp=Scores_::SECSTR; Violtyp=Viol_::SHEET;
				break;
				case Pieces_::COIL:  // external restraint
				case Pieces_::UNKNOWN:
				default:
				Scotyp=Scores_::RESTR; Violtyp=Viol_::RESTR;
			    }
			}
		    }
		    else    // general nonbond (bump etc.)
		    {
			Scotyp=Scores_::NONBD;
			Violtyp=Viol_::NONBD;
		    }
		    
		    // update the right score type
		    (*Scores)[Scotyp]+=
			Viol.rel_viol(Cad, Calow, Caup, Castrict);
		    
		    // note violation type and extent if required
		    if (Vl!=NULL)
		    {
			Viol.atom(1, "CA", i, Violtyp);
			Viol.atom(2, "CA", j);
			Vl->add_viol(Viol);
		    }
		}	// if (scoring)
		continue;   // no more adjustments
	    }
	    
	    /* Check if CA:CA virtual bond midpoints are too close.
	     * Normally, the CA:CA bumps should take care of this, 
	     * the check is done for pathological cases only where
	     * two virtual bonds "cross"
	     */
	    if (j && d>4 && (Lastflags & RINT) && Cad<AB_FAR &&
		!Restraints.hard(i, j) && !Restraints.hard(i-1, j) &&
		!Restraints.hard(i, j-1) && !Restraints.hard(i-1, j-1))
	    {
		register float Mid, Bump2;
		
		// squared bond midpoint distance
		Mid=(Dista[i][j]+Dista[i][j-1]+Dista[i-1][j]+Dista[i-1][j-1]
		    -Dista[i][i-1]-Dista[j][j-1])*0.25;
		Bump2=Restraints_::CA_BUMP; // just ONE Calpha radius away
		Bump2*=Bump2;
		if (Mid<Bump2)
		{
		    // give them a proper kick
		    register float Kick=sqrtf(Bump2/Mid), Newid;
		    Newid=Cad*Kick;
		    Idist[i][j]=limit_iddist(Newid, Restraints, i, j);
		    Newid=sqrtf(Dista[i-1][j-1])*Kick;
		    Idist[i-1][j-1]=limit_iddist(Newid, Restraints, i-1, j-1);
		    Newid=sqrtf(Dista[i-1][j])*Kick;
		    Idist[i-1][j]=limit_iddist(Newid, Restraints, i-1, j);
		    Newid=sqrtf(Dista[i][j-1])*Kick;
		    Idist[i][j-1]=limit_iddist(Newid, Restraints, i, j-1);
		    Strimat[i][j]=Strimat[i-1][j-1]=
			Strimat[i-1][j]=Strimat[i][j-1]=Restraints_::STRA;
			
		    if (Lastflags & SCORE)
			(*Scores)[Scores_::NONBD]+=
			    Viol.rel_viol(sqrtf(Mid), 2.0*Restraints_::CA_BUMP, 
				9999.9, Restraints_::STRA);
		    continue;	// no more adjustments here
		}
	    }
	    
	    // keep actual distance if not set otherwise
	    if (Strimat[i][j]==0.0)
	    {
		Idist[i][j]=Cad; 
		Strimat[i][j]=(d>=3)? 0.1: Castrict; // enforce 1st,2nd nb, lightweight otherwise
	    }
	    
	    /* don't bother if betas are too close in sequence,
	     * any previous violation was too strong or if
	     * the CAs are too far away
	     */
	    if (d<3 || Strimat[i][j]>Restraints_::STRB || Cad>BB_FAR)
		continue;
	    
	    /* the adjustment factor will be applied to alphas even
	     * for beta-violations. The displacements in these
	     * cases are scaled by the ratio of the beta-distances
	     * to the corresponding CA:CA distance. Violations belong
	     * to the NONBD category.
	     */
	    Idb=0.0; b=0; Castrict=Restraints_::STRB;
	    
	    // beta-beta check
	    if (Lok.get_bit(i) && Lok.get_bit(j))
	    {
		Bumpbb=Polymer.bumpb(i-1)+Polymer.bumpb(j-1);  // idx shifted 
		D2=Fakebeta.bb(i, j); Bumpbb2=Bumpbb*Bumpbb;
		Bmax=Caup+Fakebeta.ab(i, i)+Fakebeta.ab(j, j);
		
		if (D2>0.0 && (D2<Bumpbb2 || D2>Bmax*Bmax))  // violated
		{
		    D=sqrtf(D2);
		    Id=make_iddist(D, Bumpbb, Bmax);	// alpha
		    Idb+=Id; ++b;   // make up average ideal dist

		    if (Lastflags & SCORE)  // do the scoring
		    {
			(*Scores)[Scores_::NONBD]+=
			    Viol.rel_viol(D, Bumpbb, Bmax, Castrict);
			if (Vl!=NULL)
			{
			    Viol.atom(1, "SCC", i, Viol_::NONBD);
			    Viol.atom(2, "SCC", j);
			    Vl->add_viol(Viol);
			}
		    }
		}
	    }
	    
	    // alpha[i]-beta[j] check
	    if (!b && Cad<AB_FAR && Lok.get_bit(j))
	    {
		D2=Fakebeta.ab(i, j); 
		Bumpab=Polymer.bumpab(j-1);	// returns squared distlim
		Bmax=Caup+Fakebeta.ab(j, j);
		if (D2>0.0 && (D2<Bumpab || D2>Bmax*Bmax))
		{
		    D=sqrtf(D2); Bumpab=sqrtf(Bumpab);
		    Id=make_iddist(D, Bumpab, Bmax); // alpha now
		    Idb+=Id; ++b;
		    
		    if (Lastflags & SCORE)  // do the scoring
		    {
			(*Scores)[Scores_::NONBD]+=
			    Viol.rel_viol(D, Bumpab, Bmax, Castrict);
			if (Vl!=NULL)
			{
			    Viol.atom(1, "CA", i, Viol_::NONBD);
			    Viol.atom(2, "SCC", j);
			    Vl->add_viol(Viol);
			}
		    }
		}
	    }
	    
	    // beta[i]-alpha[j] check
	    if (!b && Cad<AB_FAR && Lok.get_bit(i))
	    {
		D2=Fakebeta.ab(j, i); 
		Bumpab=Polymer.bumpab(i-1);	// returns squared distlim
		Bmax=Caup+Fakebeta.ab(i, i);
		if (D2>0.0 && (D2<Bumpab || D2>Bmax*Bmax))
		{
		    D=sqrtf(D2); Bumpab=sqrtf(Bumpab);
		    Id=make_iddist(D, Bumpab, Bmax); // alpha now
		    Idb+=Id; ++b;

		    if (Lastflags & SCORE)  // do the scoring
		    {
			(*Scores)[Scores_::NONBD]+=
			    Viol.rel_viol(D, Bumpab, Bmax, Castrict);
			if (Vl!=NULL)
			{
			    Viol.atom(1, "SCC", i, Viol_::NONBD);
			    Viol.atom(2, "CA", j);
			    Vl->add_viol(Viol);
			}
		    }
		}
	    }
	    
	    // construct average ideal distance (for the alpha pair)
	    if (!b) continue;	// no beta-violations
	    Idist[i][j]=limit_iddist(Idb, Restraints, i, j);	// unsquared ideal "alpha" dist now
	    Strimat[i][j]=Castrict;	// with a moderate strictness
	    if (Maxstrict<Castrict)
		Maxstrict=Castrict;

	}	/* for i */
    }	    /* for d */
    
    // normalise strictness (largest is 1.0)
    if (Maxstrict>DBL_EPSILON) Strimat/=Maxstrict;
    
    // set up the spectral gradient
    if (Lastflags & SPECGRAD)
	Sp.weight(Strimat);
    
    // update score
    if (Lastflags & SCORE)
	Scores->update();
}
// END of ideal_dist()

// ---- Violation assessment ----

/* reset_viol(): sets the score normalisation factors (the appropriate
 * sums of various restraint weights) in Scores.
 * Call once before a simulation run.
 */
void Steric_::reset_viol(const Restraints_& Restraints, int Size, Scores_& Scores) const
{
    // bonds/bumps first
    Scores[Scores_::BOND].norm((Size-1)*Restraints_::STR1+(Size-2)*Restraints_::STR2);
    Scores[Scores_::NONBD].norm(double((Size-3)*(Size-4))/2.0*(Restraints_::STRA+Restraints_::STRB));
    
    // external restraints afterwards
    float Rwgt=0.0;
    if (Restraints.restr_no())
    {
	Clist1_<Restr_> Rlist(Restraints.ext_restr());  // iterator
	for (Rlist.begin(); Rlist!=NULL; Rlist++)
	    Rwgt+=Rlist->strict();
    }
    Scores[Scores_::RESTR].norm(Rwgt);
}
// END of reset_viol()

// ---- Adjustment ----

/* adjust_dist(): adjusts steric clashes in 'dist' space. Replaces the
 * entries in Dista by the corresponding entries in Idist
 * according to the cluster structure and the check choice.
 * The extent of the adjustment is controlled by Strict:
 * 0.0 means no adj, 1.0 total adj (default). 
 */
void Steric_::adjust_dist(Trimat_& Dista, const Pieces_& Pieces, 
    int Checkflags, float Strict) const
{
    if (Strict<=0.0) return;
    
    register unsigned int d, i, j, Rno=Dista.rno()-2;
    register float Str, D2;
    register int Cluno;
    
    // checkflag check
    if (!(Checkflags & ALL))   // no clus info, reset to ALL
    {
	cerr<<"\n? Steric_::adjust_dist(): No cluster check flags set, ALL assumed\n";
	Checkflags=ALL;
    }
    else Checkflags&=ALL;   // use the BETWEEN/WITHIN bits only
    
    // adjust all pairs within the current Check choice with >0.0 strictness
    for (d=1; d<Rno+2; d++)
	for (i=d; i<Rno+2; i++)
	{
	    j=i-d;
	    if (Checkflags!=ALL)
	    {
		Cluno=Pieces.members(i, j);
		if (Checkflags==WITHIN && Cluno<0 || Checkflags==BETWEEN && Cluno>=0)
		    continue;
	    }
	    Str=Strict*Strimat[i][j];
	    if (Str<=0.0) continue;	// no adjustment (zero strictness)
	    
	    D2=Idist[i][j]; D2*=D2;	// square ideal
	    if (Str>=1.0) Dista[i][j]=D2;	// full adjustment
	    else   // partial adjustment
		Dista[i][j]=(1.0-Str)*Dista[i][j]+Str*D2;
	}
}
// END of adjust_dist()

/* adjust_xyz(): tries to manoeuvre C-alphas to positions in Euclidean
 * space where bond, bump and separation constraints are not violated.
 * C-alpha coordinates are in Model. In the first overlaid version, 
 * Maxiter is the maximal number of spectral gradient iterations,
 * Eps is the relative stress change.
 * Returns the final stress or a negative float on error.
 * Also sets the Noconv int variable to non-0 if the Spectral Gradient
 * did not converge.
 * In the second overlaid version, the cluster layout is in Pieces
 * and Checkflags tells the routine what to update (uses Willie's
 * simple but efficient "pairwise displacement" optimisation).
 * No value returned.
 */
float Steric_::adjust_xyz(Points_& Model, int Maxiter, float Eps, int& Noconv)
{
    // sanity checks
    Noconv=1;
    if (!(Lastflags & (ALL | SPECGRAD)))
    {
	cerr<<"\n? Steric_::adjust_xyz(): Should have specified ALL|SPECGRAD in previous ideal_dist()\n";
	return(0.0);
    }
    if (!Maxiter)
    {
	cerr<<"\n? Steric_::adjust_xyz(SPECGRAD): No iteration number was specified, set to 10\n";
	Maxiter=10;
    }
    
    float Stress=0.0;
    int Iter=Maxiter;
    
    Stress=Sp.iterate(Idist, Model, Iter, Eps);
    Noconv=0;
    if (Stress<0.0)
    {
	cerr<<"\n? Steric_::adjust_xyz(SPECGRAD): Some horrible error occurred\n";
	return(Stress);
    }
    if (Iter<0)
    {
	cerr<<"\n Steric_::adjust_xyz(SPECGRAD): Don\'t panic :-)\n";
	Noconv=1;
    }
    return(Stress);
}

void Steric_::adjust_xyz(const Trimat_& Dista, Points_& Model,
    const Pieces_& Pieces, int Checkflags) const
{
    /* don't do anything if there's only 1 cluster and BETWEEN was prescribed */
    if (Pieces.clu_no()<=1 && (Checkflags & ALL)==BETWEEN) return;
    
    // store original Model mask and switch all vectors ON
    Bits_ Oldmask=Model.mask(true);
    
    register unsigned int d, i, j, Rno=Model.len()-2, Dim=Model.dim();
    if (!Dim)
    {
	cerr<<"\n? Steric_::adjust_xyz(): Dim mismatch among points\n";
	cerr<<"Len="<<(Rno+2)<<", mask is:\n"<<(Model.mask());
	    
	Model.mask(Oldmask);
	return;
    }

    // displacement vectors and adjustment weighting
    static Points_ Displ, Maxdispl, Newmodel;
    static Array_<float> Adjwgt, Maxdisplen2;
    
    Displ.len_dim(Rno+2, Dim); Maxdispl.len_dim(Rno+2, Dim);
    Newmodel.len_dim(Rno+2, Dim);
    for (i=0; i<Rno+2; i++)
    { Displ[i].set_values(); Maxdispl[i].set_values(); }
    Adjwgt.len(Rno+2); Adjwgt.set_values(0.0);
    Maxdisplen2.len(Rno+2); Maxdisplen2.set_values(0.0);
    
    int Cluno, Violno=0;
    register float Factor, Str, Dsplen2;
    Vector_ Half(Dim), Dvec(Dim);   // midpoint and current displacement
    
    // checkflag check
    if (!(Checkflags & ALL))   // no clus info, reset to ALL
    {
	cerr<<"\n? Steric_::adjust_xyz(): No cluster check flags set, ALL assumed\n";
	Checkflags|=ALL;
    }
    
    // scan distances, adjust violations 
    for (d=1; d<((Checkflags & BOND)? 3: Rno+2); d++)
    {
	for (i=d; i<Rno+2; i++)
	{
	    j=i-d;
	    
	    // check if i,j are good pairs for the adjustment
	    if ((Checkflags & ALL)!=ALL)
	    {
		Cluno=Pieces.members(i, j);
		if ((Checkflags & WITHIN) && Cluno<0 || (Checkflags & BETWEEN) && Cluno>=0)
		    continue;
	    }
	    
	    Str=Strimat[i][j];	// strictness 
	    if (Str<=0.0) continue; // no displacement
	    
	    Factor=(Dista[i][j]<DBL_EPSILON)? 10.0: Idist[i][j]/sqrtf(float(Dista[i][j]));
	    if (Factor<=0.0 || (Factor>0.99 && Factor<1.01))
		continue;	// no violation
	    
	    // limit extent of adjustment
	    if (Factor<0.1) Factor=0.1;
	    else if (Factor>10.0) Factor=10.0;

	    /* get weighted average displacement for each point
	     * and find the maximal displacement
	     */
	    Half=Model[i]; Half+=Model[j];
	    Half*=0.5;		// Half=(Model[i]+Model[j])/2.0;
	    Dvec=Model[i]-Half; 
	    Dvec*=Str*(Factor-1.0); // for weighting
	    
	    // store maximal displacement (premul by Str)
	    Dsplen2=Dvec.vec_len2();	// squared norm will do
	    if (Dsplen2>Maxdisplen2[i])
	    {
		Maxdispl[i]=Dvec; Maxdisplen2[i]=Dsplen2;
	    }
	    if (Dsplen2>Maxdisplen2[j])
	    {
		Maxdispl[j]=Dvec; Maxdisplen2[j]=Dsplen2;
	    }

	    // average displacements (weighted by strictness)
	    Displ[i]+=Dvec; Adjwgt[i]+=Str;
	    Displ[j]-=Dvec; Adjwgt[j]+=Str;
	    Violno++;
	}		/* for i */
    }	    /* for d */

    // everything was OK
    if (!Violno)
    {
	Model.mask(Oldmask);
	return;
    }
    
    // apply the displacements to the model
    for (i=0; i<Rno+2; i++)
    {
	if (Adjwgt[i]>DBL_EPSILON)  // Vector_ division has this limit, too
	{
	    Dsplen2=Displ[i].vec_len2();
	    if (25.0*Dsplen2<Maxdisplen2[i])	    // frustrated: avg. 5 times less than max
	    {
		Maxdispl[i]/=Adjwgt[i];
		Newmodel[i]=Model[i]+Maxdispl[i];   // use maximal displacement to "jump"
	    }
	    else
	    {
		Displ[i]/=Adjwgt[i];
		Newmodel[i]=Model[i]+Displ[i];	    // less frustrated, use avg. displacement
	    }
	}
	else Newmodel[i]=Model[i];
    }

    /* If the adjustment was to be done between clusters, then translate
     * and rotate the model clusters as rigid bodies towards Newmodel
     */
    if ((Checkflags & ALL) == BETWEEN)
    {
	Vector_ Mctr(Dim), Dctr(Dim), W(Rno+2);
	Hirot_ Hr;
	Bits_ Clumask;
	register unsigned int ci;
	bool Rotate;
	
	// adjust clusters (uses non-overlap+full-coverage implicitly)
	for (ci=0; ci<Pieces.clu_no(); ci++)
	{
	    Clumask=Pieces.clus(ci);
	    Rotate=bool(Clumask.on_no()>Dim);	// makes sense to rotate
	    Model.mask(Clumask); Newmodel.mask(Clumask);
	    
	    if (Rotate) // prepare weighted rotation
	    {
		// larger displacements have larger weight
		for (i=0; i<Clumask.on_no(); i++)
		    W[i]=0.01+Displ[i].vec_len2();
	    	Mctr=Model.centroid(W);
	    	Dctr=Newmodel.centroid(W);
	    }
	    else    // no weights, simple shift at end
	    {
	    	Mctr=Model.centroid();
	    	Dctr=Newmodel.centroid();
	    }
	    Model-=Mctr;
	    
	    if (Rotate)	// for clusters larger than the current simplex
	    {
		Newmodel-=Dctr;
		Hr.best_rot(Model, Newmodel, W);
		Model*=Hr.rot_matrix();
	    }
	    
	    // translate the model cluster to the new centroid
	    Model+=Dctr;
	}
    }
    else    // not BETWEEN: traditional non-rigid adjustment
	Model=Newmodel;

    Model.mask(Oldmask);
}
// END of adjust_xyz()

// ==== END OF FUNCTIONS Steric.c++ ====
