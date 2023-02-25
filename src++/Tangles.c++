// ==== PROJECT DRAGON: METHODS Tangles.c++ ====

/* Secondary-structure-based tangle detection and elimination. */

// SGI C++, IRIX 6.2, 14. Aug. 1996. Andris Aszodi

// ---- MODULE HEADER ----

#include "Tangles.h"

// ==== METHODS ====

// ---- Update ----

/* update_pieces(): updates the internal data members so that tangle
 * detection can be carried out on the segments in Pieces now.
 * Must be called whenever Pieces changes (no automatic consistency yet).
 */
void Tangles_::update_pieces(const Pieces_& Pieces)
{
    register unsigned int Cluno=Pieces.clu_no();
    
    /* NOTE: the routine cannot detect tangles if there are no
     * secondary structures present in the segment list.
     */
    if (Cluno<=1 || Pieces.secs().len()==0)
    {
	cerr<<"\n? Tangles_::update_pieces(): Cannot do detangling\n";
	return;
    }
    
    Displ.len(Cluno); Displ.mask(true);
    Ctrs.len(Cluno); Ctrs.mask(true);
    Dnos.len(Cluno); Dnos.set_values(0);
    Tmask.len(Cluno); Tmask.set_values(false);
}
// END of update_pieces()

// ---- Tangle detection and elimination ----

/* tangle_detect(): checks whether the structure in Xyz is entangled,
 * given the segment layout in Pieces. (The calling object is assumed
 * to have called the update_pieces() method after the last change
 * to Pieces.)
 * Return value: 1 if there was(were) tangle(s), 0 if OK.
 */
unsigned int Tangles_::tangle_detect(const Pieces_& Pieces, Points_& Xyz)
{
    if (Pieces.clu_no()<=1 || Pieces.secs().len()==0) return(0);    // no check
    
    Bits_ Oldmask=Xyz.mask(true);
    unsigned int Violno=find_tangles(Pieces, Xyz);	// check only
    Xyz.mask(Oldmask);
    return(Violno);
}
// END of tangle_detect()

/* tangle_elim(): checks and optionally adjusts tangled conformations.
 * The check is carried out by testing whether parts of the chain
 * penetrate tetrahedra put on secondary structure segments.
 * The adjustment is carried out Iter times which
 * moves the entangled segments away from each other (scaled by Tadj).
 * Pieces holds the secondary structure and intervening coil list, 
 * Xyz holds the coordinates and its masking status will be restored
 * upon return.
 * Return value: the no. of entanglements (0 means OK). Also the
 * number of iterations actually done is returned in Iter.
 */
unsigned int Tangles_::tangle_elim(const Pieces_& Pieces, Points_& Xyz,
	double Tadj, unsigned int& Iter)
{
    if (Pieces.clu_no()<=1 || Pieces.secs().len()==0)
    {
	Iter=0; return(0);    // no check
    }

    Bits_ Oldmask=Xyz.mask(true);
    unsigned int Violno=0;
    register unsigned int i;
    
    unsigned int Dim=Xyz.dim();
    if (!Dim)
    {
	cerr<<"\n? Tangles_::tangle_elim(): Dim mismatch among points or no active points\n";
	Iter=0; return(0);
    }
    
    if (!Iter) Iter=1;	    // do at least one cycle
    Tmask.set_values(false);	// nobody is tangled yet: no centroids prepared
    Displ.dim(Dim); Ctrs.dim(Dim);  // prepare for adjustment
    
    for (i=0; i<Iter; i++)
    {
	Violno=find_tangles(Pieces, Xyz, 1);	    // detect and store elim info
	if (!Violno) break;	    // OK
	
	adjust_tangles(Pieces, Xyz, Tadj);
    }
    Iter=i;
    
    Xyz.mask(Oldmask);
    return(Violno);
}
// END of tangle_elim()

/* find_tangles(): attempts to detect tangles between the segments
 * stored in Pieces. The corresponding coordinates are in Points
 * which must be fully activated on entry. 
 * If Adjust==0 (default), then only the tangle detection is performed
 * and the violating pairs are not saved. In this case, 
 * find_tangles() returns right after the first tangle has been seen
 * and does not perform an exhaustive check.
 * Return value: the number of entanglements found. Private
 */
unsigned int Tangles_::find_tangles(const Pieces_& Pieces, Points_& Xyz, 
	int Adjust)
{
    Clist1_<Sstr_> Slist=Pieces.secs();    // secstr list iterator
    unsigned int Slen=Slist.len(), Cluno=Pieces.clu_no();	// no. of secstr items and clusters
    
    if (Adjust) Viols.clear();    // no entanglements yet
    
    /* Walk along the secstr list, check the other secstr segments
     * and the coils for entanglement. NOTE: The logic of this cycle is
     * based on the assumption that the segment masks for coils in Pieces
     * follow the secstr masks, and that the secstr masks are stored
     * in the same order as the secstr objects in the list returned
     * by Pieces.secs(). 
     */
    register unsigned int Si, Ti, Gi, Violno=0;
    Array_<Thidx_> Thedra;	// tetrahedron index array
    Violpair_ Vp;	    // temp for list growth
    Bits_ Clash(Cluno);	    // keep track of clashes
    Bits_ Sheetmask;	    // stores mask for current sheet (overlap check)
    bool Issheet;	// true if current secstr is sheet (for overlap checks)
    
    for (Slist.begin(), Si=0; Si<Slen; Si++, Slist++)	    // all secondary structure segments
    {
	// obtain array of tetrahedra superimposed on current secstr element
	Thedra=(*Slist)->get_thedra();
	if (Issheet=(*Slist)->is_beta())	// = intended
	    Sheetmask=(*Slist)->mask(Xyz.len());
	
	// collect clashes here for each segment (don't test pairs twice)
	Clash.set_values(false);
	
	for (Ti=0; Ti<Thedra.len(); Ti++)   // scan all tetrahedra
	{
	    // perform SVD on the current tetrahedron
	    if (make_thedron(Xyz, Thedra[Ti]))
		continue;   // cannot make it (flat?)
		
	    // check all segments from Si+1 onward
	    for (Gi=Si+1; Gi<Cluno; Gi++)
	    {
		if (Clash.get_bit(Gi))
		    continue;	// Si and Gi have already clashed before
		
		/* check if Gi overlaps with Si if Si is a sheet.
		 * Don't do tangle checks for overlapping sheets
		 */
		if (Issheet && 
		    (Sheetmask.on_no()+Pieces.clus(Gi).on_no())!=(Sheetmask | Pieces.clus(Gi)).on_no())
			continue;
		
		if (contain_segment(Pieces.clus(Gi), Xyz, Thedra[Ti].P1))
		{
		    // violation!
		    if (!Adjust)	// only test was requested
			return(1);	// report immediately w/o further checks
		    else
		    {
			Vp.Idx1=Gi; Vp.Idx2=Si;
			Viols+=Vp;  // add new violation pair to list
			Clash.set_bit(Gi, true);  // make a note of the clash
			Violno++;
		    }
		}	// if containment
	    }	    // for segments
	}	// for tetrahedra
    }	    // for secstr elements
    
    return(Violno);	// no containment, OK
}
// END of find_tangles()

/* adjust_tangles(): performs an adjustment in Euclidean space,
 * by moving entangled segments (the indices of which are listed 
 * in Viols) away from each other. (Only centroid translations are
 * performed.) Segment info is provided by Pieces, the movements
 * are controlled by Tadj>=0.0; Tadj==0.5 means that the segments
 * will be moved by 0.5 A. Private
 */
void Tangles_::adjust_tangles(const Pieces_& Pieces, Points_& Xyz, double Tadj)
{
    if (!Viols) return;	    // no violations
    
    unsigned int Cluno=Pieces.clu_no(), Dim=Xyz.dim();
    static Bits_ Oldmask;
    Oldmask=Xyz.mask();	    // update on each entry...
    
    // set flags in Vmask for all entangled clusters in this adjustment round
    Bits_ Vmask(Cluno); Vmask.set_values(false);
    for (Viols.begin(); Viols!=NULL; Viols++)
    {
	Vmask.set_bit(Viols->Idx1);
	Vmask.set_bit(Viols->Idx2);
    }
    
    /* NOTE: built-in masking is not used in accessing Displ and Ctrs:
     * for simplicity's sake we walk over Vmask and skip over
     * untangled clusters
     */
    
    /* Prepare the displacements and centroids. Tmask bits are true
     * for those segments which have been entangled in the current round
     * of adjustments and have valid centroids stored in Ctrs. Vmask
     * bits are true for recently entangled segments. A centroid is
     * calculated only for those segments which are true in Vmask but
     * false in Tmask (the adjustment moves the centroids so no
     * re-calculation is necessary within an adjustment cycle in 
     * tangle_elim()).
     */
    register unsigned int i, ix1, ix2;
    
    for (i=0; i<Cluno; i++)
    {
	if (!Vmask.get_bit(i)) continue;    // skip untangled
	
	Displ[i].set_values();	// zero displacements
	Xyz.mask(Pieces.clus(i));   // mask with i-th cluster (tangled)
	if (!Tmask.get_bit(i))
	    Ctrs[i]=Xyz.centroid(); // get i-th centroid (first entanglement seen)
    }
    Tmask|=Vmask;   // update flags for centroids

    // generate pairwise displacements
    static Vector_ H;	// pointing from Idx1:Idx2 midpoint to Idx1
    double Lh;		// the length of H
    
    Dnos.set_values(0);
    H.dim(Dim);
    
    // process all entangled pairs in Viols
    for (Viols.begin(); Viols!=NULL; Viols++)
    {
	H=Ctrs[ix1=Viols->Idx1];
	H-=Ctrs[ix2=Viols->Idx2];
	H/=2.0;
	Lh=H.vec_norm(); Lh+=Tadj;
	H*=Lh;
	Displ[ix1]+=H; Displ[ix2]-=H;	// move in opposite directions
	Dnos[ix1]++; Dnos[ix2]++;	// count viols per segment
    }
    
    // apply displacements to coords and centroids
    static Vector_ Adj;
    Adj.dim(Dim);
    for (i=0; i<Cluno; i++)
    {
	if (!Vmask.get_bit(i) || !Dnos[i])
	    continue;    // skip untangled or unadjusted
	    
	Adj=Displ[i]/Dnos[i];	// make adjustment vector
	Xyz.mask(Pieces.clus(i));   // apply to i-th cluster
	Xyz+=Adj;   // move coordinates
	Ctrs[i]+=Adj;	// move centroid
    }
    
    Xyz.mask(Oldmask);	// reset original mask
}
// END of adjust_tangles()

/* contain_segment(): checks whether a tetrahedron (decomposition in Thsvd)
 * contains a bit of the segment represented by its membership mask
 * Segmask. The coordinates in Xyz are supposed to have been masked
 * to "fully active" before the call. Oidx holds the index of the origin
 * of the tetrahedron.
 * Return value: 1 if containment was detected, 0 if not,  -1 on error.
 * Private
 */
int Tangles_::contain_segment(const Bits_& Segmask, const Points_& Xyz, unsigned int Oidx) const
{
    // a little check
    if (Xyz.len()!=Xyz.active_len())
    {
	cerr<<"\n? contain_segment(): Coords array not fully active\n";
	return(-1);
    }
    
    static Vector_ Sprev(4), Snext(4);
    register unsigned int k;
    int Start=0;
    
    /* The following cycle walks along the whole chain, and
     * skips discontinuities (useful for sheet clusters).
     */
    for (k=0; k<Xyz.len(); k++)
    {
	if (Segmask.get_bit(k))
	{
	    if (!Start)	    // start new contiguous region
	    {
		Start=1;
		make_svect(Xyz[k], Xyz[Oidx], Sprev);
		continue;
	    }
	    else
	    {
		make_svect(Xyz[k], Xyz[Oidx], Snext);
		if (th_viol(Sprev, Snext))
		    return(1);	// containment detected, no more needs to be done
		Sprev=Snext;	// otherwise, step along
	    }
	}
	else Start=0;	// reset
    }
    return(0);	// if we got here, then there's no containment
}
// END of contain_segment()

// ---- Tetrahedra ----

/* make_thedron(): constructs a tetrahedron out of the vectors
 * in Xyz indexed by the P1..P4 indices in Thidx, 
 * the origin is at Thidx.P1 .
 * Performs an SVD and puts the result in Thsvd.
 * Return value: 1 if an error occurred (iteration limit 
 * exceeded or the 4 vectors don't span a 3D space), 0 if
 * everything was fine. Private
 */
int Tangles_::make_thedron(const Points_& Xyz, const Thidx_& Thidx)
{
    unsigned int Dim=Xyz.dim();
    
    // paranoia
    if (Dim<3)
    {
	cerr<<"\n? make_thedron(): Dim: "<<Dim<<"<3\n";
	return(1);
    }
    
    static Matrix_ A(3); // always 3 columns: row no. varies
    static Vector_ Colvec;  // used as a temporary
    
    // init A: col vectors are [p1]-[p0], ...
    A.set_size(Dim, 3); Colvec.dim(Dim);
    Colvec=Xyz[Thidx.P2]; Colvec-=Xyz[Thidx.P1]; A.col(Colvec, 0);
    Colvec=Xyz[Thidx.P3]; Colvec-=Xyz[Thidx.P1]; A.col(Colvec, 1);
    Colvec=Xyz[Thidx.P4]; Colvec-=Xyz[Thidx.P1]; A.col(Colvec, 2);
    
    // do the singular value decomposition
    if (Thsvd.make_decomp(A))
	return(1);  // on error like iter limit exceeded
    
    if (3>Thsvd.rank_cond(SVD_EPSILON))
    {
	cerr<<"\n? make_thedron(): linear dependence\n";
	return(1);
    }
    
    return(0);	// OK
}
// END of make_thedron()

/* make_svect: the 4 coordinates of S will give the lin.comb.
 * factors with which the point Vec can be made up from the 4 position
 * vectors of a tetrahedron. Thsvd contains the SVD-decomposed
 * "tetrahedron matrix", Orig holds the coordinates of the
 * zeroth apex of the tetrahedron. Private
 */
void Tangles_::make_svect(const Vector_& Vec, const Vector_& Orig, Vector_& S) const
{
    static Vector_ Point, Sol3(3);
    
    Point=Vec; Point-=Orig; // Vec-Orig is the rhs. of the lincomb eqn.
    Sol3=Thsvd.lin_solve(Point);	// get the 3 indep. coords
    
    /* For traditional reasons, the 3 indep. coords go into the
     * coords 1..3 of S, and S[0] will be 1-(S1+S2+S3). This ordering
     * is probably not necessary...
     */
    register unsigned int i;
    register double S0=0.0;
    for (i=0; i<3; i++)
    {
	if (fabs(Sol3[i])<SVD_EPSILON) S[i+1]=0.0;
	else
	{
	    S[i+1]=Sol3[i];
	    S0+=Sol3[i];
	}
    }
    S[0]=1.0-S0;
}
// END of make_svect()

/* th_viol: Sprev and Snext are the lin.comb. coefficients
 * for the i-1 th (prev) and ith (next) points in a tetrahedron, respectively.
 * We want to know if there is a region in the segment between
 * prev and next which is contained in the tetrahedron. The
 * diagnostic for this condition is that some linear combination
 * of Sprev and Snext has all 4 coefficients in the range
 * 0.0 ... 1.0. Zmin and Zmax will bracket this region.
 * Return value: 1 if 0<Zmin<Zmax<1 (the prev--next segment
 * is at least partially contained by the tetrahedron), or
 * 0 if there is no containment. Private
 */
int Tangles_::th_viol(const Vector_& Sprev, const Vector_& Snext) const
{
    register double S12, Z0, Z1, Zlo, Zhi, Zmin, Zmax;
    register unsigned int i;
    
    Zmin=0.0; Zmax=1.0;
    for (i=0; i<4; i++)
    {
	if (Sprev[i]<0.0 && Snext[i]<0.0 ||
	    Sprev[i]>1.0 && Snext[i]>1.0)
		return(0);	// "same side out":-> out
	S12=Snext[i]-Sprev[i];
	if (fabs(S12)<SVD_EPSILON)
	{ Z0=0.0; Z1=1.0; }
	else
	{
	    Z0=-Sprev[i]/S12;
	    Z1=(1.0-Sprev[i])/S12; /* ?? used to be *S12 */
	}
	if (Z0<Z1) { Zlo=Z0; Zhi=Z1; }
	else { Zlo=Z1; Zhi=Z0; }
	if (Zlo>Zmin) Zmin=Zlo;
	if (Zhi<Zmax) Zmax=Zhi;
    }	    /* for i */

	return(0.0<=Zmin && Zmin<Zmax && Zmax<=1.0);
}
// END of th_viol()

// ==== END OF METHODS Tangles.c++ ====
