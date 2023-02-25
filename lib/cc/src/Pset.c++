// ==== FUNCTIONS Pset.c++ ====

/* The Point Set Class. A point set is an array of vectors
 * representing points in an Euclidean space. The points can be
 * switched on/off. Useful for storing atomic coordinates.
 */

// SGI C++ 3.2.1, IRIX 5.2, 19. Jan. 1995. Andris

// ---- HEADER ----

#include "Pset.h"

// ==== Pset_ MEMBER FUNCTIONS ====

// ---- Constructors ----

/* Init to contain N points (default 1) of D dimensions (default 3).
 * All points will be active.
 */
Pset_::Pset_(unsigned int N, unsigned int D):
    Points(N), Dim(D), Active(N, true)
{ for (unsigned int i=0; i<N; i++) Points[i].dim(D); }

/* Init by an NxD matrix Mat to hold N points of D dimensions.
 * The point coordinates are assumed to be held in the rows
 * of the matrix. All points will be active.
 */
Pset_::Pset_(const Matrix_& Mat)
{
    unsigned int i, N=Mat.rno();
    Points.len(N); Dim=Mat.cno();
    Active.len(N); Active.set_values(true);
    for (i=0; i<N; i++) Points[i]=Mat.row(i);
}

// ---- Matrix conversion ----

/* Matrix_(): converts a point set to an NxD rectangular matrix
 * with all points as rows.
 */
Pset_::operator Matrix_() const
{
    Matrix_ Mat(Points.len(), Dim);
    for (unsigned int i=0; i<Points.len(); i++)
	Mat.row(Points[i], i);
    return(Mat);
}

// ---- Access ----

/* len(Size): sets the number of points to Size. If Size==0
 * then no action is taken. The new points are switched OFF
 * by default. Returns old size
 */
unsigned int Pset_::len(unsigned int Size)
{
    unsigned int Osize=Points.len();
    if (!Size || Size==Osize) return(Osize);	// no action
    
    Points.len(Size); Active.len(Size);
    for (unsigned int i=Osize; i<Size; i++) // adjust vector dimensions
	Points[i].dim(Dim);		    // if the set grew
    return(Osize);
}
// END of len(Size)

/* dim(D): sets the dimension of ALL (not only of the active)
 * points to D. Returns old dimension. D==0 is silently ignored.
 */
unsigned int Pset_::dim(unsigned int D)
{
    if (!D || Dim==D) return(Dim);
    unsigned int Odim=Dim;
    for (unsigned int i=0; i<Points.len(); i++)
	Points[i].dim(D);
    Dim=D;
    return(Odim);
}
// END of dim(D)

// ---- Arithmetics ----

/* NOTE: All arithmetic methods operate on active points only. */

/* Pset*=S: scaling by the scalar S. Returns the calling object */
Pset_& Pset_::operator*=(double S)
{
    for (unsigned int i=0; i<Points.len(); i++)
    {
	if (!flag(i)) continue;	// skip inactives
	Points[i]*=S;
    }
    return(*this);
}
// END of Pset*=S

/* centroid(): calculates and returns the centroid of the active
 * points. All points have equal weights. If no active points were
 * found then the null-vector will be returned.
 */
Vector_ Pset_::centroid() const
{
    unsigned int N=Points.len(), Ano=N, i;
    Vector_ Ctr(Dim);
    for (i=0; i<N; i++)
    {
	if (!flag(i)) { Ano--; continue; }
	Ctr+=Points[i];
    }
    if (Ano) Ctr/=Ano;
    return(Ctr);
}
// END of centroid()

/* Pset+=Vec: translates the active points in Pset by a vector Vec.
 * Nothing happens if Vec has wrong dimension.
 * Returns the calling object.
 */
Pset_& Pset_::operator+=(const Vector_& Vec)
{
    if (Vec.dim()!=Dim)
    {
	cerr<<"? Pset_::Pset+=Vec: Dim mismatch\n"; return(*this);
    }
    for (unsigned int i=0; i<Points.len(); i++)
    {
	if (!flag(i)) continue;	// skip inactives
	Points[i]+=Vec;
    }
    return(*this);
}
// END of Pset+=Vec

/* Pset-=Vec: translates the active points in Pset by a vector -Vec.
 * (Centers the point set on Vec in other words.)
 * Returns the calling object.
 */
Pset_& Pset_::operator-=(const Vector_& Vec)
{
    if (Vec.dim()!=Dim)
    {
	cerr<<"? Pset_::Pset-=Vec: Dim mismatch\n"; return(*this);
    }
    for (unsigned int i=0; i<Points.len(); i++)
    {
	if (!flag(i)) continue;	// skip inactives
	Points[i]-=Vec;
    }
    return(*this);
}
// END of Pset-=Vec

/* Pset*=Sqmat: premultiplies the active points in Pset by a
 * square matrix. (It would be cumbersome to define a Sqmat*Pset
 * function instead.) Nothing happens if the dimensions don't match.
 * Returns calling object.
 */
Pset_& Pset_::operator*=(const Sqmat_& Sqmat)
{
    if (Sqmat.cno()!=Dim)
    {
	cerr<<"? Pset_::Pset*=Sqmat: Dim mismatch\n"; return(*this);
    }
    for (unsigned int i=0; i<Points.len(); i++)
    {
	if (!flag(i)) continue;	// skip inactives
	Points[i]=Sqmat*Points[i];
    }
    return(*this);
}
// END of Pset*=Sqmat

// ---- Distance matrices ----

/* dist_mat(),dist_mat2(): calculate ALL the interpoint distances and the
 * squared interpoint distances, respectively. Dmat is a Trimat_
 * the size of which will be adjusted silently to the number of 
 * ALL points. 
 */
void Pset_::dist_mat(Trimat_& Dmat) const
{
    // dim mismatch
    unsigned int Len=Points.len();
    if (Dmat.rno()!=Len) Dmat.set_size(Len);
    
    register unsigned int i, j;
    for (i=0; i<Len; i++)
    {
	Dmat[i][i]=0.0;	    // self-distance
	for (j=0; j<i; j++)
	    Dmat[i][j]=diff_len(Points[i], Points[j]);
    }
}
// END of dist_mat()

void Pset_::dist_mat2(Trimat_& Dmat) const
{
    // dim mismatch
    unsigned int Len=Points.len();
    if (Dmat.rno()!=Len) Dmat.set_size(Len);
    
    register unsigned int i, j;
    for (i=0; i<Len; i++)
    {
	Dmat[i][i]=0.0;	    // self-distance
	for (j=0; j<i; j++)
	    Dmat[i][j]=diff_len2(Points[i], Points[j]);
    }
}
// END of dist_mat2()

// ---- Output ----

ostream& operator<<(ostream& Out, const Pset_& Pset)
{
    Out<<Pset.len()<<"x"<<Pset.dim()<<" point set\n";
    Out<<Pset.active();
    for (unsigned int i=0; i<Pset.len(); i++) 
	if (Pset.flag(i)) Out<<Pset[i];
    return(Out);
}

// ==== END OF FUNCTIONS Pset.c++ ====
