// ==== MEMBER FUNCTIONS Points.c++ ====

/* A maskable array of Vector_ objects for storing and
 * manipulating point coordinates. Essentially,  an improved
 * version of the Pset_ class.
 */

// SGI C++, 30. Apr. 1998. Andris Aszodi

// ---- CLASS HEADER ----

#include "Points.h"

// ==== Points_ MEMBER FUNCTIONS ====

// ---- Constructors ----

/* Initialises the array to hold N (default=1) vectors, each
 * D-dimensional (default=3). All vectors will be active.
 */
Points_::Points_(unsigned int N, unsigned int D): Maskarr_<Vector_>(N)
{
    for (unsigned int i=0; i<active_len(); i++) Idx[i]->dim(D);
}

/* Initialises the array with the bitmap in Initmask. It will have 
 * Initmask.len() items having D (default=3) dimension each.
 * Activation pattern is as in Initmask, of course.
 */
Points_::Points_(const Bits_& Initmask, unsigned int D): Maskarr_<Vector_>(Initmask)
{
    for (unsigned int i=0; i<active_len(); i++) Idx[i]->dim(D);
}

// ---- Dimension ----

/* dim_range(Low,High): returns the smallest and largest dimension
 * of the active vectors. Both are set to 0 if there are no active 
 * items.
 */
void Points_::dim_range(unsigned int& Low, unsigned int& High) const
{
    if (!active_len()) { Low=High=0; return; }
    
    unsigned int i, D;
    Low=UINT_MAX; High=0;
    for (i=0; i<active_len(); i++)
    {
	D=Idx[i]->dim();
	if (D<Low) Low=D;
	if (D>High) High=D;
    }
}
// END of dim_range()

/* dim_low(), dim_high(): return the lowest and highest dimension 
 * of the set of active points or 0 if ther were no active points.
 * (If both are needed, the dim_range() method above would be slightly
 * more efficient.)
 */
unsigned int Points_::dim_low() const
{
    if (!active_len()) return(0);
    
    unsigned int i, D, Low=UINT_MAX;
    for (i=0; i<active_len(); i++)
    {
	D=Idx[i]->dim();
	if (D<Low) Low=D;
    }
    return(Low);
}
// END of dim_low()

unsigned int Points_::dim_high() const
{
    if (!active_len()) return(0);
    
    unsigned int i, D, High=0;
    for (i=0; i<active_len(); i++)
    {
	D=Idx[i]->dim();
	if (D>High) High=D;
    }
    return(High);
}
// END of dim_high()

/* dim(): returns the dimension of the active vectors if they have the
 * same dimension or 0 if the dimensions are different or there are no
 * active vectors.
 * dim(D): sets the dimension of all active vectors to D. Returns old dimension
 * or 0 if they were different or there are no active vectors.
 */
unsigned int Points_::dim() const
{
    unsigned int Low, Up;
    dim_range(Low, Up);
    return(((!Low && !Up) || Low!=Up)? 0: Low);
}

unsigned int Points_::dim(unsigned int D)
{
    unsigned int i, Oldim=dim();
    for (i=0; i<active_len(); i++) Idx[i]->dim(D);
    return(Oldim);
}
// END of dim()

// ---- Arithmetics ----

/* Points*=Scalar: multiplies all active vectors by Scalar in place. 
 * Returns calling object.
 */
Points_& Points_::operator*=(double Scalar)
{
    for (unsigned int i=0; i<active_len(); i++) (*(Idx[i]))*=Scalar;
    return(*this);
}
// END of operator *= (scalar)

/* Points*=Matrix: premultiplies all active vectors by a square matrix
 * Matrix in place. Warnings are printed in case of dimension mismatches
 * and only those vectors which have matching dimensions will be modified.
 * Returns calling object.
 */
Points_& Points_::operator*=(const Sqmat_& Matrix)
{
    for (unsigned int i=0; i<active_len(); i++)
	*(Idx[i])=Matrix*(*(Idx[i]));
    return(*this);
}
// END of operator *=(Matrix)

/* Points+=Vector, Points-=Vector: adds/subtracts Vector to all active
 * vectors. If the dimensions don't match, then nothing happens (some
 * warnings are printed). Returns calling object.
 */
Points_& Points_::operator+=(const Vector_& Vector)
{
    for (unsigned int i=0; i<active_len(); i++)
	(*(Idx[i]))+=Vector;
    return(*this);
}

Points_& Points_::operator-=(const Vector_& Vector)
{
    for (unsigned int i=0; i<active_len(); i++)
	(*(Idx[i]))-=Vector;
    return(*this);
}

// END of operator +=,-= (Vector)

/* centroid(W): calculates the weighted centroid of the active points,
 * W should be at least as long as the number of currently active points.
 * centroid(): If the argument is omitted, then uniform weighting is assumed.  *  * Since the dimensions need not be equal, first the maximal dimension is obtained
 * and then all active vectors in the set are temporarily "promoted"
 * to this dimension by padding them with 0-s at the end (similar to the
 * Vector_::dim(D) method). If there were no active points then a 
 * warning is printed and a 3D null-vector returned. Otherwise the 
 * centroid is returned (with the maximal dimension).
 */
Vector_ Points_::centroid(const Vector_& W) const
{
    unsigned int Maxdim=dim_high();
    if (!Maxdim)
    {
	cerr<<"? Points_::centroid(W): No active points, default null-vector returned\n";
	return(Vector_(3));
    }
    
    unsigned int N=active_len();
    if (N>W.dim())
    {
	cerr<<"\n? Points_::centroid(W): Weight vector has too few elements ("
		<<W.dim()<<"<"<<N<<")\n";
	return(Vector_(3));
    }
    
    double Wsum=0.0;
    unsigned int i;
    for (i=0; i<N; i++)
    {
	if (W[i]<0.0)
	{
	    cerr<<"\n? Points_::centroid(W): W["<<i<<"]="<<W[i]<<", not allowed\n";
	    return(Vector_(3));
	}
	else Wsum+=W[i];
    }

    Vector_ Sum(Maxdim);
    Vector_ *Vec;
    unsigned int j, D;
    
    for (i=0; i<N; i++)
    {
	Vec=Idx[i]; D=Vec->dim();   // D<=Maxdim
	for (j=0; j<D; j++) Sum[j]+=W[i]*(*Vec)[j];
    }
    Sum/=Wsum;
    return(Sum);
}

Vector_ Points_::centroid() const
{
    unsigned int Maxdim=dim_high();
    if (!Maxdim)
    {
	cerr<<"? Points_::centroid(): No active points, default null-vector returned\n";
	return(Vector_(3));
    }
    
    Vector_ Sum(Maxdim);
    Vector_ *Vec;
    unsigned int i, j, D, N=active_len();
    
    for (i=0; i<N; i++)
    {
	Vec=Idx[i]; D=Vec->dim();   // D<=Maxdim
	for (j=0; j<D; j++) Sum[j]+=(*Vec)[j];
    }
    Sum/=N;
    return(Sum);
}
// END of centroid()

// ---- Distance matrices ----

/* dist_mat(Dist), dist_mat2(Dist2): construct the interpoint distance
 * and squared distance matrices using the active points only. The matrices
 * will not be touched if there are no active points or if the active
 * point dimensions don't match. The matrix size will be adjusted
 * silently if necessary. Return calling object.
 */
const Points_& Points_::dist_mat(Trimat_& Dist) const
{
    if (!dim())
    {
	cerr<<"? Points_::dist_mat(): No active points or dim mismatch within object\n";
	return(*this);
    }
    
    unsigned int i, j, N=active_len();
    Vector_ *Vecp;
    
    Dist.set_size(N);
    for (i=0; i<N; i++)
    {
	Vecp=Idx[i];
	Dist[i][i]=0.0;
	for (j=0; j<i; j++)
	    Dist[i][j]=diff_len(*Vecp, *(Idx[j]));
    }
    return(*this);
}
// END of dist_mat()

const Points_& Points_::dist_mat2(Trimat_& Dist2) const
{
    if (!dim())
    {
	cerr<<"? Points_::dist_mat2(): No active points or dim mismatch within object\n";
	return(*this);
    }
    
    unsigned int i, j, N=active_len();
    Vector_ *Vecp;
    
    Dist2.set_size(N);
    for (i=0; i<N; i++)
    {
	Vecp=Idx[i];
	Dist2[i][i]=0.0;
	for (j=0; j<i; j++)
	    Dist2[i][j]=diff_len2(*Vecp, *(Idx[j]));
    }
    return(*this);
}
// END of dist_mat2()

// ==== END OF MEMBER FUNCTIONS ====

// ---- Output ----

ostream& operator<<(ostream& Out, const Points_& Points)
{
    unsigned int i, N=Points.len(), Nact=Points.active_len();
    
    Out<<N<<" point";
    if (N!=1) Out<<'s';
    Out<<", "<<Nact<<" active\n"<<Points.mask();
    
    for (i=0; i<Nact; i++) Out<<Points[i];
    return(Out);
}
// END of operator <<

/* pdb_list(): produces a very crude PDB listing from Points provided it is
 * 3-dimensional. Amino acid code is 'G', one CA-atom per vector.
 * Return value: 0 if the dimension is outside the range [1..3], 
 * the actual dimension otherwise.
 */
int Points_::pdb_list(ostream& Out) const
{
    int D=dim();
    if (D<1 || D>3)
    {
	cerr<<"\n? Points_::pdb_list(): Data not 1D...3D\n";
	return(0);
    }
    
    for (unsigned int i=0; i<len(); i++)
    {
	Out<<"ATOM  "<<setw(5)<<(i+1)<<" "<<" CA  GLY  "<<setw(4)<<(i+1)<<"    ";
	Out.setf(ios::fixed, ios::floatfield);
	Out.setf(ios::showpoint);
	Out.precision(3);
	Out<<setw(8)<<(*this)[i][0]<<setw(8)
	    <<(D>1? (*this)[i][1]: 0.0)<<setw(8)
	    <<(D>2? (*this)[i][2]: 0.0);
	Out.precision(2);
	Out<<setw(6)<<1.0<<setw(6)<<1.0<<endl;
    }
    Out.width(0); Out.precision(6);
    return(D);
}
// END of pdb_list()

// ==== END OF FUNCTIONS Points.c++ ====
