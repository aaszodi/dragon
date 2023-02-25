// ==== FUNCTIONS Hirot.c++ ====

/* Contains an implementation of McLachlan's RMS rotation algorithm
 * (aka "Procrustes rotation" in statistics)
 * adapted to high-dimensional problems. 
 */

// SGI C++ 4.0, IRIX 5.3, 29. May 1996. Andris

//---- STANDARD HEADERS ----

#include <math.h>

// ---- UTILITY HEADERS ----

#include "Hirot.h"
#include "Lu.h"

// ==== Hirot_ METHODS ====

// ---- Rotations ----

/* best_rot(): The RMS rotation algorithm adapted to
 * >=3D spaces. The point set X is rotated onto Y. Both are assumed to have
 * the same dimensionality and the same number of active points.
 * The activity masks need not be equal, but the first active
 * point in X will correspond to the first active in Y etc.
 * The active points are assumed to have been centered before the call.
 * W is a weight vector (assumed to be all-equal if omitted).
 * NOTE:  <3D point sets are not treated correctly.
 * Return value: the sign of the determinant of the transform matrix.
 */
int Hirot_::best_rot(const Points_& X, const Points_& Y, const Vector_& W)
{
    // checks
    int Actno=check_data(X, Y, W);
    if (Actno<=0)
    {
	cerr<<"\n? Hirot_::best_rot(X, Y, W): Inconsistent input data\n";
	return(-1);
    }
    
    make_mixtensor(X, Y, W);
    return(get_rot(X.dim()));
}
// END of best_rot()

int Hirot_::best_rot(const Points_& X, const Points_& Y)
{
    // checks
    int Actno=check_data(X, Y);
    if (Actno<=0)
    {
	cerr<<"\n? Hirot_::best_rot(X, Y): Inconsistent input data\n";
	return(-1);
    }
    
    make_mixtensor(X, Y);
    return(get_rot(X.dim()));
}
// END of best_rot()

/* best_rotflip(): the same as hi_rot() but there's no check for
 * the det sign. This routine returns a unitary transform matrix that
 * flips the conformations if necessary to achieve optimal superposition.
 */
void Hirot_::best_rotflip(const Points_& X, const Points_& Y, const Vector_& W)
{
    // checks
    int Actno=check_data(X, Y, W);
    if (Actno<=0)
    {
	cerr<<"\n? Hirot_::best_rotflip(X, Y, W): Inconsistent input data\n";
	return;
    }
    
    make_mixtensor(X, Y, W);
    get_rotflip(X.dim());
}

void Hirot_::best_rotflip(const Points_& X, const Points_& Y)
{
    // checks
    int Actno=check_data(X, Y);
    if (Actno<=0)
    {
	cerr<<"\n? Hirot_::best_rotflip(X, Y): Inconsistent input data\n";
	return;
    }
    
    make_mixtensor(X, Y);
    get_rotflip(X.dim());
}
// END of best_rotflip()

/* det_sign(): returns the sign of the determinant of the Mixtensor.
 * Returns 1 for "pure" rotations, -1 for "improper" rotations
 * (i. e. with inversion), 0 if the structures are "flat" and
 * the tensor is rank deficient. Conditions the SVD of Mixtensor
 */
int Hirot_::det_sign()
{
    if (Rank<0)	    // no decomp yet
    {
	cerr<<"\n? Hirot_::det_sign(): Not initialised\n";
	return(0);
    }
    
    // rank deficient?
    if (Rank<Mixtensor.rno()) return(0);
    
    // Check for improper rotations (the determinant of Mixtensor is -1)
    double Detu;
    Lu_ Lu(Mixtensor.rno());
    
    Lu.decomp(Mixtensor);
    Detu=Lu.det();
    return((Detu>0.0)? 1: -1);
}
// END of det_sign()

// ---- RMS and transformation ----

/* get_rms(): returns the weighted RMS value between the point sets
 * Y and Rot*X (Rot is in the calling object). 
 * Both sets should have the same no. of active points
 * and dimensions. W holds the weights for the pairs (if not uniform).
 * Returns -1.0 on error.
 */
double Hirot_::get_rms(const Points_& X, const Points_& Y, const Vector_& W) const
{
    register unsigned int Actno=X.active_len();
    
    if (!Actno || Actno!=Y.active_len() || W.dim()<Actno ||
	X.dim()!=Y.dim() || X.dim()!=Rot.rno())
    {
	cerr<<"\n? Hirot_::get_rms(X, Y, W): Dim problem\n";
	return(-1.0);
    }
    
    register double Err=0.0, Wsum=0.0;
    register unsigned int k;
    
    for (k=0; k<Actno; k++)
    {
	Err+=W[k]*diff_len2(Y[k], Rot*X[k]);	// weighted square dist
	Wsum+=W[k];
    }
    if (fabs(Wsum)<DBL_EPSILON)
	cerr<<"\n? Hirot_::get_rms(X, Y, W): W almost null-vector\n";
    else Err/=Wsum;
    return(sqrt(Err));
}

double Hirot_::get_rms(const Points_& X, const Points_& Y) const
{
    register unsigned int Actno=X.active_len();
    
    if (!Actno || Actno!=Y.active_len() || 
	X.dim()!=Y.dim() || X.dim()!=Rot.rno())
    {
	cerr<<"\n? Hirot_::get_rms(X, Y): Dim problem\n";
	return(-1.0);
    }
    
    register double Err=0.0;
    register unsigned int k;
    
    for (k=0; k<Actno; k++)
	Err+=diff_len2(Y[k], Rot*X[k]);	// square dist

    Err/=Actno;	    // sure Actno>0
    return(sqrt(Err));
}
// END of get_rms()

/* apply_transform(): replaces X with Rot*X and returns the
 * weighted RMS value between Y and Rot*X or -1.0 on error.
 * If W is not specified then uniform weights are used.
 * Rot is inside the calling object.
 */
double Hirot_::apply_transform(Points_& X, const Points_& Y, const Vector_& W) const
{
    register unsigned int Actno=X.active_len();
    
    if (!Actno || Actno!=Y.active_len() || W.dim()<Actno ||
	X.dim()!=Y.dim() || X.dim()!=Rot.rno())
    {
	cerr<<"\n? Hirot_::apply_transform(X, Y, W): Dim problem\n";
	return(-1.0);
    }
    
    register double Err=0.0, Wsum=0.0;
    register unsigned int k;
    
    X*=Rot;  // apply the transform matrix
    for (k=0; k<Actno; k++)
    {
	Err+=W[k]*diff_len2(Y[k], X[k]);	// weighted square dist
	Wsum+=W[k];
    }
    if (fabs(Wsum)<DBL_EPSILON)
	cerr<<"\n? Hirot_::apply_transform(X, Y, W): W vector almost null-vector\n";
    else Err/=Wsum;
    return(sqrt(Err));
}

double Hirot_::apply_transform(Points_& X, const Points_& Y) const
{
    register unsigned int Actno=X.active_len();
    
    if (!Actno || Actno!=Y.active_len() || 
	X.dim()!=Y.dim() || X.dim()!=Rot.rno())
    {
	cerr<<"\n? Hirot_::apply_transform(X, Y): Dim problem\n";
	return(-1.0);
    }
    
    register double Err=0.0;
    register unsigned int k;
    
    X*=Rot;  // apply the transform matrix
    for (k=0; k<Actno; k++)
	Err+=diff_len2(Y[k], X[k]);
    Err/=Actno;
    return(sqrt(Err));
}
// END of apply_transform()

// ---- Rotation matrix generation ----

/* get_rot(), get_rotflip(): finish up the job by generating the
 * transform matrix Rot from the SVD-decomposed Mixtensor.
 * get_rot() returns the sign of the determinant of Mixtensor
 * which is -1 for improper rotations ("flips") and 0 for 
 * rank-deficient cases ("flat rotations"). Both private
 */
int Hirot_::get_rot(unsigned int Dim)
{
    // do the SVD
    if (Svd.make_decomp(Mixtensor)) // error in decomposition
    {
	cerr<<"\n? Hirot_::get_rot(): Cannot decompose Mixtensor, no rotation\n";
	Rot.set_size(Dim);  // make it the Dim x Dim unit matrix
	Rot.diag_matrix();
	Rank=0;     // pretend a "flat" rotation
	return(1);
    }
    Rank=Svd.rank_cond();	// default rank conditioning
    
    /* NOTE: if the rank of Mixtensor is less than the dimension of the
     * point sets then the resulting "rotation" matrix will be
     * singular (det(Mixtensor)==0).
     */
    int Dsign=det_sign();
    register unsigned int i, j, k, Smpos;
    register double Temp;
    
    // find the smallest positive singular value if det(Mixtensor)==-1
    if (Dsign<0)
    {
	register double Sv;
	Temp=HUGE_VAL;
	for (i=0; i<Dim; i++)
	{
	    if ((Sv=Svd.w()[i])==0.0) continue;
	    if (Sv<Temp)
	    {
		Temp=Sv; Smpos=i;
	    }
	}
    }
    
    // generate the transform matrix in hyperspace
    Rot.set_size(Dim);
    for (i=0; i<Dim; i++)
	for (j=0; j<Dim; j++)
	{
	    Temp=0.0;
	    for (k=0; k<Dim; k++)
	    {
		Temp+=((Dsign<0 && k==Smpos)? -1: 1)*Svd.v()[i][k]*Svd.u()[j][k];
	    }
	    Rot[i][j]=Temp;
	}

    return(Dsign);
}
// END of get_rot()

void Hirot_::get_rotflip(unsigned int Dim)
{
    if (Svd.make_decomp(Mixtensor)) // error in decomposition
    {
	cerr<<"\n? Hirot_::get_rotflip(): Cannot decompose Mixtensor, no rotation\n";
	Rot.set_size(Dim);  // make it the Dim x Dim unit matrix
	Rot.diag_matrix();
	Rank=0;
	return;
    }
    Rank=Svd.rank_cond();

    /* generate the transform matrix in hyperspace */
    register double Temp;
    register unsigned int i, j, k;
    
    Rot.set_size(Dim);
    for (i=0; i<Dim; i++)
	for (j=0; j<Dim; j++)
	{
	    Temp=0.0;
	    for (k=0; k<Dim; k++) Temp+=Svd.v()[i][k]*Svd.u()[j][k];
	    Rot[i][j]=Temp;
	}
}
// END of get_rotflip()

// ---- Auxiliaries ----

/* check_data(): checks the integrity of the input data.
 * Returns -1 if the points have different dimensions, 
 * 0 if there are no active points in X, 
 * -2 if the no. of active points in X and Y don't match, 
 * -3 if W is shorter than the no. of active points (applies only
 * to the X, Y, W version).
 * Returns the no. of active points (>0) if OK. Private
 */
int Hirot_::check_data(const Points_& X, const Points_& Y)
{
    // fatal problems
    int Dim=X.dim();	// may be 0 if dims vary among actives
    if (!Dim || X.dim()!=Dim || Y.dim()!=Dim)
    {
	cerr<<"\n? Hirot_::check_data(): Point sets have different dimensions\n";
	return(-1);
    }
    unsigned int Actno=X.active_len();
    if (!Actno)
    {
	cerr<<"\n? Hirot_::check_data(): No active points in X\n";
	return(0);
    }
    if (Actno!=Y.active_len())
    {
	cerr<<"\n? Hirot_::check_data(): Active point no. mismatch, "
	    <<Actno<<"!="<<Y.active_len()<<endl;
	return(-2);
    }
    
    return(Actno);
}

int Hirot_::check_data(const Points_& X, const Points_& Y, const Vector_& W)
{
    int Actno=check_data(X, Y);
    if (Actno<=0) return(Actno);

    if (W.dim()<Actno)
    {
	cerr<<"\n? Hirot_::check_data(X, Y, W): Weight vector too short, "
	    <<W.dim()<<"<"<<Actno<<endl;
	return(-3);
    }
    return(Actno);
}
// END of check_data()

/* make_mixtensor(): constructs the DimxDim matrix 
 * which is at the heart of the method and hasn't got a name.
 * No dim checks are done. Private
 */
void Hirot_::make_mixtensor(const Points_& X, const Points_& Y, const Vector_& W)
{
    register unsigned int i, j, k, Dim=X.dim(), Actno=X.active_len();
    
    Mixtensor.set_size(Dim); Mixtensor.set_values(); // resize and zero
    for (k=0; k<Actno; k++)
    {
	for (i=0; i<Dim; i++)
	    for (j=0; j<Dim; j++)
		Mixtensor[i][j]+=(W[k]*Y[k][j]*X[k][i]);
    }
}

void Hirot_::make_mixtensor(const Points_& X, const Points_& Y)
{
    register unsigned int i, j, k, Dim=X.dim(), Actno=X.active_len();
    
    Mixtensor.set_size(Dim); Mixtensor.set_values(); // resize and zero
    for (k=0; k<Actno; k++)
    {
	for (i=0; i<Dim; i++)
	    for (j=0; j<Dim; j++)
		Mixtensor[i][j]+=(Y[k][j]*X[k][i]);
    }
}
// END of make_mixtensor()

// ==== END OF FUNCTIONS Hirot.c++ ====
