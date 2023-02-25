// ==== FUNCTIONS Lu.c++ ====

/* LU-decomposition and linear equation solver
 * routines for square matrices.
 */

// SGI C++ 4.0, IRIX 5.3, 25. Oct. 1995. Andris

// ---- HEADER ---- 

#include "Lu.h"

// ---- DEFINITIONS ----

#ifdef FLT_MIN
#define LU_EPSILON (10.0*FLT_MIN)
#else
#define LU_EPSILON (1.0e-30)
#endif

// ==== Lu_ METHODS ====

// ---- Decomposition ----

/* decomp(): performs an LU-decomposition on the square matrix A.
 * Return value: the permutation sign (0 if A was singular).
 */
int Lu_::decomp(const Sqmat_& A)
{
    register unsigned int i, j, k, imax, n=A.rno();
    register double Large, Pivot, Tmp, Tmp2;
        
    /* get implicit scaling: if a row contains 0-s only,
     * then the matrix is singular which will be indicated
     * by setting Psign=0. Precision is controlled by
     * the constant LU_EPSILON (see Definitions above).
     */
    Array_<double> Scal(n);
    Psign=1;
    Lu=A;
    for (i=0; Psign && i<n; i++)
    {
	Large=0.0;
	for (j=0; j<n; j++)
	    if ((Tmp=fabs(Lu[i][j]))>Large) Large=Tmp;
	if (Large<LU_EPSILON)	// (almost) singular
	{ Psign=0; break; }
	Scal[i]=1.0/Large;
    }
    
    // loop over columns
    double *Tmpa=new double [n];    // C array for fast memcpy
    Perm.len(n);
    for (j=0; Psign && j<n; j++)
    {
	for (i=0; i<j; i++)
	{
	    Tmp=Lu[i][j];
	    for (k=0; k<i; k++) Tmp-=Lu[i][k]*Lu[k][j];
	    Lu[i][j]=Tmp;
	}
	
	// find largest pivot
	Large=0.0;
	for (i=j; i<n; i++)
	{
	    Tmp=Lu[i][j];
	    for (k=0; k<j; k++) Tmp-=Lu[i][k]*Lu[k][j];
	    Lu[i][j]=Tmp;
	    
	    if ((Tmp2=Scal[i]*fabs(Tmp))>=Large)    /* best so far */
	    { Large=Tmp2; imax=i; }
	}
	
	// interchange rows?
	if (j!=imax)
	{
	    memcpy(Tmpa, Lu[imax], n*sizeof(double));
	    memcpy(Lu[imax], Lu[j], n*sizeof(double));
	    memcpy(Lu[j], Tmpa, n*sizeof(double));
	    Psign*=(-1);    // parity change
	    Scal[imax]=Scal[j];
	}
	Perm[j]=imax;
	
	// get the pivot
	Pivot=Lu[j][j];
	if (fabs(Pivot)<LU_EPSILON)	// singularity
	{ Psign=0; break; }
	
	// divide by the pivot
	if (j<n-1)
	    for (i=j+1; i<n; i++) Lu[i][j]/=Pivot;
    }	// for j
    
    delete [] Tmpa; 
    return(Psign);
}
// END of decomp()

/* det: calculates the determinant from the decomposition results
 * stored inside. It is assumed that decomp() has been called beforehand.
 */
double Lu_::det() const
{
    if (!Psign) return(0.0);	// matrix is singular
    
    register unsigned int i;
    register double Det=Psign;
    
    // the "sum-of-logarithms" trick is NOT used here
    for (i=0; i<Lu.rno(); i++) Det*=Lu[i][i];
    return(Det);
}
// END of det()

// ---- Solution ----

/* solve(): solves the linear equation A*x=B and returns x.
 * The calling object is supposed to have been primed with
 * decomp(A) before this call. B is the "right-hand-side"
 * vector which is preserved. Checks for dim mismatches and returns
 * B untouched if not everything is OK.
 */
Vector_ Lu_::solve(const Vector_& B) const
{
    register unsigned int n=Lu.rno();

    // simple checks for dim mismatches
    if (n!=B.dim())
    {
	cerr<<"\n? Lu_::solve(): Dim mismatch\n";
	return(B);
    }
    
    register int i;
    register unsigned int j, ip;
    register double Tmp;
    Vector_ X(B);

    // permute forward
    for (i=0; i<n-1; i++)
	if ((ip=Perm[i])!=i)
	{ Tmp=X[ip]; X[ip]=X[i]; X[i]=Tmp; }
    
    // forward substitution
    for (i=0; i<n; i++)
    {
	Tmp=X[i];
	for (j=0; j<i; j++) Tmp-=Lu[i][j]*X[j];
	X[i]=Tmp;
    }
    
    // back substitution
    for (i=n-1; i>=0; i--)
    {
	Tmp=X[i];
	for (j=i+1; j<n; j++) Tmp-=Lu[i][j]*X[j];
	X[i]=Tmp/Lu[i][i];
    }
    return(X);
}
// END of solve()

/* lineq(): solves the nxn linear equation A*X=B and returns X.
 * Iterative improvement is done a la Num. Recipes (max. Maxit times, 
 * default 0).
 * This routine should be used if only one right-hand side vector B
 * is present. For several right-hand sides, use the decomposition
 * and solver routines separately. 
 * Return value: 1 means OK, 0 means A was singular, 
 * -1 means a dim mismatch occurred. 
 */
int Lu_::lineq(const Sqmat_& A, const Vector_& B, Vector_& X, unsigned int Maxit)
{
    // preparation
    unsigned int n=A.rno();
    if (n!=B.dim())
    {
	cerr<<"? lu_lineq(): Dim mismatch\n";
	return(-1);
    }
    
    // LU-decomposition
    decomp(A);
    if (!Psign)
    {
	cerr<<"? lineq(): matrix is singular\n";
	return(0);
    }
    
    // solution (sets the dimension of X)
    X=solve(B);
    
    // iterative improvement
    if (Maxit)
    {
	Vector_ Dx(n);
	for (unsigned int i=0; i<Maxit; i++)
	{
	    Dx=solve(A*X-B);
	    if (Dx.vec_len()<LU_EPSILON) break;
	    else X-=Dx;
	}
    }
    
    // OK, result has been calculated at last
    return(1);
}
// END of lineq()

#undef LU_EPSILON

// ==== END OF FUNCTIONS Lu.c++ ====

