// ==== FUNCTIONS Vmutils.c++ ====

/* Various utilities for vectors and matrices which are
 * not member functions of the vector and matrix classes
 * but nevertheless deserve to be included in the VM library.
 */

// SGI C++ 3.2.1, IRIX 5.2, 21. Nov. 1994. Andris

// ---- HEADER ----

#include "Vmutils.h"

// ==== FUNCTIONS ====

// ---- Transpose products ----

/* trans_prod(): calculates the product X'X of a matrix X.
 * The result is a Trimat_ because it is symmetric. Safe
 * indexing is used which may incur a slight performance
 * penalty but it is necessary for generality here as X
 * is a base class ref which may be either a rectangular
 * Matrix_ or a Trimat_. See below for specialised versions.
 */
Trimat_ trans_prod(const Matbase_& X)
{
    unsigned int i, j, k, R=X.rno(), C=X.cno();
    Trimat_ Prod(C);
    double Temp;
    
    for (i=0; i<C; i++)
	for (j=0; j<=i; j++)
	{
	    Temp=0.0;
	    for (k=0; k<R; k++) Temp+=X(k,i)*X(k,j);
	    Prod[i][j]=Temp;
	}
    return(Prod);
}
// END of trans_prod()

/* trans_mprod(): calculates the product X'X of a rectangular Matrix_ X.
 * Result is a symmetric Trimat_. This routine should be slightly
 * faster for rectangulars than trans_prod() above but it cannot
 * be used for Trimat_s.
 */
Trimat_ trans_mprod(const Matrix_& X)
{
    unsigned int i, j, k, R=X.rno(), C=X.cno();
    Trimat_ Prod(C);
    const double *Xk;
    double Temp;
    
    for (i=0; i<C; i++)
	for (j=0; j<=i; j++)
	{
	    Temp=0.0; 
	    for (k=0; k<R; k++) { Xk=X[k]; Temp+=Xk[i]*Xk[j]; }
	    Prod[i][j]=Temp;
	}
    return(Prod);
}
// END of trans_mprod()

/* trans_wprod(): Calculates the product X'WX where X is a general
 * matrix and W is a diagonal matrix (represented by a Vector_).
 * Safe indexing is used (cf comments on trans_prod()). In case
 * of a dim mismatch between X and W then X'X is calculated.
 */
Trimat_ trans_wprod(const Matbase_& X, const Vector_& W)
{
    unsigned int i, j, k, R=X.rno(), C=X.cno();
    if (W.dim()!=R)
    {
	cerr<<"? trans_wprod(): X'WX dim mismatch, X'X returned\n";
	return(trans_prod(X));
    }
    
    Trimat_ Prod(C);
    double Temp;
    
    for (i=0; i<C; i++)
	for (j=0; j<=i; j++)
	{
	    Temp=0.0;
	    for (k=0; k<R; k++) Temp+=W[k]*X(k,i)*X(k,j);
	    Prod[i][j]=Temp;
	}
    return(Prod);
}
// END of trans_wprod()

/* trans_mwprod(): calculates the product X'WX of a rectangular matrix
 * X and a diagonal matrix represented by the Vector_ W. Less general
 * but slightly faster than trans_wprod() above.
 */
Trimat_ trans_mwprod(const Matrix_& X, const Vector_& W)
{
    unsigned int i, j, k, R=X.rno(), C=X.cno();
    if (W.dim()!=R)
    {
	cerr<<"? trans_mwprod(): X'WX dim mismatch, X'X returned\n";
	return(trans_mprod(X));
    }
    
    Trimat_ Prod(C);
    const double *Xk;
    double Temp;
    
    for (i=0; i<C; i++)
	for (j=0; j<=i; j++)
	{
	    Temp=0.0; 
	    for (k=0; k<R; k++) { Xk=X[k]; Temp+=W[k]*Xk[i]*Xk[j]; }
	    Prod[i][j]=Temp;
	}
    return(Prod);
}
// END of trans_mwprod()

// ==== END OF FUNCTIONS Vmutils.c++ ====

