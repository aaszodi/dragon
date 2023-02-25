// ==== FUNCTIONS Ql.c++ ====

/* Diagonalisation of symmetric real matrices (Trimat_ class).
 * Based on Numerical Recipes, with some modifications and
 * C++ -isation. Implements the QL-algorithm with implicit shifts.
 * Also contains a simple iteration by "popular demand".
 */

// SGI C++ 4.0, IRIX 5.3, 2. May 1996. Andris

// ---- HEADER ----

#include "Ql.h"

// ---- DEFINITIONS ----

#ifdef FLT_MIN
#define QL_EPSILON (10.0*FLT_MIN)
#else
#define QL_EPSILON (1.0e-30)
#endif
#define RND0(x) (x=(fabs(x)<QL_EPSILON)? 0.0: (x))

// ---- TYPES ---- 

/* type of function returning int: required by qsort() */
typedef int (*Compfnc_)(const void*,const void*);

/* type for eigenvalue sorting */
typedef struct
{
    double Eig;	/* an eigenvalue along with its */
    int Idx;	/* original position */
} Eigidx_ ;

// ---- PROTOTYPES ----

/* eig_cmp: compares two Eigidx structs for qsort(). This is in C */
static int eig_cmp(const Eigidx_ *E1, const Eigidx_ *E2);

/* tred2(): Housholder tridiagonalisation. a is a real, symmetric matrix
 * (size n*n). The main diagonal of the tridiag. output is returned in d,
 * the second diagonal in e with e[1]==0.0. On return, a contains the
 * transformation matrix "Q". 'd' and 'e' are assumed to have the
 * correct lengths adjusted before the call.
 */
static void tred2(Sqmat_& a, double *d, double *e);

/* tqli(): QL algorithm with implicit shifts on tridiagonal matrices.
 * The main diagonal is in d, the second diagonal is in e, with e[1]
 * ignored. If d and e were obtained from a general symmetric
 * real matrix by Housholder transformation by tred2(), then z should
 * contain the square transformation matrix "Q" on input; otherwise it should
 * be the unit matrix. On output, d contains the eigenvalues and z the
 * eigenvectors, with the k-th _column_ corresponding to the k-th eigenvalue.
 * Itno supplies the maximum allowable no. of iterations (an addition
 * by A.A.). Return value: 0 if OK, 1 if iteration limit has been
 * exceeded (added by A.A. to replace the nrerror() error message function
 * originally used in Numerical Recipes routines).
 */
static int tqli(double *d, double *e, Sqmat_& z, int Itno);

// ==== FUNCTIONS ====

// ---- WRAPPER ----

/* eigen_ql(): a 'wrap' function driving the Housholder and QL routines.
 * Takes a lower triangular matrix Mat as input (preserved) and
 * produces the eigenvalues in Eval and the eigenvectors in Evec.
 * The sizes of these are adjusted silently if necessary.
 * Eigenvalues are sorted in decreasing order and the corresponding
 * eigenvectors are the _col_vectors_ of Evec.
 * Index shifts are performed to hack around the Fortranese
 * [1..N] convention of Numerical Recipes.
 * Return value: 0 if OK, 1 if iteration limit was exceeded in tqli().
 */ 
int eigen_ql(const Trimat_& Mat, Vector_& Eval, Sqmat_& Evec)
{
    const int ITERNO=30;    // max. number of iterations
    unsigned int Size=Mat.rno();

    /* Size checks. Eval and Evec may not have the right dimensions,
     * if this happens then the size is adjusted.
     */
    Eval.dim(Size); Evec.set_size(Size);
    
    /* Qmat is a full square matrix initialised with Mat and will
     * contain the eigenvectors in the end as _columns_
     */
    Sqmat_ Qmat(Mat);	// tri->sq conversion init 
    Qmat.ftn_idx();	// Fortranese [1..N] index shift

    /* create conventional arrays to hold the main and second
     * diagonals of the tridiagonalised Qmat, also perform
     * index shifting.
     */
    double *Diag1, *Diag2;
    Diag1=new double [Size]; Diag1--;	// alloc and FORTRAN shift
    Diag2=new double [Size]; Diag2--;

    tred2(Qmat,Diag1,Diag2);	// Housholder tridiagonalisation

    /* apply shifted QL-transforms: get eigenvalues in Diag1,
     * eigenvectors in Qmat's _columns_. No sorting yet
     */
    int Err=tqli(Diag1,Diag2,Qmat,ITERNO);
    if (Err)
	cerr<<"? eigen_ql(): Iteration limit ("<<ITERNO<<") exceeded\n";

    /* shift back addresses to C-style indexing */
    Qmat.c_idx();
    Diag1++; Diag2++;
    delete [] Diag2;	// not needed any more

    /* Sort eigenvalues in decreasing order. The unsorted
     * eigenvalues are copied from Diag1 into a Size-long
     * conventional array of Eigidx_ structures, the sorting
     * is done via qsort() from the C standard library
     */
    Eigidx_ *Evs=new Eigidx_ [Size];
    unsigned int i;
    for (i=0; i<Size; i++)
    {
	Evs[i].Eig=RND0(Diag1[i]);  // small eigenvalues are rounded to 0.0
	Evs[i].Idx=i;	// store index as well
    }
    delete [] Diag1;
    qsort(Evs,Size,sizeof(Eigidx_),(Compfnc_)eig_cmp);	// C stdlib quicksort
    
    /* permute eigenvectors according to the new order
     * of the eigenvalues and copy them into the _cols_
     * of the Evec matrix. Also copy eigenvalues to Eval
     */
    unsigned int j, k;
    for (i=0; i<Size; i++)
    {
	Eval[i]=Evs[i].Eig;	// copy back eigenvalues
	k=Evs[i].Idx;
	for (j=0; j<Size; j++)
	    Evec[j][i]=Qmat[j][k];	/* copy eigenvectors */
    }

    /* output and cleanup */
    delete [] Evs;
    return(Err);
}
// END of eigen_ql()

/* ---- AUXILIARY ROUTINES TO eigen_ql() ---- */

/* eig_cmp(): compares two Eigidx structs for qsort(). This is in C */
static int eig_cmp(const Eigidx_ *E1, const Eigidx_ *E2)
{
    return((E2->Eig>E1->Eig)? 1: (E2->Eig<E1->Eig)? -1: 0);
}
// END of eig_cmp()

/* ---- LARGE EIGENVALUE/EIGENVECTOR METHODS ---- */

/* eigen_positer(): tries to find the Poseno largest positive eigenvalues
 * and corresponding eigenvectors of Mat. Mat, Eval, Evec are
 * as in eigen_ql(). Poseno<=size of Mat. The actual number
 * of positive eigenvalues are returned.
 */
int eigen_positer(int Poseno, const Trimat_& Mat, Vector_& Eval, Sqmat_& Evec)
{
    const double EPS=1e-6;
    const int MAXITERNO=100;

    int Evalno, Iterno, i, j, Posevalno=0, Size=Mat.rno();
    if (Poseno>Size) Poseno=Size;
    
    Trimat_ Matrix(Mat);    // original preserved
    Vector_ Vec(Size), Oldvec(Size);
    double Ev, Oldev;
    
    Eval.set_values(); Evec.set_values();   // zero output
    for (Evalno=0; Evalno<Size; Evalno++)
    {
	// init as random unit vector
	for (i=0; i<Size; i++)
	    Vec[i]=2.0*(drand48()-1.0);
	Ev=Vec.vec_norm();
	Iterno=0;
	do
	{
	    Oldvec=Vec; Oldev=Ev;
	    Vec=Matrix*Oldvec;
	    Ev=Oldvec*Vec;
	    Vec.vec_norm();
	}
	while (++Iterno<MAXITERNO && fabs(Ev-Oldev)>EPS*Oldev);
	
	if (Iterno>=MAXITERNO)	// no convergence :-(
	{ Evalno--; continue; }
	
	if (Ev>0)   // positive eigenvalue found
	{
	    Eval[Posevalno]=Ev;
	    Evec.col(Vec, Posevalno);
	    ++Posevalno;
	    if (Posevalno==Poseno)  // no more positives required
		return(Posevalno);  // don't do anything else
	}
	
	// deflation
	for (i=0; i<Size; i++)
	    for (j=0; j<=i; j++)
		Matrix[i][j]-=Ev*Vec[i]*Vec[j];
    }
    return(Posevalno);	// less than needed were found
}
// END of eigen_positer()

/* eigen_poscheb(): generates the first Poseno eigenvalues and
 * eigenvectors. Syntax is the same as that of eigen_positer()
 * but the algorithm used is the Chebyshev iteration (more robust).
 */
int eigen_poscheb(int Poseno, const Trimat_& Mat, Vector_& Eval, Sqmat_& Evec)
{
    /* This Chebyshev iteration is based on the algorithm given
     * in Crippen & Havel, Distance Geometry and Molecular Conformation, 
     * 1988., p.311.
     */
    
    const double EPS=1e-6;
    const int MAXITERNO=100;

    int Evalno, Iterno, i, j, Posevalno=0, Size=Mat.rno();
    if (Poseno>Size) Poseno=Size;
    
    Trimat_ Matrix(Mat);    // original preserved
    
    Vector_ Q0(Size), Q1(Size), Q2(Size), Mq1(Size);
    double Scale, Ev, Oldev, Realev, Len2;
    
    Eval.set_values(); Evec.set_values();   // zero output
    for (Evalno=0; Evalno<Size; Evalno++)
    {
	// scale the matrix
	Scale=Matrix.get_trace()/Size;
	Matrix/=Scale;
	
	// init as random unit vector
	for (i=0; i<Size; i++)
	    Q1[i]=2.0*(drand48()-1.0);
	Q1.vec_norm();
	
	// init the Chebyshev sequence
	Q2=Matrix*Q1; Mq1=Matrix*Q2;
	Ev=Q2.vec_len();
	Iterno=0;
	do
	{
	    Q0=Q1; Q1=Q2; Oldev=Ev;
	    Q2=2.0*Mq1-Q0;
	    Mq1=Matrix*Q2;
	    Len2=Q2.vec_len2();
	    Ev=(Mq1*Q2)/Len2;
	}
	while (++Iterno<MAXITERNO && fabs(Ev-Oldev)>EPS*Oldev);
	
	if (Iterno>=MAXITERNO)	// no convergence :-(
	{ Evalno--; continue; }
	
	Realev=Ev*Scale;    // get back "real" eigenvalue
	Q2.vec_norm();
	if (Realev>0.0)   // positive eigenvalue found
	{
	    Eval[Posevalno]=Realev;
	    Evec.col(Q2, Posevalno);
	    ++Posevalno;
	    if (Posevalno==Poseno)  // no more positives required
		return(Posevalno);  // don't do anything else
	}
	
	// deflation (of the scaled thing)
	for (i=0; i<Size; i++)
	    for (j=0; j<=i; j++)
		Matrix[i][j]-=Ev*Q2[i]*Q2[j];
	Matrix*=Scale;	// get back to unscaled
    }
    return(Posevalno);	// less than needed were found
}
// END of eigen_poscheb()

/* ---- EIGENROUTINES ---- */

/* NOTE: the following routines come from Numerical Recipes.
 * Originally, these were C functions (written in FORTRAN ;-> )
 * and the [1..N] indexing convention is still retained.
 * The conventional double array parameters are adjusted in
 * the calling routine to the [1..N] convention in the usual way, 
 * while the matrix array accesses are adjusted through the
 * ftn_idx() member function. See the corresponding comments
 * in eigen_ql() above. 
 * The routines now operate on square matrix objects and therefore
 * the "size" parameters are not passed explicitly. No rewriting
 * was done in the function bodies: these are valid C code still.
 * (Some editing and "register" local vars were added by A.A. in
 * the old C implementation.)
 */

/* tred2(): Housholder tridiagonalisation. a is a real, symmetric matrix
 * (size n*n). The main diagonal of the tridiag. output is returned in d,
 * the second diagonal in e with e[1]==0.0. On return, a contains the
 * transformation matrix "Q". 'd' and 'e' are assumed to have the
 * correct lengths adjusted before the call.
 */
static void tred2(Sqmat_& a, double *d, double *e)
{
    register int l,k,j,i, n=a.rno();	// get size
    register double scale,hh,h,g,f;

    for (i=n;i>=2;i--) {
	l=i-1;
	h=scale=0.0;
	if (l > 1) {
	    for (k=1;k<=l;k++)
		scale += fabs(a[i][k]);
	    if (scale < QL_EPSILON)
		e[i]=a[i][l];
	    else {
		for (k=1;k<=l;k++) {
		    a[i][k] /= scale;
		    h += a[i][k]*a[i][k];
		}
		f=a[i][l];
		g = (RND0(f)>0.0) ? -sqrt(h) : sqrt(h);
		e[i]=scale*g;
		h -= f*g;
		a[i][l]=f-g;
		f=0.0;
		for (j=1;j<=l;j++) {
		/* Next statement can be omitted if eigenvectors not wanted */
		    a[j][i]=a[i][j]/h;
		    g=0.0;
		    for (k=1;k<=j;k++)
			g += a[j][k]*a[i][k];
		    for (k=j+1;k<=l;k++)
			g += a[k][j]*a[i][k];
		    e[j]=g/h;
		    f += e[j]*a[i][j];
		}
		hh=f/(h+h);
		for (j=1;j<=l;j++) {
		    f=a[i][j];
		    e[j]=g=e[j]-hh*f;
		    for (k=1;k<=j;k++)
			a[j][k] -= (f*e[k]+g*a[i][k]);
		}
	    }
	} else
	    e[i]=a[i][l];
	d[i]=h;
    }
    /* Next statement can be omitted if eigenvectors not wanted */
    d[1]=e[1]=0.0;
    /* Contents of this loop can be omitted if eigenvectors not
	    wanted except for statement d[i]=a[i][i]; */
    for (i=1;i<=n;i++) {
	l=i-1;
	if (RND0(d[i])!=0.0) {	/* !=0.0 added */
	    for (j=1;j<=l;j++) {
		g=0.0;
		for (k=1;k<=l;k++)
		    g += a[i][k]*a[k][j];
		for (k=1;k<=l;k++)
		    a[k][j] -= g*a[k][i];
	    }
	}
	d[i]=(fabs(a[i][i])<QL_EPSILON)? 0.0: a[i][i];
	a[i][i]=1.0;
	for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
    }
}
// END of tred2()

/* tqli(): QL algorithm with implicit shifts on tridiagonal matrices.
 * The main diagonal is in d, the second diagonal is in e, with e[1]
 * ignored. If d and e were obtained from a general symmetric
 * real matrix by Housholder transformation by tred2(), then z should
 * contain the square transformation matrix "Q" on input; otherwise it should
 * be the unit matrix. On output, d contains the eigenvalues and z the
 * eigenvectors, with the k-th _column_ corresponding to the k-th eigenvalue.
 * Itno supplies the maximum allowable no. of iterations (an addition
 * by A.A.). Return value: 0 if OK, 1 if iteration limit has been
 * exceeded (added by A.A. to replace the nrerror() error message function
 * originally used in Numerical Recipes routines).
 */
static int tqli(double *d, double *e, Sqmat_& z, int Itno)
{
    register int m,l,iter,i,k, n=z.rno();
    register double s,r,ra,p,g,f,c,b;
    register float dd;

    for (i=2;i<=n;i++) e[i-1]=e[i];
    e[n]=0.0;
    for (l=1;l<=n;l++) {
	iter=0;
	do {
	    for (m=l;m<=n-1;m++) {
		dd=(float)(fabs(d[m])+fabs(d[m+1]));
		if ((float)fabs(e[m])+dd ==dd) break; 
	    }
	    if (m != l) {
		if (iter++ >= Itno)	/* too many iters */
		    return(1);
		g=(d[l+1]-d[l])/(2.0*e[l]);
		r=sqrt((g*g)+1.0);
		ra=(RND0(g)<0.0) ? -fabs(r) : fabs(r);
		g=d[m]-d[l]+e[l]/(g+ra);
		s=c=1.0;
		p=0.0;
		for (i=m-1;i>=l;i--) {
		    f=s*e[i];
		    b=c*e[i];
		    if (fabs(f) >= fabs(g)) {
			c=g/f;
			r=sqrt((c*c)+1.0);
			e[i+1]=f*r;
			c *= (s=1.0/r);
		    } else {
			s=f/g;
			r=sqrt((s*s)+1.0);
			e[i+1]=g*r;
			s *= (c=1.0/r);
		    }
		    g=d[i+1]-p;
		    r=(d[i]-g)*s+2.0*c*b;
		    p=s*r;
		    d[i+1]=g+p;
		    g=c*r-b;
		    /* Next loop can be omitted if eigenvectors not wanted */
		    for (k=1;k<=n;k++) {
			f=z[k][i+1];
			z[k][i+1]=s*z[k][i]+c*f;
			z[k][i]=c*z[k][i]-s*f;
		    }
		}
		d[l]=d[l]-p;
		e[l]=g;
		e[m]=0.0;
	    }
	} while (m != l);
    }
    return(0);
}
/* END of tqli */

#undef QL_EPSILON
#undef RND0

// ==== END OF FUNCTIONS Ql.c++ ====
