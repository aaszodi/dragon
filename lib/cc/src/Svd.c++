// ==== FUNCTIONS Svd.c++ ====

/* Singular value decomposition based on the algorithm
 * in Numerical Recipes.
 * The two matrices and the weight vector produced by SVD
 * are bundled into a little class.
 */

// SGI C++, IRIX 6.2, 20. June 1998. Andris

/* ----	HEADER ---- */

#include "Svd.h"

// ---- PROTOTYPES ----

extern "C"
{
    static double max(double a, double b)
    {
	return((a>b)? a: b);
    }
    
    static double sign(double a, double b)
    {
	return((b>=0.0)? fabs(a): -fabs(a));
    }
}

// ==== Static initialisation ====

const Safety_ Svd_::Safe; 	// see Safety.h for defaults

// ==== Svd_ MEMBER FUNCTIONS ====

// ---- AUXILIARIES ----

/* utb(): the vector B is premultiplied by the transpose of the
 * U matrix in the decomposition and the result is placed into
 * Utb. Dim checks are not carried out. Private
 */
inline
void Svd_::utb(const Vector_& B, Vector_& Utb) const
{
    register unsigned int i, j;
    register double Temp;
    
    for (j=0; j<C; j++)
    {
	Temp=0.0;
	for (i=0; i<R; i++)
	    Temp+=U[i][j]*B[i];
	Utb[j]=Temp;
    }
}
// END of utb()

// ---- Constructors ----

/* Set up SVD for a Row x Col matrix where Row>=Col. If Row<Col,
 * then the rows of U will be padded to make it ColxCol
 * with a warning. If either Row or Col is 0, will be changed
 * to 3 (cf matrix/vector ctors) and a warning printed.
 * The default size is also 3x3.
 * The members will be initialised to 0.0 (all elements).
 */
Svd_::Svd_(unsigned int Row, unsigned int Col) :
    U(((Row<Col)? Col: Row), Col), // row padding
    W(Col), 
    V(Col)
{
    // size check warnings (actually the members are init'd already)
    if (!Row) { cerr<<"\n? Svd_: Row==0, was set to 3\n"; Row=3; }
    if (!Col) { cerr<<"\n? Svd_: Col==0, was set to 3\n"; Col=3; }
    C=Col; Rorig=Row;
    if (Row<Col)
    {
	cerr<<"\n? Svd_: Row="<<Row<<", padded to "<<Col<<endl;
	Row=Col;
    }
    R=Row;
}

// ---- Size ----

/* set_size(): modifies the sizes of the SVD components to
 * accommodate a Row x Col matrix. If Row < Col, then the
 * rows of U will be padded to give a Col x Col matrix.
 * If Row or Col is 0, then 3 will be used instead and a warning printed.
 * Without arguments, the size is set to 3x3.
 */
void Svd_::set_size(unsigned int Row, unsigned int Col)
{
    if (R==Row && C==Col) return;	// no change
    
    if (!Row) { cerr<<"\n? Svd_::set_size(): Row==0, was set to 3\n"; Row=3; }
    if (!Col) { cerr<<"\n? Svd_::set_size(): Col==0, was set to 3\n"; Col=3; }

    C=Col; Rorig=Row;
    if (Row<Col)
    {
	cerr<<"\n? Svd_::set_size(): Row="<<Row<<", padded to "<<Col<<endl;
	Row=Col;
    }
    R=Row;
    
    U.set_size(R, C); W.dim(C); V.set_size(C);
}
// END of set_size()

// ---- SVD routines ----

/* make_decomp(): the SVD according to Num. Recipes, done
 * on the general matrix A (which is preserved). The
 * three member objects U, W, V will be set in the calling
 * object. Return value: 1 if the iteration limit (hard-coded
 * in svd_core()) is exceeded, 0 if OK, -1 if a dimension
 * mismatch occurred.
 */
int Svd_::make_decomp(const Matrix_& A)
{
    // adjust size if necessary
    set_size(A.rno(), A.cno());
    
    // copy A to U for decomposition
    if (A.rno()>=A.cno()) U=A;
    else    // copy rows only
    {
	for (unsigned int i=0; i<A.rno(); i++)
	    U.row(A.row(i), i);
    }
    
    // actual decomposition: returns max. iterno if exceeded, 0 if OK
    int Err=svd_core();
    if (Err)
    {
	cerr<<"\n? Svd_::make_decomp(): No convergence in "<<Err<<" iteration(s)\n";
    }
    return(Err);
}
// END of make_decomp()

/* rank_cond(): checks the N singular values W[] of a matrix 
 * after SVD. The condition number in *Cond
 * (ratio of the largest and Safe.small()est singular value) is also
 * calculated if Cond!=NULL. The singular values which are Safe.small()er than
 * Eps times the largest are set to 0.0.
 * Return value: the rank of the matrix.
 */
unsigned int Svd_::rank_cond(double Eps, double *Cond)
{
    register double Wmax=-HUGE_VAL, Wmin=HUGE_VAL;
    register unsigned int i;
    
    // get the largest and Safe.small()est singular value
    for (i=0; i<C; i++)
    {
	if (W[i]>Wmax) Wmax=W[i];
	if (W[i]<Wmin) Wmin=W[i];
    }
    
    // calc the condition number: set to HUGE_VAL if Wmin==0.0
    if (Cond!=NULL)
	*Cond=(Wmin==0.0)? HUGE_VAL: Wmax/Wmin;
    
    /* set all singular values which are Safe.small()er than Eps*Wmax
     * to zero: this is the conditioning. Calc the rank
     */
    Wmax*=fabs(Eps);
    unsigned int Rank, Maxrank;
    Rank=Maxrank=(Rorig>=C)? C: Rorig;
    for (i=0; i<Maxrank; i++)
	if (W[i]<Wmax) 
	{
	    W[i]=0.0; Rank--;
	}
    return(Rank);
}
// END of rank_cond()

/* lin_solve(): back-substitution routine for solving linear equations
 * AX=B. A should be SV-decomposed into U, W and V' by make_decomp()
 * and the weight vector should be "conditioned" (Safe.small() entries
 * zeroed) by rank_cond() prior to the call to this routine.
 * Return value: the X vector.
 */
Vector_ Svd_::lin_solve(const Vector_& B) const
{
    // dim mismatches
    unsigned int Bdim=B.dim();
    if (R==Rorig && Bdim!=R || R>Rorig && !(Bdim==R || Bdim==Rorig))
    {
	cerr<<"\n? Svd_::lin_solve: B has wrong dimensions, nullvector returned\n";
	cerr<<"R="<<R<<", Rorig="<<Rorig<<", Bdim="<<Bdim<<endl;
	Vector_ X(C);	// a null-vector
	return(X);
    }
    
    Vector_ Wub(C);	// this will be W*U'*B
    // the "fewer eqns than unknowns" case
    if (R>Rorig && Bdim==Rorig)
    {
	cerr<<"\n? Svd_::lin_solve: fewer eqns ("<<Rorig
	    <<") than unknowns ("<<R<<"), B zero-padded\n";
	Vector_ Bpad(B); Bpad.dim(R);   // pad B w/ zeroes
	utb(Bpad, Wub);
    }
    else	// normal cases
	utb(B, Wub);
    
    for (unsigned int j=0; j<C; j++)
	if (W[j]==0.0) Wub[j]=0.0; else Wub[j]/=W[j];
    
    // get solution vector and return
    Vector_ X=V*Wub;
    return(X);
}
// END of lin_solve()

// ---- SVD CORE ----

/* svd_core(): Performs the SVD of the member matrix U (not preserved).
 * (It is assumed to have been initialised before the call.)
 * Based on the ANSI C code of the svdcmp() routine in Numerical Recipes, 
 * 2nd edition (1992). Re-written to C++ in dbl prec so that it operates on
 * the Svd_ member matrices. Indexing is not changed: the FORTRAN
 * indexing convention is used within and the matrix ptrs are
 * restored on exit. Returns error status: 0 if OK, >0 on error,
 * SVD_ITMAX in particular if iteration limit (hardcoded) exceeded. Private
 */
int Svd_::svd_core()
{
	const int SVD_ITMAX=30, m=R, n=C;	// max. itno, row and column size
	register int flag,i,its=1,j,jj,k,l,nm;
	int Retval=0;
	register double anorm=0.0,c,f,g=0.0,h,s,scale=0.0,recscale, x,y,z;
	register double *rv1, *Warr;
	
	/* FORTRAN indexing: the matrices are converted to [1..m][1..n]
	 * here. There is no corresponding member function in my Vector_
	 * class so temp vectors are allocated as C arrays and
	 * the idx shift is done on them. W gets its value at the end
	 * at the expense of a copy.
	 */
	U.ftn_idx(); V.ftn_idx();
	rv1=new double [n]; rv1--;
	Warr=new double [n]; Warr--;
	
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(U[k][i]);
			if (scale>Safe.small())
			{
				recscale=1.0/scale;
				for (k=i;k<=m;k++) {
					U[k][i] *= recscale;
					s += U[k][i]*U[k][i];
				}
				f=U[i][i];
				g = -sign(sqrt(s),f);
				h=Safe.safe_div(1.0,f*g-s,__LINE__);
				U[i][i]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += U[k][i]*U[k][j];
					f=s*h;
					for (k=i;k<=m;k++) U[k][j] += f*U[k][i];
				}
				for (k=i;k<=m;k++) U[k][i] *= scale;
			}
		}
		Warr[i]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(U[i][k]);
			if (scale>Safe.small())
			{
				recscale=1.0/scale;
				for (k=l;k<=n;k++) {
					U[i][k] *= recscale;
					s += U[i][k]*U[i][k];
				}
				f=U[i][l];
				g = -sign(sqrt(s),f);
				h=Safe.safe_div(1.0,f*g-s,__LINE__);
				U[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=U[i][k]*h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += U[j][k]*U[i][k];
					for (k=l;k<=n;k++) U[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++) U[i][k] *= scale;
			}
		}
		anorm=max(anorm,(fabs(Warr[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (fabs(g)>Safe.small()) {
				g=1.0/g;
				for (j=l;j<=n;j++)
					V[j][i]=(U[i][j]/U[i][l])*g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += U[i][k]*V[k][j];
					for (k=l;k<=n;k++) V[k][j] += s*V[k][i];
				}
			}
			for (j=l;j<=n;j++) V[i][j]=V[j][i]=0.0;
		}
		V[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=(m<n)? m: n;i>=1;i--) {
		l=i+1;
		g=Warr[i];
		for (j=l;j<=n;j++) U[i][j]=0.0;
		if (fabs(g)>Safe.small()) {
			g=1.0/g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += U[k][i]*U[k][j];
				f=Safe.safe_div(s,U[i][i],__LINE__);
				f*=g;
				for (k=i;k<=m;k++) U[k][j] += f*U[k][i];
			}
			for (j=i;j<=m;j++) U[j][i] *= g;
		} else for (j=i;j<=m;j++) U[j][i]=0.0;
		(U[i][i])++;
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=SVD_ITMAX;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if ((fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((fabs(Warr[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((fabs(f)+anorm) == anorm) break;
					g=Warr[i];
					h=Safe.pythag(f,g);
					Warr[i]=h;
					h=Safe.safe_div(1.0,h,__LINE__);
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=U[j][nm];
						z=U[j][i];
						U[j][nm]=y*c+z*s;
						U[j][i]=z*c-y*s;
					}
				}
			}
			z=Warr[k];
			if (l == k) {
				if (z < 0.0) {
					Warr[k] = -z;
					for (j=1;j<=n;j++) V[j][k] = -V[j][k];
				}
				break;
			}
			
			// trouble, leave...
			if (its == SVD_ITMAX) 
			{
			    cerr<<"\n? Svd_::svd_core(): No convergence\n";
			    Retval=SVD_ITMAX;	// non-zero
			    goto ERROR;
			}
			    
			x=Warr[l];
			nm=k-1;
			y=Warr[nm];
			g=rv1[nm];
			h=rv1[k];
			f=Safe.safe_div((y-z)*(y+z)+(g-h)*(g+h),(2.0*h*y),__LINE__);
			g=Safe.pythag(f,1.0);
		    	if (fabs(x)<Safe.small())
			{
		    	    cerr<<"\n? Svd_::svd_core(): Impending div by zero (x="<<x<<")\n";
			    Retval=its;
			    goto ERROR;
			}
			f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=Warr[i];
				h=s*g;
				g=c*g;
				z=Safe.pythag(f,h);
		    	    	if (z<Safe.small())
			    	{
		    	    	    cerr<<"\n? Svd_::svd_core(): Impending div by zero (z="<<z<<")\n";
			    	    Retval=its;
			    	    goto ERROR;
			    	}
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=V[jj][j];
					z=V[jj][i];
					V[jj][j]=x*c+z*s;
					V[jj][i]=z*c-x*s;
				}
				z=Safe.pythag(f,h);
				Warr[j]=z;
				if (z>Safe.small()) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=U[jj][j];
					z=U[jj][i];
					U[jj][j]=y*c+z*s;
					U[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			Warr[k]=x;
		}
	}
	
	// get back to C indexing and clean up
	// on error, non-0 is returned
	ERROR:
	U.c_idx(); V.c_idx();
	Warr++; rv1++;
	delete [] rv1;
	if (!Retval) W=Vector_(Warr, n);    // was OK, copy weight vector
	delete [] Warr;
	return(Retval);
}
#undef PYTHAG
// END of svd_core()

// ---- GLOBAL FUNCTIONS ----

/* <<: the overloaded output operator. Prints a neat listing to Out */
ostream& operator<<(ostream& Out, const Svd_& Svd)
{
    Out<<Svd.Rorig<<"x"<<Svd.C<<" singular decomposition";
    unsigned int Ex=Svd.R-Svd.Rorig;
    if (Ex) Out<<" ("<<Ex<<" row"<<((Ex==1)? " ":"s ")<<"added)\n";
    else Out<<'\n';
    Out<<"Singular values:\n"<<Svd.W;
    Out<<"The U matrix:\n"<<Svd.U;
    Out<<"The V matrix:\n"<<Svd.V;
    return(Out);
}
// END of <<

// ==== END OF FUNCTIONS Svd.c++ ====
