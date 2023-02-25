// ==== PROJCT DRAGON: METHODS Rsmdiag.c++ ====

/* Real symmetric matrix diagonalisation class.
 * Use if only selected eigenvectors are wanted.
 * Based on NETLIB routines.
 */

// SGI C++, IRIX 6.2, 21. June 1998. Andris Aszodi

// ---- MODULES ----

#include "Rsmdiag.h"

// ---- Static initialisation ----

const Safety_ Rsmdiag_::Safe;	// see Safety.h for default settings

// ==== METHODS ====

// ---- Diagonalisation ----

/* get_evals(): obtain all eigenvalues of a symmetric matrix Mat
 * and put them into the Evals vector (size set if necessary).
 * The eigenvalues are sorted in decreasing order.
 * Return value: 0 if OK, k>0 if the k:th eigenvalue failed to
 * converge.
 */
int Rsmdiag_::get_evals(const Trimat_& Mat, Vector_& Evals)
{
    int Size=Mat.rno();
    
    // init
    int Err=0;
    
    c_idx();
    Qmat=Mat;	// working copy
    set_size(Size); // adjust array lengths
    ftn_idx();	// set to F77 indexing
    
    // transform Qmat to tridiagonal form
    tred_1();
    
    // obtain all eigenvalues of the tridiagonal matrix (QL transform)
    Err=imt_qlv();
    if (Err)
    {
	cerr<<"\n? Rsmdiag_::get_evals(): "<<Err<<". eigenvalue not found\n";
	// ####
	cout<<"The matrix was:\n"<<Mat;
	exit(EXIT_FAILURE);
	// ####
	return(Err);
    }
    
    // copy eigenvalues from [1..Size] temp array
    Evals.dim(Size);
    for (register unsigned i=0; i<Size; ++i)
	Evals[i]=w[i+1]; 
    return(0);
}
// END of get_evals()

/* get_evecs(): obtain the first Evno eigenvectors. It is assumed but
 * not checked that a corresponding get_evals() has been executed
 * beforehand. The eigenvectors will be placed into the first Evno
 * columns of the Evecs matrix whose size will be set silently.
 * Return value: 0 if OK, r>0 if the r:th eigenvector failed
 * to converge.
 */
int Rsmdiag_::get_evecs(int Evno, Sqmat_& Evecs)
{
    Evecs.set_size(Qmat.rno());
    Evecs.ftn_idx();
    
    // generate the first Evno eigenvectors of the tridiagonal matrix
    int Err=inv_iter(Evno, Evecs);
    if (Err)
    {
	cerr<<"\n? Rsmdiag_::get_evecs(): "<<(-Err)<<". eigenvector not found\n";
	// #####
	cout<<"Evno="<<Evno<<", Evecs was:"<<Evecs;
	exit(EXIT_FAILURE);
	// #####
    }
    else
	tr_bak1(Evno, Evecs);	// transform tridiagonal eigenvectors into original

    Evecs.c_idx();
    return(-Err);
}
// END of get_evecs()

// ---- Diag from NETLIB ----

/* NOTE: FORTRAN [1..N] style indexing is used throughout. Perform
 * base ptr shifts before invoking the routines.
 */

/* tred_1(): transforms Qmat into tridiagonal form.
 * (d is the main diagonal, e the subdiagonal, e2 the squared subdiagonal).
 * The transformation information is written into Qmat on return.
 */
void Rsmdiag_::tred_1()
{
    register double f, g, h, scale, d__1;
    register int n=Qmat.rno(), i, j, k, l, ii, jp1;

    for (i = 1; i <= n; ++i)
    {
	d[i] = Qmat[n][i];
	Qmat[n][i] = Qmat[i][i];
    }
    for (ii = 1; ii <= n; ++ii)
    {
	i = n + 1 - ii;
	l = i - 1;
	h = 0.;
	scale = 0.;
	if (l < 1) {
	    goto L130;
	}
	/* scale row */
	for (k = 1; k <= l; ++k)
	    scale += fabs(d[k]);

	if (scale > Safe.small()) {
	    goto L140;
	}

	for (j = 1; j <= l; ++j)
	{
	    d[j] = Qmat[l][j];
	    Qmat[l][j] = Qmat[i][j];
	    Qmat[i][j] = 0.;
	}

L130:
	e[i] = e2[i] = 0.;
	continue;

L140:
	d__1=1.0/scale;
	for (k = 1; k <= l; ++k)
	{
	    d[k] *= d__1;  /* was /= */
	    h += d[k] * d[k];
	}

	e2[i] = scale * scale * h;
	f = d[l];
	d__1 = sqrt(h);
	g = (f>=0.0)? -fabs(d__1): fabs(d__1); /* was g = -d_sign(&d__1, &f); */
	e[i] = scale * g;
	h -= f * g;
	d[l] = f - g;
	if (l == 1) {
	    goto L285;
	}
/*     .......... form a*u .......... */
	for (j = 1; j <= l; ++j) e[j] = 0.;

	for (j = 1; j <= l; ++j)
	{
	    f = d[j];
	    g = e[j] + Qmat[j][j] * f;
	    jp1 = j + 1;
	    if (l >= jp1)
	    {
		for (k = jp1; k <= l; ++k)
		{
		    d__1=Qmat[k][j];
		    g += d__1 * d[k];	
		    e[k] += d__1 * f;
		}
	    }
	    e[j] = g;
	}
/*     .......... form p .......... */
	f = 0.;

	d__1=Safe.safe_div(1.0,h,__LINE__);
	for (j = 1; j <= l; ++j)
	{
	    e[j] *= d__1;
	    f += e[j] * d[j];
	}

	h = Safe.safe_div(f,2*h,__LINE__);
/*     .......... form q .......... */
	for (j = 1; j <= l; ++j)
	    e[j] -= h * d[j];

/*     .......... form reduced a .......... */
	for (j = 1; j <= l; ++j)
	{
	    f = d[j];
	    g = e[j];

	    for (k = j; k <= l; ++k)
		Qmat[k][j] -= f * e[k] + g * d[k];
	}

L285:
	for (j = 1; j <= l; ++j)
	{
	    f = d[j];
	    d[j] = Qmat[l][j];
	    Qmat[l][j] = Qmat[i][j];
	    Qmat[i][j] = f * scale;
	}
    }	    /* for ii */

} /* tred1_ */

/* imt_qlv(): finds the eigenvalues of an n x n real symmetric tridiagonal
 * matrix (main diagonal in d, off-diagonal in e, squares of off-diag
 * values in e2, as supplied by tred_1()). The eigenvalues are returned
 * in w, sorted in descending order.
 * Return value: 0 if OK, j>0 if the j:th eigenvalue was not found
 * within MAXITERNO iterations.
 */
int Rsmdiag_::imt_qlv()
{
    const int MAXITERNO=30;
    
    register double b, c, f, g, p, r, s, tst1, tst2;
    register int n=Qmat.rno(), i, j, k, l, m, ii, tag, mml;
    int ierr=0;
    
    double *rv1=(double *) calloc(n, sizeof(double));	/* aux array */
    --rv1;  /* adjust ftn idx */

    k = tag = 0;

    for (i = 1; i <= n; ++i)
    {
	w[i] = d[i];
	if (i != 1) {
	    rv1[i - 1] = e[i];
	}
    }

    e2[1] = rv1[n] = 0.;

    for (l = 1; l <= n; ++l) {
	j = 0;
/*     .......... look for Safe.small() sub-diagonal element .......... */
L105:
	for (m = l; m <= n; ++m) {
	    if (m == n) {
		goto L120;
	    }
	    tst1=fabs(w[m])+fabs(w[m+1]);
	    tst2 = tst1 + fabs(rv1[m]);
	    if (tst2 == tst1) {
		goto L120;
	    }
/*     .......... guard against underflowed element of e2 ........
.. */
	    if (fabs(e2[m + 1]) < Safe.small()) {
		goto L125;
	    }
/* L110: */
	}

L120:
	if (m <= k) {
	    goto L130;
	}
	if (m != n) {
	    e2[m + 1] = 0.;
	}
L125:
	k = m;
	++tag;
L130:
	p = w[l];
	if (m == l) {
	    goto L215;
	}
	if (j == MAXITERNO) {
	    ierr=l; goto L1001;	/* no convergence */
	}
	++j;
/*     .......... form shift .......... */
	g = (w[l + 1] - p) / (rv1[l] * 2.);
	r = Safe.pythag(g, 1.0);
	g = w[m] - p + rv1[l] / (g + (g>=0.0)? fabs(r): -fabs(r) /*d_sign(&r, &g)*/ );
	s = c = 1.;
	p = 0.;
	mml = m - l;
/*     .......... for i=m-1 step -1 until l do -- .......... */
	for (ii = 1; ii <= mml; ++ii) {
	    i = m - ii;
	    f = s * rv1[i];
	    b = c * rv1[i];
	    r = Safe.pythag(f, g);
	    rv1[i + 1] = r;
	    if (fabs(r) < Safe.small()) {
		goto L210;
	    }
	    s = f / r;
	    c = g / r;
	    g = w[i + 1] - p;
	    r = (w[i] - g) * s + c * 2. * b;
	    p = s * r;
	    w[i + 1] = g + p;
	    g = c * r - b;
/* L200: */
	}

	w[l] -= p;
	rv1[l] = g;
	rv1[m] = 0.;
	goto L105;
/*     .......... recover from underflow .......... */
L210:
	w[i + 1] -= p;
	rv1[m] = 0.;
	goto L105;
/*     .......... order eigenvalues .......... */
L215:
	/* order eigenvalues: into _descending_ order
	 * as opposed to the ascending order in the
	 * original NETLIB implementation
	 */
	for (i=1; i<l; ++i)
	{
	    if (p>=w[i])    /* p will be the i:th */
	    {
		for (ii=l; ii>i; --ii)	/* shift away rest */
		{
		    w[ii]=w[ii-1]; Index[ii]=Index[ii-1];
		}
		break;
	    }
	}
	w[i]=p; Index[i]=tag;

/* L290: */
    }

L1001:
    ++rv1; free(rv1);
    e2[1]=2.0;	/* descending order */
    return(ierr);
}
/* END of imt_qlv() */

/* inv_iter(): generates the eigenvectors belonging to the first m
 * eigenvalues (supplied in w) of the n x n real symmetric tridiagonal
 * matrix (diagonal in d, off-diagonal in e, squared off-diagonal in
 * e2, e2[1]==0.0 if the eigenvalues are in ascending order, 2.0 if
 * in descending order: imt_qlv() takes care of this). ind holds
 * the index array for the submatrices belonging to the eigenvalues.
 * On return, z contains the first m eigenvectors.
 * Return value: 0 if OK, -r if the r:th eigenvector failed to converge.
 */
int Rsmdiag_::inv_iter(int m, Sqmat_& z) const
{
    /* f2c-generated variables */
    double d__3, d__4;
    
    double *rv1, *rv2, *rv3, *rv4, *rv6;    /* aux arrays */
    register double norm, u, v, order, x0, x1, uk, xu, eps2, eps3, eps4;
    register int n=z.rno(), i, j, p, q, r, s, group, ii, jj, ip, tag, its, ierr=0;

    /* This subroutine is a translation of the inverse iteration technique */
    /* in the algol procedure tristurm by Peters and Wilkinson. */
    /* Handbook for Auto. Comp., Vol.ii-Linear Algebra, 418-439(1971). */

    /* check for silly values */
    if (!m) return(-1);
    m=abs(m);
    
    /* Alloc aux arrays as one big array */
    rv1=(double *) calloc(5*n, sizeof(double)); --rv1;
    rv2=rv1+n; rv3=rv2+n; rv4=rv3+n; rv6=rv4+n;

    tag = 0;
    order = 1. - e2[1];
    q = 0;
/*     .......... establish and process next submatrix .......... */
L100:
    p = q + 1;

    for (q = p; q <= n; ++q) {
	if (q == n) {
	    goto L140;
	}
	if (fabs(e2[q + 1])<Safe.small()) {
	    goto L140;
	}
/* L120: */
    }
/*     .......... find vectors by inverse iteration .......... */
L140:
    ++tag;
    s = 0;

    for (r = 1; r <= m; ++r) {
	if (Index[r] != tag) {
	    goto L920;
	}
	its = 1;
	x1 = w[r];
	if (s != 0) {
	    goto L510;
	}
/*     .......... check for isolated root .......... */
	xu = 1.;
	if (p != q) {
	    goto L490;
	}
	rv6[p] = 1.;
	goto L870;
L490:
	norm = fabs(d[p]);
	ip = p + 1;

	for (i = ip; i <= q; ++i) {
/* L500: */
/* Computing MAX */
	    d__3 = norm; d__4 = fabs(d[i]) + fabs(e[i]);
	    norm = (d__3>=d__4)? d__3: d__4;	/* max(d__3,d__4); */
	}
/*     .......... eps2 is the criterion for grouping, */
/*                eps3 replaces zero pivots and equal */
/*                roots are modified by eps3, */
/*                eps4 is taken very Safe.small() to avoid overflow .........
. */
	eps2 = norm * .001;
	eps3 = DBL_EPSILON*norm;    // simpler than original: 9-V-98
	uk = (double) (q - p + 1);
	eps4 = uk * eps3;
	uk = eps4 / sqrt(uk);
	s = p;
L505:
	group = 0;
	goto L520;
/*     .......... look for close or coincident roots .......... */
L510:
	if (fabs(x1 - x0) >= eps2) goto L505;

	++group;
	if (order * (x1 - x0) <= 0.) {
	    x1 = x0 + order * eps3;
	}
/*     .......... elimination with interchanges and */
/*                initialization of vector .......... */
L520:
	v = 0.;

	for (i = p; i <= q; ++i) {
	    rv6[i] = uk;
	    if (i == p) {
		goto L560;
	    }
	    if (fabs(e[i]) < fabs(u)) goto L540;

/*     .......... warning -- a divide check may occur here if */
/*                e2 array has not been specified correctly ......
.... */
	    xu = Safe.safe_div(u,e[i],__LINE__);
	    rv4[i] = xu;
	    rv1[i - 1] = e[i];
	    rv2[i - 1] = d[i] - x1;
	    rv3[i - 1] = 0.;
	    if (i != q) {
		rv3[i - 1] = e[i + 1];
	    }
	    u = v - xu * rv2[i - 1];
	    v = -xu * rv3[i - 1];
	    goto L580;
L540:
	    xu = Safe.safe_div(e[i],u,__LINE__);
	    rv4[i] = xu;
	    rv1[i - 1] = u;
	    rv2[i - 1] = v;
	    rv3[i - 1] = 0.;
L560:
	    u = d[i] - x1 - xu * v;
	    if (i != q) {
		v = e[i + 1];
	    }
L580:
	    ;
	}

	if (fabs(u) < Safe.small()) {
	    u = eps3;
	}
	rv1[q] = u;
	rv2[q] = 0.;
	rv3[q] = 0.;
/*     .......... back substitution */
/*                for i=q step -1 until p do -- .......... */
L600:
	for (ii = p; ii <= q; ++ii) {
	    i = p + q - ii;
	    rv6[i] = Safe.safe_div(rv6[i] - u * rv2[i] - v * rv3[i],rv1[i],__LINE__);
	    v = u;
	    u = rv6[i];
/* L620: */
	}
/*     .......... orthogonalize with respect to previous */
/*                members of group .......... */
	if (group == 0) {
	    goto L700;
	}
	j = r;

	for (jj = 1; jj <= group; ++jj) {
L630:
	    --j;
	    if (Index[j] != tag) {
		goto L630;
	    }
	    xu = 0.;

	    for (i = p; i <= q; ++i) 
		xu += rv6[i] * z[i][j];
	    for (i = p; i <= q; ++i) 
		rv6[i] -= xu * z[i][j];

/* L680: */
	}

L700:
	norm = 0.;
	for (i = p; i <= q; ++i) {
/* L720: */
	    norm += fabs(rv6[i]);
	}

	if (norm >= 1.) {
	    goto L840;
	}
/*     .......... forward substitution .......... */
	if (its == 5) {
	    goto L830;
	}
	if (norm>Safe.small()) {
	    goto L740;
	}
	rv6[s] = eps4;
	++s;
	if (s > q) {
	    s = p;
	}
	goto L780;
L740:
	xu = eps4 / norm;

	for (i = p; i <= q; ++i) rv6[i] *= xu;

/*     .......... elimination operations on next vector */
/*                iterate .......... */
L780:
	for (i = ip; i <= q; ++i) {
	    u = rv6[i];
/*     .......... if rv1(i-1) .eq. e(i), a row interchange */
/*                was performed earlier in the */
/*                triangularization process .......... */
	    if (rv1[i - 1] != e[i]) {
		goto L800;
	    }
	    u = rv6[i - 1];
	    rv6[i - 1] = rv6[i];
L800:
	    rv6[i] = u - rv4[i] * rv6[i - 1];
/* L820: */
	}

	++its;
	goto L600;
/*     .......... set error -- non-converged eigenvector .......... */
L830:
	ierr = -r;
	xu = 0.;
	goto L870;
/*     .......... normalize so that sum of squares is */
/*                1 and expand to full order .......... */
L840:
	u = 0.;
	for (i = p; i <= q; ++i) {
/* L860: */
	    u = Safe.pythag(u, rv6[i]);
	}

	xu = Safe.safe_div(1.0, u, __LINE__);

L870:
	for (i = 1; i <= n; ++i) z[i][r] = 0.;
	for (i = p; i <= q; ++i) z[i][r] = rv6[i] * xu;
	x0 = x1;
L920:
	;
    }

    if (q < n) {
	goto L100;
    }
    ++rv1; free(rv1);
    return(ierr);	/* OK */
}
/* END of inv_iter() */

/* tr_bak1(): forms the first m eigenvectors of a real symmetric matrix
 * by back-transforming the eigenvectors of the tridiagonal matrix
 * constructed by tred_1(). n is the matrix size, "a" contains
 * the lower triangle holding the orthogonal transforms as produced
 * by tred_1(), "e" holds the subdiagonal of the tridiagonal matrix
 * made by tred_1(), m is the number of eigenvectors desired, 
 * z on input holds the eigenvectors of the tridiagonal matrix
 * as constructed by inv_iter(). On output, the transformed 
 * eigenvectors will be in z's first m columns.
 */
void Rsmdiag_::tr_bak1(int m, Sqmat_& z) const
{
    /* Local variables */
    register int i, j, k, l, n=Qmat.rno();
    register double s;

    /*     This subroutine is a translation of the algol procedure trbak1, */
    /*     Num. Math. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson. */
    /*     Handbook for Auto. Comp., Vol.ii-Linear Algebra, 212-226(1971). */

    /* Function Body */
    m=abs(m); n=abs(n);
    if (!m || n<=1) return; /* rubbish */

    for (i = 2; i <= n; ++i)
    {
	l = i - 1;
	if (fabs(e[i])<Safe.small()) break;

	for (j = 1; j <= m; ++j)
	{
	    s = 0.0;

	    for (k = 1; k <= l; ++k)
		s += Qmat[i][k] * z[k][j];
	    
	    /* divisor below is negative of h formed in tred1. 
	     * double division avoids possible underflow
	     */
	    /* s = s / Qmat[i][l] / e[i]; */
	    s=Safe.safe_div(s,Qmat[i][l],__LINE__);
	    s=Safe.safe_div(s,e[i],__LINE__);
	    for (k = 1; k <= l; ++k)
		z[k][j] += s * Qmat[i][k];
	}
    }
}
/* END of tr_bak1() */

/* ---- Auxiliaries ---- */

/* c_idx(), ftn_idx(): adjust all ptrs to C or FORTRAN-style
 * indexing. This is necessary because the diagonalisation routines
 * have been translated from f77 code. The member Ftnidx stores
 * the actual indexing state. Cf. the c_idx(), ftn_idx() methods
 * in various matrix classes.
 */
void Rsmdiag_::c_idx()    // sets ptrs to C-style [0..N-1] indexing
{
    if (Ftnidx)
    {
	Qmat.c_idx();
	if (d!=NULL) { ++d; ++e; ++e2; ++w; }
	if (Index!=NULL) ++Index;
	Ftnidx=false;
    }
}
// END of c_idx()

void Rsmdiag_::ftn_idx()    // sets ptrs to FORTRAN-style [1..N] indexing
{
    if (!Ftnidx)
    {
	Qmat.ftn_idx();
	if (d!=NULL) { --d; --e; --e2; --w; }
	if (Index!=NULL) --Index;
	Ftnidx=true;
    }
}
// END of ftn_idx()   

/* set_size(): adjusts the length of the temp arrays to Size.
 * Size==0 is legal; all arrays are deleted if this is specified.
 * Note that the double arrays are in fact allocated as one
 * long array and the ptrs are adjusted accordingly.
 * Indexing is temporarily reset to C style if necessary.
 * Return value: the new size.
 */
unsigned int Rsmdiag_::set_size(unsigned int Size)
{
    bool Oldftnidx=Ftnidx;
    
    c_idx();	// go to C indexing
    delete [] d; delete [] Index;   // remove old
    if (Size)
    {
	d=new double [4*Size];	// one long array
	e=d+Size; e2=e+Size; w=e2+Size;	// set other ptrs
	Index=new int [Size];
	if (Oldftnidx) ftn_idx();	// go back to f77 indexing if needed
    }
    else    // reset: also used for dtor call
    {
	d=e=e2=w=NULL; Index=NULL;
    }
    return(Size);
}
// END of set_size()

// ==== END OF METHODS Rsmdiag.c++ ====
