// ==== FUNCTIONS Siva.c++ ====

/* Singular value decomposition based on the algorithm
 * in Pal Rozsa's book (Linearis algebra es alkalmazasai).
 * The two matrices and the weight vector produced by SVD
 * are bundled into a little class but access to them is
 * public.
 */

// SGI C++ 3.2.1, IRIX 5.2, 1. Dec. 1994. Andris

/* ----	HEADER ---- */

#include "Siva.h"

// ==== Siva_ MEMBER FUNCTIONS ====

// ---- Constructors ----

/* Set up SVD for a Row x Col matrix where Row>=Col. If Row<Col,
 * then the rows of U will be padded to make it ColxCol and a
 * warning is printed. If either Row or Col is 0, will be changed
 * to 1 (cf matrix/vector ctors) and another warning printed. 
 * The members will be initialised to 0.0 (all elements).
 */
Siva_::Siva_(unsigned int Row, unsigned int Col) :
    U(((Row<Col)? Col: Row), Col), // row padding
    W(Col), 
    V(Col)
{
    // size check warnings (actually the members are init'd already)
    if (!Row) { cerr<<"? Siva_: Row==0, was set to 1\n"; Row=1; }
    if (!Col) { cerr<<"? Siva_: Col==0, was set to 1\n"; Col=1; }
    C=Col; Rorig=Row;
    if (Row<Col)
    {
	cerr<<"? Siva_: ("<<Row<<"x"<<Col<<") size specified, set to ("
	    <<Col<<"x"<<Col<<")\n";
	Row=Col;
    }
    R=Row;
}

// ---- SVD routines ----

/* make_decomp(): the SVD according to P. Rozsa, done
 * on the general matrix A (which is preserved). The
 * three member objects U, W, V will be set in the calling
 * object. Return value: 1 if the iteration limit (hard-coded
 * in eigen_ql()) is exceeded, 0 if OK, -1 if a dimension
 * mismatch occurred.
 */
int Siva_::make_decomp(const Matrix_& A)
{
    // dim mismatch check
    if (Rorig!=A.rno() || C!=A.cno())
    {
	cerr<<"? Siva_::make_decomp(): dimension mismatch\n";
	return(-1);
    }
    
    // generate the A'A symmetric matrix
    Trimat_ Ata=trans_mprod(A);
    
    /* diagonalise A'A: W will hold the eigenvalues in
     * sorted (decreasing) order, V will contain the
     * eigenvectors as cols. eigen_ql() returns 1 on error
     */
    if (eigen_ql(Ata, W, V)) return(1);
    
    /* take the sqrt of the elements in W */
    for (unsigned j=0; j<C; j++)
	W[j]=(W[j]<0.0)? 0.0: sqrt(W[j]);
    
    /* get the matrix U: Ro'zsa says that A*v(j)=W[j]*u(j),
     * where u(j) and v(j) are the j-th columns of U and V, 
     * respectively. If W[j] is too small, then the corresponding
     * u(j) will be set to 0 (cf. sqrt(W) above)
     */
    Vector_ Uj(R);
    for (j=0; j<C; j++)
    {
	if (W[j]==0.0) continue;    // leave u(j) zeroed
	Uj=(A*V.col(j))/W[j];
	if (Rorig<R) Uj.dim(R);	// expand dims if Row<Col
	U.col(Uj, j);
    }
    
    return(0);	// OK
}
// END of make_decomp()

/* rank_cond(): checks the N singular values W[] of a matrix 
 * after SVD. The condition number Cond
 * (ratio of the largest and smallest singular value) is also
 * calculated. The singular values which are smaller than
 * Eps times the largest are set to 0.0.
 * Return value: the rank of the matrix.
 */
unsigned int Siva_::rank_cond(double Eps, double& Cond)
{
    double Wmax=-HUGE_VAL, Wmin=HUGE_VAL;
    unsigned int i;
    
    /* get the largest and smallest singular value */
    for (i=0; i<C; i++)
    {
	if (W[i]>Wmax) Wmax=W[i];
	if (W[i]<Wmin) Wmin=W[i];
    }
    
    /* calc the condition number: set to HUGE_VAL if Wmin==0.0 */
    Cond=(Wmin==0.0)? HUGE_VAL: Wmax/Wmin;
    
    /* set all singular values which are smaller than Eps*Wmax
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
 * and the weight vector should be "conditioned" (small entries
 * zeroed) by rank_cond() prior to the call to this routine.
 * Return value: the X vector.
 */
Vector_ Siva_::lin_solve(const Vector_& B) const
{
    // dim mismatches
    unsigned int Bdim=B.dim();
    if (R==Rorig && Bdim!=R || R>Rorig && !(Bdim==R || Bdim==Rorig))
    {
	cerr<<"? Siva_::lin_solve: B has wrong dimensions, nullvector returned\n";
	Vector_ X(C);	// a null-vector
	return(X);
    }
    
    Vector_ Wub(C);	// this will be W*U'*B
    // the "fewer eqns than unknowns" case
    if (R>Rorig && Bdim==Rorig)
    {
	cerr<<"? Siva_::lin_solve: fewer eqns ("<<Rorig
	    <<") than unknowns ("<<R<<"), B zero-padded\n";
	Vector_ Bpad(B); Bpad.dim(R);   // pad B w/ zeroes
	Wub=U.get_transpose()*Bpad;
    }
    else	// normal cases
	Wub=U.get_transpose()*B;
    
    for (unsigned int j=0; j<C; j++)
	if (W[j]==0.0) Wub[j]=0.0; else Wub[j]/=W[j];
    
    // get solution vector and return
    Vector_ X=V*Wub;
    return(X);
}
// END of lin_solve()

// ---- GLOBAL FUNCTIONS ----

/* <<: the overloaded output operator. Prints a neat listing to Out */
ostream& operator<<(ostream& Out, const Siva_& Svd)
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

// ==== END OF FUNCTIONS Siva.c++ ====
