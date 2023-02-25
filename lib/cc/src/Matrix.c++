// ==== FUNCTIONS Matrix.c++ ====

/* Double-precision class for general RxC rectangular matrices. 
 * The Vector_ class is used for various linear algebra operations.
 */
 
// SGI C++ 4.0, IRIX 5.3, 25. Oct. 1995. Andris

// ----HEADER ----

#include "Matrix.h"

// ==== Matrix_ MEMBER FUNCTIONS ====

// ---- Assignment ----

Matrix_& Matrix_::operator=(const Matrix_& Mat)
{
    if (this==&Mat) return(*this);  // Mat=Mat; nothing is done
    
    // if the dimensions happen to match, then no realloc is needed
    if (rno()!=Mat.rno() || cno()!=Mat.cno())
    {
	// tough luck! must alloc new arrays
	double **Newrows=alloc_rows(Mat.rno());
	if (Newrows==NULL) return(*this);
	Newrows[0]=alloc_elems(Mat.Eno);
	if (Newrows[0]==NULL) { delete [] Newrows; return(*this); }
	
	// OK, we've got space, go ahead
	delete [] Rows[0]; delete [] Rows;	// destroy original
	Rows=Newrows;	// get new array(s)
	R=Mat.rno(); C=Mat.cno(); Eno=Mat.Eno;	// adjust new dimensions
	init_rowptrs(C);	// adjust new row pointers
    }

    // just copy elements
    memcpy(Rows[0], Mat.Rows[0], Eno*sizeof(double));
    return(*this);
}
// END of operator= 

// ---- Size ----

/* set_size(): resets the size to Rno x Cno. The upper left
 * corner overlap is preserved, if the new row/col number is
 * less than the original then the extra rows/cols will be lost.
 * If new rows/cols are added then they will be initialised to 0.0.
 * Zero row/col numbers not permitted.
 * If the new size is 0 or the same as the old then no action is taken.
 */
void Matrix_::set_size(unsigned int Rno, unsigned int Cno)
{
    if (!Rno || !Cno || rno()==Rno && cno()==Cno) return;   // no change
    
    // allocate new arrays
    double **Newrows=alloc_rows(Rno);
    if (Newrows==NULL) return;	// out of memory
    unsigned int Rc=Rno*Cno;
    Newrows[0]=alloc_elems(Rc);
    if (Newrows[0]==NULL) { delete [] Newrows; return; }
    
    // copy elements
    memset(Newrows[0], '\0', Rc*sizeof(double));	// zero everything
    register unsigned int i, 
	Rmin=(R>Rno)? Rno: R, Cmin=(C>Cno)? Cno: C;
    for (i=0; i<Rmin; i++)
    {
	if (i) Newrows[i]=Newrows[i-1]+Cno; // init rowptrs on the fly
	memcpy(Newrows[i], Rows[i], Cmin*sizeof(double));   // copy row
    }
    for (; i<Rno; i++) Newrows[i]=Newrows[i-1]+Cno; // finish up rowptrs
    
    // replace old with new
    delete [] Rows[0]; delete [] Rows;
    Rows=Newrows; R=Rno; C=Cno; Eno=Rc;
}
// END of set_size()

// ---- Simple linear algebra ----

/* Matrix multiplication: the left operand is returned if
 * a dimension mismatch is detected. There is no in-place (*=) version.
 */
Matrix_ Matrix_::operator*(const Rectbase_& Mat) const
{
    if (Mat.rno()!=cno()) { prt_err(DIM_MISMATCH, "Mat*Mat"); return(*this); }
    Matrix_ Prod(rno(), Mat.cno());
    register unsigned int i, j, k;
    register double *Ps, *Rs;
    register double Temp;
    
    for (i=0; i<rno(); i++)
    {
	Ps=Prod.Rows[i]; Rs=Rows[i];	// fixed row pointers in this cycle
	for (j=0; j<Mat.cno(); j++)
	{
	    Temp=0.0;
	    for (k=0; k<cno(); k++) Temp+=Rs[k]*Mat[k][j];
	    Ps[j]=Temp;
	}
    }
    return(Prod);
}
// END of Mat*Mat

/* Matrix*vector multiplication: the vector is returned if a
 * dimension mismatch was detected.
 */
Vector_ Matrix_::operator*(const Vector_& Vec) const
{
    if (cno()!=Vec.dim()) { prt_err(DIM_MISMATCH, "Mat*Vec"); return(Vec); }
    Vector_ Prod(rno());
    register unsigned int i, j;
    register double *Rs;
    register double Temp;
    
    for (i=0; i<rno(); i++)
    {
	Temp=0.0; Rs=Rows[i];
	for (j=0; j<cno(); j++) Temp+=Rs[j]*Vec[j];	// quick vector access
	Prod[i]=Temp;
    }
    return(Prod);
}
// END of Mat*Vec

// ---- Arithmetics ----

/* NOTE: SGI Delta/C++ 4.0 is sensitive to the order of the
 * overloading of the multiplication operator. Although
 * it makes no sense, Matrix_::operator*(double) must follow
 * the others when compiling with DCC -smart.
 */

/* Matrix addition and subtraction (+,-,+=,-=) check for dimension
 * mismatches. If a mismatch is found then the left
 * operand is returned unchanged. Division by scalar ( / )
 * is always checked for division by zero and the unchanged matrix
 * is returned if a div by 0 is attempted. *=, /= are inherited
 * from Matbase_ (could not be declared there because of the
 * return type).
 */
Matrix_ Matrix_::operator+(const Rectbase_& Mat) const
{
    if (rno()!=Mat.rno() || cno()!=Mat.cno()) { prt_err(DIM_MISMATCH, "Mat+Mat"); return(*this); }
    Matrix_ Sum(*this);
    Sum+=Mat; return(Sum);
}

Matrix_ Matrix_::operator-(const Rectbase_& Mat) const
{
    if (rno()!=Mat.rno() || cno()!=Mat.cno()) { prt_err(DIM_MISMATCH, "Mat-Mat"); return(*this); }
    Matrix_ Dif(*this);
    Dif-=Mat; return(Dif);
}

Matrix_ Matrix_::operator*(double Factor) const
{
    Matrix_ Mat(*this);
    Mat*=Factor; return(Mat);
}

Matrix_ operator*(double Factor, const Rectbase_& Mat)	// this is a friend
{
    Matrix_ Matr(Mat);
    Matr*=Factor; return(Matr);
}

Matrix_ Matrix_::operator/(double Div) const
{
    if (Div==0.0) { prt_err(DIV_BY_ZERO, "Mat/Scal"); return(*this); }
    Matrix_ Mat(*this);
    Mat/=Div; return(Mat);
}

// ---- Transpose ----

/* get_transpose(): returns the transpose of the calling object. */
Matrix_ Matrix_::get_transpose() const
{
    Matrix_ Tr(cno(), rno());
    register unsigned int i, j;
    
    for (i=0; i<rno(); i++)
	for (j=0; j<cno(); j++)
	    Tr[j][i]=Rows[i][j];
    return(Tr);
}
// END of get_transpose()

// ==== END OF FUNCTIONS Matrix.c++ ====
