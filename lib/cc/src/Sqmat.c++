// ==== FUNCTIONS Sqmat.c++ ====

/* Double-precision square matrix class. Derived from the
 * general RxC rectangular matrix class Rectbase_ and the
 * square matrix base class Sqbase_ . 
 */

// SGI C++ 4.0, IRIX 5.3, 25. Oct. 1995. Andris 

// ---- HEADER ----

#include "Sqmat.h"

// ==== Sqmat_ MEMBER FUNCTIONS ====

// ---- Constructors ----

/* Init with a Rectbase_ (RxC) object (conversion ctor). The resulting Sqmat_
 * will contain the full matrix so that its size will be max(R, C)
 * with the extra rows/cols padded with zeroes. Also works with Sqmat_ .
 */
Sqmat_::Sqmat_(const Rectbase_& Rbase):
    Matbase_(Rbase.max_size()), Rectbase_(Rbase.max_size())
{
    register unsigned int i, j;
    for (i=0; i<Rbase.rno(); i++)
	for (j=0; j<Rbase.cno(); j++)
	    Rows[i][j]=Rbase[i][j];
}

// ---- Assignment ----

Sqmat_& Sqmat_::operator=(const Sqmat_& Sq)
{
    if (this==&Sq) return(*this);  // Sq=Sq; nothing is done
    
    // if the dimensions happen to match, then no realloc is needed
    if (rno()!=Sq.rno())
    {
	// tough luck! must alloc new arrays
	double **Newrows=alloc_rows(Sq.rno());
	if (Newrows==NULL) return(*this);
	Newrows[0]=alloc_elems(Sq.Eno);
	if (Newrows[0]==NULL) { delete [] Newrows; return(*this); }
	
	// OK, we've got space, go ahead
	delete [] Rows[0]; delete [] Rows;	// destroy original
	Rows=Newrows;	// get new array(s)
	R=Sq.rno(); Eno=Sq.Eno;	// adjust new dimensions
	init_rowptrs(R);	// adjust new row pointers
    }

    // just copy elements
    memcpy(Rows[0], Sq.Rows[0], Eno*sizeof(double));
    return(*this);
}
// END of operator= 

// ---- Size ----

/* set_size(): resets the size to Size x Size. The upper left
 * corner overlap is preserved, if the new size is
 * less than the original then the extra rows and cols will be lost.
 * If new rows/cols are added then they will be initialised to 0.0.
 * Zero row/col numbers not permitted.
 * If the new size is 0 or the same as the old then no action is taken.
 */
void Sqmat_::set_size(unsigned int Size)
{
    if (!Size || rno()==Size) return;   // no change
    
    // allocate new arrays
    double **Newrows=alloc_rows(Size);
    if (Newrows==NULL) return;	// out of memory
    unsigned int Rc=Size*Size;
    Newrows[0]=alloc_elems(Rc);
    if (Newrows[0]==NULL) { delete [] Newrows; return; }
    
    // copy elements
    memset(Newrows[0], '\0', Rc*sizeof(double));	// zero everything

    register unsigned int i, Rmin=(R>Size)? Size: R;
    for (i=0; i<Rmin; i++)
    {
	if (i) Newrows[i]=Newrows[i-1]+Size; // init rowptrs on the fly
	memcpy(Newrows[i], Rows[i], Rmin*sizeof(double));   // copy row
    }
    for (; i<Size; i++) Newrows[i]=Newrows[i-1]+Size; // finish up rowptrs
    
    // replace old with new
    delete [] Rows[0]; delete [] Rows;
    Rows=Newrows; R=Size; Eno=Rc;
}
// END of set_size()

// ---- Simple linear algebra ----

/* Matrix multiplication: the left operand is returned if
 * a dimension mismatch is detected. The code is exactly the same
 * as the corresponding Matrix_ operator; however, the return type
 * does not allow us to decl/def the whole thing in Rectbase_
 * (where it logically belongs).
 */
Matrix_ Sqmat_::operator*(const Rectbase_& Mat) const
{
    if (Mat.rno()!=cno()) { prt_err(DIM_MISMATCH, "Mat*Mat"); return(*this); }
    Matrix_ Prod(rno(), Mat.cno());
    register unsigned int i, j, k;
    register double *Ps, *Rs;
    register double Temp;
    
    for (i=0; i<rno(); i++)
    {
	Ps=Prod[i]; Rs=Rows[i];	// fixed row pointers in this cycle
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
Vector_ Sqmat_::operator*(const Vector_& Vec) const
{
    if (R!=Vec.dim()) { prt_err(DIM_MISMATCH, "Sq*Vec"); return(Vec); }
    Vector_ Prod(R);
    register unsigned int i, j;
    register double *Rs;
    register double Temp;
    
    for (i=0; i<R; i++)
    {
	Temp=0.0; Rs=Rows[i];
	for (j=0; j<R; j++) Temp+=Rs[j]*Vec[j];	// quick vector access
	Prod[i]=Temp;
    }
    return(Prod);
}
// END of Mat*Vec

// ---- Arithmetics ----
    
/* NOTE: SGI Delta/C++ 4.0 is sensitive to the order of the
 * overloading of the multiplication operator. Although
 * it makes no sense, Sqmat_::operator*(double) must follow
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

Sqmat_ Sqmat_::operator+(const Sqmat_& Mat) const
{
    if (rno()!=Mat.rno()) { prt_err(DIM_MISMATCH, "Sqmat+Sqmat"); return(*this); }
    Sqmat_ Sum(*this);
    Sum+=Mat; return(Sum);
}

Sqmat_ Sqmat_::operator-(const Sqmat_& Mat) const
{
    if (rno()!=Mat.rno()) { prt_err(DIM_MISMATCH, "Sqmat-Sqmat"); return(*this); }
    Sqmat_ Dif(*this);
    Dif-=Mat; return(Dif);
}

Sqmat_ Sqmat_::operator*(double Factor) const
{
    Sqmat_ Mat(*this);
    Mat*=Factor; return(Mat);
}

Sqmat_ operator*(double Factor, const Sqmat_& Mat)	// this is a friend
{
    Sqmat_ Matr(Mat);
    Matr*=Factor; return(Matr);
}

Sqmat_ Sqmat_::operator/(double Div) const
{
    if (Div==0.0) { prt_err(DIV_BY_ZERO, "Sqmat/Scal"); return(*this); }
    Sqmat_ Mat(*this);
    Mat/=Div; return(Mat);
}

// ---- Transpose ----

/* get_transpose(): returns the transposed calling object. */
Sqmat_ Sqmat_::get_transpose() const
{
    Sqmat_ Tr(*this);
    Tr.transpose_inplace();
    return(Tr);
}
// END of get_transpose()

/* transpose_inplace(): transposes the calling object "in place". */
void Sqmat_::transpose_inplace()
{
    register double Temp;
    register unsigned int i, j;
    for (i=0; i<R; i++)
	for (j=0; j<i; j++)
	{
	    Temp=Rows[i][j]; Rows[i][j]=Rows[j][i]; Rows[j][i]=Temp;
	}
}
// END of transpose_inplace()

// ==== END OF FUNCTIONS Sqmat.c++ ====
