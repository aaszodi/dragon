// ==== FUNCTIONS Trimat.c++ ====

/* Triangular matrix class. Symmetric square matrices are
 * stored as lower triangles (main diagonal included)
 * in this class to save memory.
 */

// SGI C++ 4.0, IRIX 5.3, 25. Oct. 1995. Andris

// ---- HEADER ----

#include "Trimat.h"

// ==== Trimat_ MEMBER FUNCTIONS ====

// ---- Constructors ----

/* Init to SizexSize (default 3x3), fill up with 0. */
Trimat_::Trimat_(unsigned int Size) : Matbase_(Size)
{
    if (!Size) Size=3;
    Eno=(Size*(Size+1))/2;
    
    // make element array: Rows[] were made by the Matbase_ ctor
    double *Elems=alloc_elems(Eno);
    if (Elems==NULL) return;
    
    memset(Elems, 0, Eno*sizeof(double));  // init to 0.0
	
    Rows[0]=Elems; init_rowptrs(Size);
}

/* Allocate a Row*Row triangular matrix and initialise it with
 * a traditional double ** style array Arr. It is the caller's 
 * responsibility to make sure that Arr has the correct dimensions.
 * Arr is supposed to have the "triangular" layout already but
 * there's no way of checking it. 
 * If Arr==NULL then the matrix elements will be set to 0.0
 */
Trimat_::Trimat_(const double **Arr, unsigned int Row)
{
    if (!Row) Row=3;	// zero rows/cols not accepted
    
    Eno=(Row*(Row+1))/2;  // store matrix dimensions and no. of elements
    
    // make element array: Rows[] were made by the Matbase_ ctor
    register double *Elems=alloc_elems(Eno);
    if (Elems==NULL) return;
    Rows[0]=Elems; init_rowptrs(Row);  
    
    // init elements
    register unsigned int i, j;
    if (Arr==NULL) memset(Elems, 0, Eno*sizeof(double));
    else 
	for (i=0; i<Row; i++)
	    for (j=0; j<=i; j++)
		memcpy(Rows[i], Arr[i], (i+1)*sizeof(double));
}
// END OF Rectbase_(Arr,Row,Col)

/* The copy constructor */
Trimat_::Trimat_(const Trimat_& Tri) : Matbase_(Tri.rno()) // 5-Jul-1997. g++
{
    Eno=(R*(R+1))/2;
    
    // make element array: Rows[] were made by the Matbase_ ctor
    register double *Elems=alloc_elems(Eno);
    if (Elems==NULL) return;
    
    memcpy(Elems, Tri.Rows[0], Eno*sizeof(double)); // copy elements	
    Rows[0]=Elems; init_rowptrs(R);  
}

/* Tri<--Sq conversion constructor. Inits with Sq's lower triangle */
Trimat_::Trimat_(const Sqmat_& Sq) : Matbase_(Sq)
{
    Eno=(R*(R+1))/2;
    
    // make element array: Rows[] were made by the Matbase_ ctor
    register double *Elems=alloc_elems(Eno);
    if (Elems==NULL) return;
    Rows[0]=Elems; init_rowptrs(R);
    
    // copy lower triangle and main diag
    register unsigned int i, j;
    for (i=0; i<R; i++)
	for (j=0; j<=i; j++)
	    Rows[i][j]=Sq[i][j];
}

// ---- Assignment ----

Trimat_& Trimat_::operator=(const Trimat_& Tri)
{
    if (this==&Tri) return(*this);  // Tri=Tri; nothing is done
    
    // if the dimensions happen to match, then no realloc is needed
    if (rno()!=Tri.rno())
    {
	// tough luck! must alloc new arrays
	double **Newrows=alloc_rows(Tri.rno());
	if (Newrows==NULL) return(*this);
	Newrows[0]=alloc_elems(Tri.Eno);
	if (Newrows[0]==NULL) { delete [] Newrows; return(*this); }
	
	// OK, we've got space, go ahead
	delete [] Rows[0]; delete [] Rows;	// destroy original
	Rows=Newrows;	// get new array(s)
	R=Tri.rno(); Eno=Tri.Eno;	// adjust new dimensions
	init_rowptrs(R);	// adjust new row pointers
    }

    // just copy elements
    memcpy(Rows[0], Tri.Rows[0], Eno*sizeof(double));
    return(*this);
}
// END of operator= 

// ---- Conversion ----

/* Converts a Trimat_ object into a Sqmat_ object, i.e.
 * makes a full symmetric matrix out of the sparse representation.
 */
Trimat_::operator Sqmat_() const
{
    Sqmat_ Sq(R);
    for (register unsigned int i=0; i<R; i++)
    {
	Sq[i][i]=Rows[i][i];	// copy main diag
	for (unsigned int j=0; j<i; j++)
	    Sq[i][j]=Sq[j][i]=Rows[i][j];   // copy lower triangle
    }
    return(Sq);
}
// END of Sqmat_ <-- Trimat_

/* get_array(): converts the calling object to a traditional
 * double ** array and sets Row and Col to its corresponding
 * dimensions. Returns a NULL ptr if allocation failed.
 * The returned array will retain the "triangular" layout
 * (similar to my earlier C "Trimat_" type).
 */
double** Trimat_::get_array(unsigned int& Row, unsigned int& Col) const
{
    register double **Arr=NULL;
    register unsigned int i, j;
    
    Row=Col=rno();   // copy dimensions
    Arr=new double* [rno()];	// alloc array of row ptrs
    if (Arr==NULL) { prt_err(NO_MEM, "get_array()"); return(NULL); }
    for (i=0; i<rno(); i++)
    {
	Arr[i]=new double [i+1];    // prepare for the triangular layout
	if (Arr[i]==NULL)	// getting hyper-paranoid
	{
	    prt_err(NO_MEM, "get_array()");
	    for (j=0; j<i; j++) delete [] Arr[j];
	    delete [] Arr; return(NULL);
	}
	memcpy(Arr[i], Rows[i], (i+1)*sizeof(double));
    }
    return(Arr);
}
// END of get_array()

// ---- Size ----

/* set_size(): Sets the size to Size x Size. If the new size is
 * larger than the old then the new rows will be set to 0.0, 
 * if less then the upper triangle will be preserved and the extra
 * rows will be lost. Zero size isn't allowed (no action).
 */
void Trimat_::set_size(unsigned int Size)
{
    if (!Size || R==Size) return;	    // no change
    
    // alloc new storage
    double **Newrows=alloc_rows(Size);
    if (Newrows==NULL) return;	// out-of-memory
    unsigned int Elemno=(Size*(Size+1))/2;
    register double *Elems=alloc_elems(Elemno);
    if (Elems==NULL) { delete [] Newrows; return; }
    
    // copy elements
    memset(Elems, '\0', Elemno*sizeof(double));	// zero
    unsigned int Pres=(R<Size)? R: Size;
    Pres=(Pres*(Pres+1))/2;	// no. of elements to be preserved
    memcpy(Elems, Rows[0], Pres*sizeof(double));
    
    // replace old with new
    delete [] Rows[0]; delete [] Rows;
    Rows=Newrows; Rows[0]=Elems;
    R=Size; Eno=Elemno; init_rowptrs(R);
}
// END of set_size()

// ---- Access ----

/* Safe access: the overlaid function call operator is used for this.
 * Error messages are printed if one of the two indices is out of
 * range and the index in question will be replaced by 0.
 * If the column index is higher than the row index (j>i)
 * then the two will be swapped (Trimat_s are symm square matrices).
 * Note that there are two versions of this function; one for
 * const access, the other for modifiable lvalues.
 */
double Trimat_::operator()(unsigned int i, unsigned int j) const
{
    if (i>=R) { prt_err(BAD_ROWRANGE, "(i, j)"); i=0; }
    if (j>=R) { prt_err(BAD_COLRANGE, "(i, j)"); j=0; }
    return((i>=j)? Rows[i][j]: Rows[j][i]); // lower triangle!
}

double& Trimat_::operator()(unsigned int i, unsigned int j)
{
    if (i>=R) { prt_err(BAD_ROWRANGE, "(i, j)"); i=0; }
    if (j>=R) { prt_err(BAD_COLRANGE, "(i, j)"); j=0; }
    return((i>=j)? Rows[i][j]: Rows[j][i]); // lower triangle!
}

/* row(Idx): returns the Idx-th row as a Vector_ object. The row
 * is "logical", i.e. the row of the symmetric square matrix 
 * corresponding to Trimat. Physically, the Idx-th row is copied
 * up to the main diagonal, then the Idx-th column down from the
 * main diag completes the operation.
 * row(Vec, Idx): sets the elements in the Idx-th logical row to the elements
 * of the vector Vec, provided Idx is legal and Vec has the right
 * number of elements.
 * Note that the col() method simply calls the corresponding row()
 * because of the symmetry.
 */
Vector_ Trimat_::row(unsigned int Idx) const
{
    if (Idx>=R) { prt_err(BAD_ROWRANGE, "get_row(Idx)"); Idx=0; }
    Vector_ Rowvec(R);
    register unsigned int i;
    
    for (i=0; i<=Idx; i++) Rowvec[i]=Rows[Idx][i];  // copy row
    for (i=Idx+1; i<R; i++) Rowvec[i]=Rows[i][Idx]; // copy col
    return(Rowvec);
}
// END of row(Idx)

void Trimat_::row(const Vector_& Vec, unsigned int Idx)
{
    if (Vec.dim()!=R) 
	{ prt_err(DIM_MISMATCH, "set_row(Idx)"); return; }
    if (Idx>=R)  { prt_err(BAD_ROWRANGE, "set_row(Idx)"); Idx=0; }
    register unsigned int i;
    
    for (i=0; i<=Idx; i++) Rows[Idx][i]=Vec[i];
    for (i=Idx+1; i<R; i++) Rows[i][Idx]=Vec[i];
}
// END of row(Vec,Idx)

/* col(Idx): returns the Idx-th column as a Vector_ object safely.
 * col(Vec,Idx): sets the elements in the Idx-th column to the elements
 * of the vector Vec, provided Idx is legal and Vec has the right
 * number of elements.
 */
Vector_ Trimat_::col(unsigned int Idx) const { return row(Idx); }
void Trimat_::col(const Vector_& Vec, unsigned int Idx) { row(Vec, Idx); }

// ---- Simple linear algebra ----

/* Matrix multiplication: the left operand is returned if
 * a dimension mismatch is detected. For simplicity, the
 * "safe indexing" is used here with the associated performance
 * penalty.
 * Returns a general RxC Matrix_.
 */
Matrix_ Trimat_::operator*(const Matbase_& Mat) const
{
    if (Mat.rno()!=rno())
    {
	prt_err(DIM_MISMATCH, "Tri*Mat");
	return(Sqmat_(*this));
    }

    Matrix_ Prod(rno(), Mat.cno());
    register unsigned int i, j, k;
    register double *Ps;
    register double Temp;
    
    for (i=0; i<rno(); i++)
    {
	Ps=Prod[i];
	for (j=0; j<rno(); j++)
	{
	    Temp=0.0;
	    for (k=0; k<cno(); k++) Temp+=(*this)(i, k)*Mat(k, j);
	    Ps[j]=Temp;
	}
    }
    return(Prod);
}
// END of Tri*Mat

/* Matrix*vector multiplication: the vector Vec is returned if a
 * dimension mismatch was detected.
 */
Vector_ Trimat_::operator*(const Vector_& Vec) const
{
    if (rno()!=Vec.dim()) { prt_err(DIM_MISMATCH, "Tri*Vec"); return(Vec); }
    Vector_ Prod(rno());
    register unsigned int i, j;
    register double Temp;
    
    for (i=0; i<rno(); i++)
    {
	Temp=0.0;
	for (j=0; j<=i; j++) Temp+=Rows[i][j]*Vec[j];
	for (j=i+1; j<cno(); j++) Temp+=Rows[j][i]*Vec[j];
	Prod[i]=Temp;
    }
    return(Prod);
}
// END of Tri*Vec

// ---- Arithmetics ----
    
/* NOTE: SGI Delta/C++ 4.0 is sensitive to the order of the
 * overloading of the multiplication operator. Although
 * it makes no sense, Trimat_::operator*(double) must follow
 * the others when compiling with DCC -smart.
 */

/* Matrix addition and subtraction (+,-) check for dimension
 * mismatches. If a mismatch is found then the left
 * operand is returned unchanged. Division by scalar ( / )
 * is always checked for division by zero and the unchanged matrix
 * is returned if a div by 0 is attempted. These routines are
 * based on the "in-place" operators (+=, -=, *=, /=) inherited
 * from Matbase_ (could not be declared there because of the
 * return type).
 */

Trimat_ Trimat_::operator+(const Trimat_& Mat) const
{
    if (R!=Mat.R) { prt_err(DIM_MISMATCH, "Mat+Mat"); return(*this); }
    Trimat_ Sum(*this);
    Sum+=Mat; return(Sum);
}

Trimat_ Trimat_::operator-(const Trimat_& Mat) const
{
    if (R!=Mat.R) { prt_err(DIM_MISMATCH, "Mat-Mat"); return(*this); }
    Trimat_ Dif(*this);
    Dif-=Mat; return(Dif);
}

Trimat_ Trimat_::operator*(double Factor) const
{
    Trimat_ Mat(*this);
    Mat*=Factor; return(Mat);
}

Trimat_ operator*(double Factor, const Trimat_& Tri)	// friend
{
    Trimat_ Mat(Tri);
    Mat*=Factor; return(Mat);
}

Trimat_ Trimat_::operator/(double Div) const
{
    if (Div==0.0) { prt_err(DIV_BY_ZERO, "Mat/Scal"); return(*this); }
    Trimat_ Mat(*this);
    Mat/=Div; return(Mat);
}

// ---- Printing ----

/* print_rows(): print the rows of a Trimat so that
 * if the matrix is big then it is cut up into chunks nicely.
 * Called by list_matrix() only: see Matbase_::list_matrix() for
 * a description of the parameters.
 */
void Trimat_::print_rows(ostream& Out, unsigned int Sizew, unsigned int Jbeg, 
	unsigned int Items, unsigned int Width, unsigned int Prec) const
{
    for (unsigned int i=0; i<rno(); i++)
    {
	Out.unsetf(ios::left);
	Out << setw(Sizew) << i << setw(0) << " | " ;   // row idx
	Out.setf(ios::left);
	/* list a row: either all items or Items which fit on one line */
	for (unsigned int j=Jbeg; j<=i && j<Jbeg+Items; j++)
	    Out << setw(Width) << setprecision(Prec) << (*this)(i, j);
	Out << endl;
    }
}
// END of print_rows()

// ---- Memory management ----

/* init_rowptrs(): initialises the row pointers of a Trimat_
 * object. It is assumed that Rows[0] has already been set to
 * point to the element array and that all other data members
 * are initialised. Called by ctors and the copy operator.
 */
void Trimat_::init_rowptrs(unsigned int Colno)
{
    if (Rows[0]==NULL) 
    {
	cerr<<"? Trimat_::init_rowptrs(): Cannot init (Rows[0]==NULL)\n";
	return;
    }
    for (register unsigned int i=1; i<Colno; i++)	    // init row index array
	Rows[i]=Rows[i-1]+i;	// longer and longer rows in lower triangle
}
// END of init_rowptrs()

// ==== END OF FUNCTIONS Trimat.c++ ====
