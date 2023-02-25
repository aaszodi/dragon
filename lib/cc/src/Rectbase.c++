// ==== FUNCTIONS Rectbase.c++ ====

/* Abstract base class for rectangular matrices (square and non-square). */
 
// SGI C++ 4.0, IRIX 5.3, 4. May 1995. Andris

// ----HEADER ----

#include "Rectbase.h"

// ==== Rectbase_ MEMBER FUNCTIONS ====

// ---- Constructors ----

/* Allocates a Row x Col matrix. Without any arguments, a 3x3
 * matrix will be created. With one argument, the Col parameter
 * defaults to 0,  which is interpreted as an instruction to create
 * a *square* Row x Row matrix. Square matrices may be created
 * by Rectbase_(N, N) or Rectbase_(N) or Rectbase_(N, 0).
 * All elements will be initialised to 0.0.
 */
Rectbase_::Rectbase_(unsigned int Row, unsigned int Col)
{
    if (!Row) Row=3;	// zero rows not accepted
    if (!Col) Col=Row;	// 0 columns->RowxRow matrix

    Eno=Row*Col;  // store matrix dimensions and no. of elements
    
    // make element array: Rows[] were made by the Matbase_ ctor
    double *Elems=alloc_elems(Eno);
    if (Elems==NULL) return;
    
    memset(Elems, 0, Eno*sizeof(double));  // init elements to 0.0 (usually works)
	
    Rows[0]=Elems; init_rowptrs(Col);
}
// END OF Rectbase_(Row,Col)

/* Allocate a Row*Col matrix as above and initialise it with
 * a traditional double ** style array Arr. It is the caller's 
 * responsibility to make sure that Arr has the correct dimensions.
 * If Arr==NULL then the matrix elements will be set to 0.0
 */
Rectbase_::Rectbase_(const double **Arr, unsigned int Row, unsigned int Col)
{
    if (!Row) Row=3;	// zero rows/cols not accepted
    if (!Col) Col=Row;
    
    Eno=Row*Col;  // store matrix dimensions and no. of elements
    
    // make element array: Rows[] were made by the Matbase_ ctor
    double *Elems=alloc_elems(Eno);
    if (Elems==NULL) return;
    Rows[0]=Elems; init_rowptrs(Col);  
    
    // init elements
    register unsigned int i;
    if (Arr==NULL) memset(Elems, 0, Eno*sizeof(double));
    else 
	for (i=0; i<Row; i++) memcpy(Rows[i], Arr[i], Col*sizeof(double));
}
// END OF Rectbase_(Arr,Row,Col)

/* the copy constructor */
Rectbase_::Rectbase_(const Rectbase_& Mat)
{
    Eno=Mat.Eno;	// set sizes 
    
    // make element array: Rows[] were made by the Matbase_ ctor
    double *Elems=alloc_elems(Eno);
    if (Elems==NULL) return;
    Rows[0]=Elems; init_rowptrs(Mat.cno());  
    
    // just copy element array
    memcpy(Elems, Mat.Rows[0], Eno*sizeof(double));
}
// END of Rectbase_(Mat)

// ---- Access ----

/* Safe access: the overlaid function call operator is used for this.
 * Error messages are printed if one of the two indices is out of
 * range and the index in question will be replaced by 0.
 * Note that there are two versions of this function; one for
 * const access, the other for modifiable lvalues.
 */
double Rectbase_::operator()(unsigned int i, unsigned int j) const
{
    if (i>=rno()) { prt_err(BAD_ROWRANGE, "(i, j)"); i=0; }
    if (j>=cno()) { prt_err(BAD_COLRANGE, "(i, j)"); j=0; }
    return(Rows[i][j]);
}
// END OF operator() (const access)

double& Rectbase_::operator()(unsigned int i, unsigned int j)
{
    if (i>=rno()) { prt_err(BAD_ROWRANGE, "(i, j)"); i=0; }
    if (j>=cno()) { prt_err(BAD_COLRANGE, "(i, j)"); j=0; }
    return(Rows[i][j]);
}
// END OF operator() (lvalue access)

/* row(Idx): returns the Idx-th row as a Vector_ object safely.
 * row(Vec,Idx): sets the elements in the Idx-th row to the elements
 * of the vector Vec, provided Idx is legal and Vec has the right
 * number of elements.
 */
Vector_ Rectbase_::row(unsigned int Idx) const
{
    if (Idx>=rno()) { prt_err(BAD_ROWRANGE, "row(Idx)"); Idx=0; }
    Vector_ Rowvec(Rows[Idx], cno());   // use the traditional array ctor
    return(Rowvec);
}
// END of row(Idx)

void Rectbase_::row(const Vector_& Vec, unsigned int Idx)
{
    if (Vec.dim()!=cno()) 
	{ prt_err(DIM_MISMATCH, "row(Vec, Idx)"); return; }
    if (Idx>=rno())  { prt_err(BAD_ROWRANGE, "row(Vec, Idx)"); Idx=0; }
    register double *Row=Rows[Idx];
    register unsigned int i;
    for (i=0; i<cno(); i++) Row[i]=Vec[i];
}
// END of row(Vec,Idx)

/* col(Idx): returns the Idx-th column as a Vector_ object safely.
 * col(Vec,Idx): sets the elements in the Idx-th column to the elements
 * of the vector Vec, provided Idx is legal and Vec has the right
 * number of elements.
 */
Vector_ Rectbase_::col(unsigned int Idx) const
{
    if (Idx>=cno()) { prt_err(BAD_COLRANGE, "col(Idx)"); Idx=0; }
    Vector_ Colvec(rno());
    register unsigned int i;
    for (i=0; i<rno(); i++) Colvec[i]=Rows[i][Idx];
    return(Colvec);
}
// END of col(Idx)

void Rectbase_::col(const Vector_& Vec, unsigned int Idx)
{
    if (Vec.dim()!=rno()) 
	{ prt_err(DIM_MISMATCH, "col(Vec, Idx)"); return; }
    if (Idx>=cno())  { prt_err(BAD_ROWRANGE, "col(Vec, Idx)"); Idx=0; }
    register unsigned int i;
    for (i=0; i<rno(); i++) Rows[i][Idx]=Vec[i];
}
// END of col(Vec,Idx)

/* get_array(): converts the calling object to a traditional
 * double ** array and sets Row and Col to its corresponding
 * dimensions. Returns a NULL ptr if allocation failed.
 */
double** Rectbase_::get_array(unsigned int& Row, unsigned int&Col) const
{
    register double **Arr=NULL;
    register unsigned int i, j;
    
    Row=rno(); Col=cno();   // copy dimensions
    Arr=new double* [rno()];
    if (Arr==NULL) { prt_err(NO_MEM, "get_array()"); return(NULL); }
    for (i=0; i<rno(); i++)
    {
	Arr[i]=new double [cno()];
	if (Arr[i]==NULL)	// getting hyper-paranoid
	{
	    prt_err(NO_MEM, "get_array()");
	    for (j=0; j<i; j++) delete [] Arr[j];
	    delete [] Arr; return(NULL);
	}
	memcpy(Arr[i], Rows[i], cno()*sizeof(double));
    }
    return(Arr);
}
// END of get_array()

// ---- Printing ----

/* print_rows(): print the rows of a rectangular matrix so that
 * if the matrix is big then it is cut up into chunks nicely.
 * Called by list_matrix() only: see Matbase_::list_matrix() for
 * a description of the parameters.
 */
void Rectbase_::print_rows(ostream& Out, unsigned int Sizew, unsigned int Jbeg, 
	unsigned int Items, unsigned int Width, unsigned int Prec) const
{
    for (unsigned int i=0; i<rno(); i++)
    {
	Out.unsetf(ios::left);
	Out << setw(Sizew) << i << setw(0) << " | " ;   // row idx
	Out.setf(ios::left);
	/* list a row: either all items or Items which fit on one line */
	for (unsigned int j=Jbeg; j<cno() && j<Jbeg+Items; j++)
	    Out << setw(Width) << setprecision(Prec) << (*this)(i, j);
	Out << endl;
    }
}
// END of print_rows()

// ---- Memory management ----

/* init_rowptrs(): initialises the row pointers of a Rectbase_
 * object. It is assumed that Rows[0] has already been set to
 * point to the element array which will contain Colno columns.
 * are initialised. Called by ctors and the copy operator.
 */
void Rectbase_::init_rowptrs(unsigned int Colno)
{
    if (Rows[0]==NULL) 
    {
	cerr<<"? Rectbase_::init_rowptrs(): Cannot init (Rows[0]==NULL)\n";
	return;
    }
    for (register unsigned int i=1; i<rno(); i++)	    // init row index array
	Rows[i]=Rows[i-1]+Colno;
}
// END of init_rowptrs()

// ==== END OF FUNCTIONS Rectbase.c++ ====
