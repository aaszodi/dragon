// ==== FUNCTIONS Matbase.c++ ====

/* Abstract base class for rectangular and triangular matrices. */

// SGI C++ 4.0, IRIX. 5.3, 25. Oct. 1995. Andris

// ---- HEADER ----

#include "Matbase.h"

// ==== Matbase_ MEMBER FUNCTIONS ====

// ---- Constructors ----

/* Init a matrix to have Row>=1 rows. The element array
 * is not allocated and Rows[] is not initialised.
 * The element number is set to 0.
 */
Matbase_::Matbase_(unsigned int Row): Eno(0)
{
    if (!Row) Row=3;
    Rows=alloc_rows(Row); R=Row; Rows[0]=NULL;
}

/* Copy constructor: allocates the Rows array to the size of Mb,
 * does not initialise it, sets R to Mb.rno(), sets Eno=0.
 */
Matbase_::Matbase_(const Matbase_& Mb): R(Mb.rno()), Eno(0)
{
    Rows=alloc_rows(R); Rows[0]=NULL;
}

// ---- Access ----

/* set_values(): sets all elements to a specified value Val.
 * Default is 0.0.
 */
void Matbase_::set_values(double Val)
{
    register unsigned int i;
    register double *Elems=Rows[0];
    for (i=0; i<Eno; i++) Elems[i]=Val;
}
// END of set_values()

// ---- FORTRAN Indexing ----

/* ftn_idx(): moves the data pointers so that 1..N
 * indexing could be used instead of 0..N-1. Note that
 * this operation is inherently unsafe (the "safe indexing"
 * method knows nothing about it!) and should be used only
 * when writing wrappers around Numerical Recipes code etc.
 * and only if the unsafe [][] indexing is used.
 * Make sure the object abused by this method is not destroyed
 * and indexing is reset via c_idx() as soon as possible.
 */
void Matbase_::ftn_idx()
{
    for (register unsigned int i=0; i<R; i++) Rows[i]--;
    Rows--;
}
// END of ftn_idx()

/* c_idx(): resets the ptrs to their normal C-style position
 * after a call to ftn_idx().
 */
void Matbase_::c_idx()
{
    Rows++;
    for (register unsigned int i=0; i<R; i++) Rows[i]++;
}
// END of c_idx()

// ---- Arithmetics ----

/* Matrix addition and subtraction (+,+=,-,-=) check for dimension
 * mismatches. If a mismatch is found then the left
 * operand is returned unchanged. Division by scalar ( /,  /=)
 * is always checked for division by zero and the unchanged matrix
 * is returned if a div by 0 is attempted. These routines operate
 * on all elements of the element array, regardless of the row/col
 * layout and therefore are defined in the base class here and
 * will simply be inherited.
 */

Matbase_& Matbase_::operator+=(const Matbase_& Mat)
{
    if (rno()!=Mat.rno() || cno()!=Mat.cno() || Eno!=Mat.Eno)
    { prt_err(DIM_MISMATCH, "Mat+=Mat"); return(*this); }
    double *E1=Rows[0], *E2=Mat.Rows[0];
    for (register unsigned int i=0; i<Eno; i++) E1[i]+=E2[i];
    return(*this);
}

Matbase_& Matbase_::operator-=(const Matbase_& Mat)
{
    if (rno()!=Mat.rno() || cno()!=Mat.cno() || Eno!=Mat.Eno)
    { prt_err(DIM_MISMATCH, "Mat-=Mat"); return(*this); }
    double *E1=Rows[0], *E2=Mat.Rows[0];
    for (register unsigned int i=0; i<Eno; i++) E1[i]-=E2[i];
    return(*this);
}

Matbase_& Matbase_::operator*=(double Factor)
{
    double *Elem=Rows[0];
    for (register unsigned int i=0; i<Eno; i++) Elem[i]*=Factor;
    return(*this);
}

Matbase_& Matbase_::operator/=(double Div)
{
    if (fabs(Div)<DBL_EPSILON) { prt_err(DIV_BY_ZERO, "Mat/=Scal"); return(*this); }
    double *Elem=Rows[0];
    Div=1.0/Div;    // "multiply-by-reciprocal" trick
    for (register unsigned int i=0; i<Eno; i++) Elem[i]*=Div;
    return(*this);
}

// ---- End of public methods interface ----

// ---- Memory management ----

/* alloc_rows(): Row pointer array allocation. Note that a local variable is
 * used with a different name to emphasise that no access is
 * made to data members since these are undefined at construction time.
 * Zero rows are not allowed.
 * Returns the array address or NULL if out of memory.
 */
double** Matbase_::alloc_rows(unsigned int Rowno)
{
    if (!Rowno) Rowno=1;
    double **Rs=new double* [Rowno];	// alloc array of row pointers
    if (Rs==NULL) 
    {
	prt_err(NO_MEM, "alloc_rows()"); 
	return(NULL);
    }
    Rs[0]=NULL;	    // highly paranoid
    return(Rs);
}
// END of alloc_rows()

/* alloc_elems(): Element array allocation. Allocates an Elemno>=1
 * long double array and returns its address or NULL if out of memory.
 */
double* Matbase_::alloc_elems(unsigned int Elemno)
{
    if (!Elemno) Elemno=1;
    double *Elems=new double [Elemno];	// alloc element array
    if (Elems==NULL) 
    {
	prt_err(NO_MEM, "alloc_elems()"); 
	return(NULL);
    }
    else return(Elems);
}
// END of alloc_elems()

// ---- Printing ----

/* list_matrix: lists calling object to stdout with entries occupying Width chars,
 * in scientific format, Prec digits precision. If a row takes up 
 * more than Linewidth chars, then the matrix is cut up nicely. 
 * If Width<Prec+8,  then it is adjusted but no warning is given.
 * Pure C++ I/O. See the redefined << in the Matbase.h header.
 * Could be called by all derived classes.
 */
void Matbase_::list_matrix(ostream& Out, 
	unsigned int Prec, unsigned int Width, 
	unsigned int Linewidth) const
{
    char Nwidth[10];  // buffer to determine the printed width of N
    ostrstream Ostr(Nwidth, sizeof(Nwidth));	// stream output buffer
    long Oldflags=Out.flags();	// store previous settings here
    int Oldprec=Out.precision();
    unsigned int j,k,N, Sizew,Chunk,Jbeg,Items,Ulinelen;

    if (Width<Prec+8) Width=Prec+8;	// adjust Width
    N=cno();
    Ostr << N << ends;	    	// get printed width of C
    Sizew=Ostr.pcount()-1;
    if (Sizew>Width) Width=Sizew;
    Items=(Linewidth-Sizew-3)/(Width+1);	// columns per chunk

    /* output will be right-justified by default */
    Out.unsetf(ios::right|ios::internal);
    Out.setf(ios::showpoint|ios::scientific);	// %e in C
    
    /* main cycle: print chunks of the matrix */    
    for (Jbeg=0,Chunk=(N-1)/Items+1; Chunk>0; Jbeg+=Items,Chunk--)
    {
	/* underline length */
	Ulinelen=(Chunk>1)? Items*(Width+1)+Sizew+3:
		(N-Jbeg)*(Width+1)+Sizew+3;
		
	/* print column heading */
	for (k=0; k<Sizew+3; k++) Out.put(' ');  // leading spaces
	
	/* print col indices, left-justified, Width width */
	Out.setf(ios::left);
	for (j=Jbeg; j<N && j<Jbeg+Items; j++)
	    Out << setw(Width) << j;
	Out << endl;
	
	/* underline column heading */
	for (k=0; k<Ulinelen; k++) Out.put('-');
	Out << endl;

	/* print chunks for all rows: this varies for rectangular
	 * and triangular matrices so a virt function call will do it
	 */
	print_rows(Out, Sizew, Jbeg, Items, Width, Prec);

	/* chunk separator */
	for (k=0; k<Ulinelen; k++) Out.put('=');
	Out << endl << endl;
    }
    
    /* clean up mess */
    Out.width(0);	    // each item uses its length as width
    Out.precision(Oldprec);	    // old precision
    Out.flags(Oldflags);    // reset flags to state before call
}
// END of list_matrix() 

// ---- Errors ----

/* prt_err(): prints an error message to cerr. Etyp is the
 * type of the error, Funcnm is a string that contains the
 * name of the function in which the error has occurred.
 */
void Matbase_::prt_err(const Errtype_ Etyp, const char *Funcnm) const
{
    cerr << "? Matrix error in " << Funcnm <<": ";
    switch(Etyp)
    {
        case NO_MEM:
        cerr << "Out of memory\n"; break;
        case DIM_MISMATCH:
        cerr << "Dimension mismatch\n"; break;
        case DIV_BY_ZERO:
        cerr << "Division by zero\n"; break;
        case BAD_ROWRANGE:
        cerr << "Row index out of range\n"; break;
        case BAD_COLRANGE:
        cerr << "Col index out of range\n"; break;
        default:
        cerr << "Error (no code available)\n"; break;
    }
}
// END of prt_err()

// ==== END OF METHODS Matbase_ ====

// ---- GLOBAL FUNCTIONS ----

/* Overloaded version of <<: always calls the proper virtual
 * function for matrix listing.
 */
ostream& operator << (ostream& Ostr, const Matbase_& Mat)
{
    Mat.list_matrix(Ostr); return(Ostr);
}
// END of <<

// ==== END OF FUNCTIONS Matbase.c++ ====

