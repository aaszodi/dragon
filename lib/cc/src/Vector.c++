// ==== FUNCTIONS Vector.c++ ====

/* Double-precision vector class for simple vector algebra.
 */

// SGI C++ 4.0, IRIX 5.3, 28. Sept. May 1995. Andris

// ---- HEADER ----

#include "Vector.h"

// ==== Vector_ DEFINITIONS ==== 

// ---- Constructors ----

/* init to N-dimensional null-vector.
 * Default dim is 3, 0 dim not permitted
 */
Vector_::Vector_(unsigned int N)
{
    if (!N) N=1;
    X=new double [Dim=N]; 
    if (X==NULL) prt_err(NO_MEM, "Vector(N)");	// just screaming
    else memset(X, 0, Dim*sizeof(double));
}

/* init to values contained in the array Arr[] which is assumed
 * to be N items long. If Arr==NULL or N==0, then a N- or 1-dimensional
 * null-vector is created, respectively.
 */
Vector_::Vector_(const double Arr[], unsigned int N)
{
	// handle pathological cases here
    if (!N) N=1; 
    X=new double [Dim=N];
    if (X==NULL) prt_err(NO_MEM, "Vector(Arr, N)");
    else if (Arr==NULL)	// N-dimensional null-vector
	memset(X, 0, Dim*sizeof(double));
    else
	memcpy(X, Arr, Dim*sizeof(double));
}

/* Init by another Vector_ object. */
Vector_::Vector_(const Vector_& Vec)
{
    X=new double[Dim=Vec.Dim]; 
    if (X==NULL) prt_err(NO_MEM, "Vector(Vec)");	// just screaming
    else memcpy(X, Vec.X, Dim*sizeof(double));
}

// ---- Access ----

/* The operator [] will provide "unsafe" access, the () operator
 * will give "safe" (index range checked) access. Bad ranges will 
 * cause a warning to be printed and the index reset to 0.
 * Note that there are const and non-const versions of these functions.
 */

const double& Vector_::operator()(unsigned int Idx) const
{
    if (Idx>=Dim) { prt_err(BAD_IDX, "()"); Idx=0; }
    return(X[Idx]);
}

double& Vector_::operator()(unsigned int Idx)
{
    if (Idx>=Dim) { prt_err(BAD_IDX, "()"); Idx=0; }
    return(X[Idx]);
}

/* get_array(): constructs and returns a conventional array
 * which contains the coordinates. Returns the dimension in Len.
 * If allocation fails, NULL is returned (Len will not be changed).
 */
double* Vector_::get_array(unsigned int& Len) const
{
    double *Arr=new double [Dim];
    if (Arr==NULL) { prt_err(NO_MEM, "get_array()"); return(NULL); }
    memcpy(Arr, X, Dim*sizeof(double)); Len=Dim;
    return(Arr);
}
// END of get_array()

/* set_values(): sets all coordinates to Val (default=0.0) */
Vector_& Vector_::set_values(double Val) 
{
    register unsigned int i;
    for (i=0; i<Dim; i++) X[i]=Val;
    return(*this);
}
// END of set_values()

// ---- Dimension change ----

/* dim(): with no parameters, returns the current dimension.
 * With the parameter N, sets the dimension to N. If N==Dim or there is
 * no memory available, no action is taken. if N<Dim, then
 * only the first N coordinates are kept. If N>Dim, then the
 * old Dim coordinates will be padded up by 0.0.
 */
void Vector_::dim(unsigned int N)
{
    if (N==Dim) return;	// do nothing
    
    double *New=new double [N];	// allocate storage for new items
    if (New==NULL) { prt_err(NO_MEM, "set_dim()"); return; }
    
    memcpy(New, X, ((N>Dim)? Dim: N)*sizeof(double));	// copy old
    delete [] X; X=New; 	    // replace old data members with new
    if (N>Dim) memset(X+Dim, 0, (N-Dim)*sizeof(double)); // pad with zeroes
    Dim=N;	    // now you may adjust the dimension 
    return;
}
// END of dim(N)

// ---- Arithmetics ----

/* Assignment, vector addition and subtraction, postfix
 * multiplication and division by scalar. All operations
 * are available in in-place versions as well. The assignment
 * is destructive, it resets the dimension of the target.
 * In case of a dim mismatch, the += and -= operators 
 * do not modify the target, the + and - operators return
 * the left operand.
 */

Vector_& Vector_::operator=(const Vector_& Vec)
{
    if (this==&Vec) return(*this);  // x=x
    
    double *New=new double [Vec.Dim];
    if (New==NULL) { prt_err(NO_MEM, "="); return(*this); }
    Dim=Vec.Dim;
    delete [] X; X=New; memcpy(X, Vec.X, Dim*sizeof(double));
    return(*this);
}
// END of operator=

Vector_ Vector_::operator+(const Vector_& Vec) const
{
    if (Dim!=Vec.Dim) { prt_err(DIM_MISMATCH, "Vec+Vec"); return(*this); }
    Vector_ Temp(*this);
    Temp+=Vec; return(Temp);
}

Vector_& Vector_::operator+=(const Vector_& Vec)
{
    if (Dim!=Vec.Dim) { prt_err(DIM_MISMATCH, "Vec+=Vec"); return(*this); }
    unsigned int i;
    for (i=0; i<Dim; i++) X[i]+=Vec[i];
    return(*this);
}

Vector_ Vector_::operator-(const Vector_& Vec) const
{
    if (Dim!=Vec.Dim) { prt_err(DIM_MISMATCH, "Vec-Vec"); return(*this); }
    Vector_ Temp(*this);
    Temp-=Vec; return(Temp);
}

Vector_& Vector_::operator-=(const Vector_& Vec)
{
    if (Dim!=Vec.Dim) { prt_err(DIM_MISMATCH, "Vec-=Vec"); return(*this); }
    unsigned int i;
    for (i=0; i<Dim; i++) X[i]-=Vec[i];
    return(*this);
}

Vector_ Vector_::operator*(double Scal) const
{
    Vector_ Temp(*this);
    Temp*=Scal; return(Temp);
}

Vector_ operator*(double Scal, const Vector_& Vec)
{
    Vector_ Temp=Vec;
    Temp*=Scal; return(Temp);
}

Vector_& Vector_::operator*=(double Scal)
{
    register unsigned int i;
    for (i=0; i<Dim; i++) X[i]*=Scal;
    return(*this);
}

Vector_ Vector_::operator/(double Scal) const
{
    Vector_ Temp(*this);
    if (fabs(Scal)<DBL_EPSILON) prt_err(DIV_BY_ZERO, "Vec/Scal");
    else Temp/=Scal;
    return(Temp);
}

Vector_& Vector_::operator/=(double Scal)
{
    if (fabs(Scal)<DBL_EPSILON) prt_err(DIV_BY_ZERO, "Vec/=Scal");
    else
    {	// "multiply-by-reciprocal" trick
	Scal=1.0/Scal;
	for (register unsigned int i=0; i<Dim; i++) X[i]*=Scal;
    }
    return(*this);
}

// ---- Vectorial products ----

/* the scalar product: if the dimensions of the vectors differ,
 * then the lower dimension is used (as if the lower-dim
 * vector were padded up w/ zeroes).
 */
double Vector_::operator*(const Vector_& Vec) const
{
    register double Prd=0.0;
    register unsigned int i, D=(Dim<Vec.Dim)? Dim: Vec.Dim;
    for (i=0; i<D; i++) Prd+=X[i]*Vec[i];
    return(Prd);
}
// END of operator* (scalar product)

/* cross_prod(): performs the cross-product of two 3D vectors.
 * If the arguments are not 3D, then a 3D null-vector is returned.
 */
Vector_ cross_prod(const Vector_& Vec1, const Vector_& Vec2)
{
    Vector_ Temp(3);
    if (Vec1.Dim==3 && Vec2.Dim==3)
    {
	Temp.X[0]=Vec1.X[1]*Vec2.X[2]-Vec1.X[2]*Vec2.X[1];
	Temp.X[1]=Vec1.X[2]*Vec2.X[0]-Vec1.X[0]*Vec2.X[2];
	Temp.X[2]=Vec1.X[0]*Vec2.X[1]-Vec1.X[1]*Vec2.X[0];
    }
    else Temp.prt_err(Vector_::DIM_MISMATCH, "Vec x Vec"); 
    return(Temp);
}
// END of cross_prod()

// ---- Modulus ----

/* vec_len2(): calculates the square of the Euclidean norm. */
double Vector_::vec_len2() const
{
    register double Len=0.0;
    register unsigned int i;
    for (i=0; i<Dim; i++) Len+=X[i]*X[i];
    return(Len);
}
// END of vec_len2()

/* vec_norm(): normalises the calling object to a unit vector.
 * Returns the original length. If the length<DBL_EPSILON
 * then the calling object will be set to an exact null-vector
 * and 0.0 will be returned.
 */
double Vector_::vec_norm()
{
    register double L=vec_len();
    if (L<DBL_EPSILON)
    { L=0.0; set_values(); }	// exact null-vector
    else (*this)/=L;		// normalise
    return(L);
}
// END of vec_norm()

/* diff_len2(): returns the squared length of Vec1-Vec2. On dim
 * mismatch a warning is printed and 0.0 returned.
 */
double diff_len2(const Vector_& Vec1, const Vector_& Vec2)
{
    if (Vec1.Dim!=Vec2.Dim)
    {
	Vec1.prt_err(Vector_::DIM_MISMATCH, "|Vec1-Vec2|^2"); return(0.0);
    }
    register unsigned int i;
    register double D, L=0.0;
    for (i=0; i<Vec1.Dim; i++)
    {
	D=Vec1.X[i]-Vec2.X[i];
	L+=D*D;
    }
    return(L);
}
// END of diff_len2()

// ---- Printing ----

/* Overloaded version of <<: note that this is a GLOBAL
 * function. Calls list_vector() with the default arguments.
 */
ostream& operator << (ostream& Ostr, const Vector_& Vec)
{
    Vec.list_vector(Ostr); return(Ostr);
}
// END of <<

/* list_vector(): lists the calling object to stream Out
 * in a neat column format using scientific notation with
 * Prec digits precision (default=2). This function is called by the
 * overloaded << operator with the default arguments.
 */
void Vector_::list_vector(ostream& Out, unsigned int Prec) const
{
    char Nwidth[10];  // buffer to determine the printed width of N
    ostrstream Ostr(Nwidth, sizeof(Nwidth));	// stream output buffer
    long Oldflags=Out.flags();	// store previous settings here
    int Oldprec=Out.precision(Prec);
    unsigned int i, k, N, Width, Sizew, Ulinelen;

    N=dim();
    Ostr << N << ends;	    	// get printed width of indices
    Sizew=Ostr.pcount()-1;
    Width=Sizew+Prec+5;		// coordinate field width
    Ulinelen=Width+Sizew+4;	// underlining length (decoration)

    // underline column heading
    for (k=0; k<Ulinelen; k++) Out.put('-');
    Out << endl;
    
    // output will be right-justified by default
    Out.unsetf(ios::right|ios::internal);
    Out.setf(ios::showpoint|ios::scientific);	// %e in C
    
    // print the values in rows
    for (i=0; i<N; i++)
    {
	Out.unsetf(ios::left);
	Out << setw(Sizew) <<i<< setw(0) <<" | " ;   // row idx
	Out.setf(ios::left);
	Out<<setw(Width)<<X[i]<<endl;
    }
    
    // double underline at the end
    for (k=0; k<Ulinelen; k++) Out.put('=');
    Out << endl << endl;
    
    // clean up mess
    Out.width(0);	    // each item uses its length as width
    Out.precision(Oldprec);	    // old value
    Out.flags(Oldflags);    // reset flags to orig state
}
// END of list_vector()
    
// ---- Errors ----

/* prt_err(): prints an error message to cerr. Etyp is the
 * type of the error, Funcnm is a string that contains the
 * name of the function in which the error has occurred.
 */
void Vector_::prt_err(const Errtype_ Etyp, const char *Funcnm) const
{
    cerr << "? Vector_::" << Funcnm <<": ";
    switch(Etyp)
    {
        case NO_MEM:
        cerr << "Out of memory\n"; break;
        case DIM_MISMATCH:
        cerr << "Dimension mismatch\n"; break;
        case DIV_BY_ZERO:
        cerr << "Division by zero\n"; break;
        case BAD_IDX:
        cerr << "Index out of range\n"; break;
        default:
        cerr << "Error (no code available)\n"; break;
    }
}
// END of prt_err()

// ==== END OF FUNCTIONS Vector.c++ ====
