#ifndef __VECTOR_CLASS__
#define __VECTOR_CLASS__

// ==== HEADER Vector.h ====

/* Double-precision vector class for simple vector algebra.
 */

// SGI C++ 4.0, IRIX 5.3, 18. May 1995. Andris

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <iostream.h>
#include <strstream.h>
#include <iomanip.h>

// ==== CLASSES ====

/* Class Vector_ : double-precision vector algebra class.
 * The coordinates are kept in a simple double array.
 * Vectors can be converted to conventional arrays.
 */
class Vector_
{
    // data
    double *X;  // coordinate array
    unsigned int Dim;   // dimensions
    
    // functions
    public:
        // constructors
    /* init to N-dimensional null-vector.
     * Default dim is 3, 0 dim not permitted
     */
    Vector_(unsigned int N=3);
    
    /* init to values contained in the array Arr[] which is assumed
     * to be N items long. If Arr==NULL or N==0, then a N- or 1-dimensional
     * null-vector is created, respectively.
     */
    Vector_(const double Arr[], unsigned int N);
    
    /* Init by another Vector_ object. */
    Vector_(const Vector_& Vec);
    
        // destructor
    ~Vector_() { delete [] X; }
    
        // access
    /* The operator [] will provide "unsafe" access, the () operator
     * will give "safe" (index range checked) access. Bad ranges will 
     * cause a warning to be printed and the index reset to 0.
     * Note that there are const and non-const versions of these functions.
     */

    double& operator[](unsigned int Idx) { return(X[Idx]); }
    const double& operator[](unsigned int Idx) const { return(X[Idx]); };
    const double& operator()(unsigned int Idx) const;
    double& operator()(unsigned int Idx);
    
    /* get_array(): constructs and returns a conventional array
     * which contains the coordinates. Returns the dimension in Len.
     * If allocation fails, NULL is returned (Len will not be changed).
     */
    double* get_array(unsigned int& Len) const;
    
    /* set_values(): sets all coordinates to Val (default=0.0) */
    Vector_& set_values(double Val=0.0);
    
        // dimension access
    /* dim(): with no parameters, returns the current dimension.
     * With the parameter N, sets the dimension to N. If N==Dim or there is
     * no memory available, no action is taken. if N<Dim, then
     * only the first N coordinates are kept. If N>Dim, then the
     * old Dim coordinates will be padded up by 0.0.
     */
    unsigned int dim() const { return(Dim); }
    void dim(unsigned int N);
    
        // arithmetics
    /* Assignment, vector addition and subtraction, postfix
     * multiplication and division by scalar. All operations
     * are available in in-place versions as well. The assignment
     * is destructive, it resets the dimension of the target.
     * In case of a dim mismatch, the += and -= operators 
     * do not modify the target, the + and - operators return
     * the left operand.
     */
    Vector_& operator=(const Vector_& Vec);
    Vector_ operator+(const Vector_& Vec) const;
    Vector_& operator+=(const Vector_& Vec);
    Vector_ operator-(const Vector_& Vec) const;
    Vector_& operator-=(const Vector_& Vec);
    Vector_ operator*(double Scal) const;
    friend Vector_ operator*(double Scal, const Vector_& Vec);
    Vector_& operator*=(double Scal);
    Vector_ operator/(double Scal) const;
    Vector_& operator/=(double Scal);
    
        // vectorial products
    /* the scalar product: if the dimensions of the vectors differ,
     * then the lower dimension is used (as if the lower-dim
     * vector were padded up w/ zeroes).
     */
    double operator*(const Vector_& Vec) const;
    
    /* cross_prod(): performs the cross-product of two 3D vectors.
     * If the arguments are not 3D, then a 3D null-vector is returned.
     */
    friend Vector_ cross_prod(const Vector_& Vec1, const Vector_& Vec2);
    
        // modulus
    /* vec_len(), vec_len2(): calculate the Euclidean norm and the
     * squared Euclidean norm,  respectively.
     */
    double vec_len() const { return(sqrt(vec_len2())); }
    double vec_len2() const;
    
    /* vec_norm(): normalises the calling object to a unit vector.
     * Returns the original length. If the length<DBL_EPSILON
     * then the calling object is considered a null-vector, will be set to
     * exact length 0.0 and 0.0 will be returned.
     */
    double vec_norm();
    
    /* diff_len(),diff_len2(): return the length of Vec1-Vec2 and
     * the squared length,  respectively. On dim
     * mismatch a warning is printed and 0.0 returned.
     */
    friend double diff_len(const Vector_& Vec1, const Vector_& Vec2)
	{ return(sqrt(diff_len2(Vec1, Vec2))); }
    friend double diff_len2(const Vector_& Vec1, const Vector_& Vec2);
    
        // printing
    /* list_vector(): lists the calling object to stream Out
     * in a neat column format using scientific notation with
     * Prec digits precision (default=2). This function is called by the
     * overloaded << operator with the default arguments.
     */
    void list_vector(ostream& Ostr, unsigned int Prec=2) const;
    
    private:
        // errors
    enum Errtype_ { NO_MEM, DIV_BY_ZERO, BAD_IDX, DIM_MISMATCH };
        // error messages
    /* prt_err(): prints an error message to cerr. Etyp is the
     * type of the error, Funcnm is a string that contains the
     * name of the function in which the error has occurred.
     */
    void prt_err(const Errtype_ Etyp, const char *Funcnm) const;
};
// END OF CLASS Vector_

/* Overloaded version of <<: note that this is a GLOBAL
 * function. Calls list_vector() with the default arguments.
 */
ostream& operator << (ostream& Out, const Vector_& Vec);

// ==== END OF HEADER Vector.h ====

#endif
