#ifndef __CPLX_CLASS__
#define __CPLX_CLASS__

// ==== HEADER Cplx.h ====

/* Double-precision complex arithmetics class. */

// SGI C++ 4.0, IRIX 5.3, 25. Oct. 1995. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <iostream.h>

/* ==== CLASSES ==== */

#define UNITY Complex_(1.0)
#define IMAG Complex_(0.0, 1.0)

/* Class Complex_ : double-precision class with the
 * usual arithmetics overloaded + a few transcendental
 * functions. Error checking is sloppy.
 */
class Complex_
{
    // data
    private:
    
    double Re, Im;  // real and imaginary parts
    
    // methods
    public:
    
	// constructors
    Complex_(double Real=0.0, double Imag=0.0): Re(Real), Im(Imag) {}
    
	// member access
    friend double real(const Complex_& C) { return(C.Re); }
    friend double imag(const Complex_& C) { return(C.Im); }
    
	// equality and arithmetics
    friend int operator==(const Complex_& C1, const Complex_& C2);
    Complex_& operator-();   // unary minus
    friend Complex_ operator+(const Complex_& C1, const Complex_& C2);
    friend Complex_ operator-(const Complex_& C1, const Complex_& C2);
    friend Complex_ operator*(const Complex_& C1, const Complex_& C2);
    friend Complex_ operator/(const Complex_& C1, const Complex_& C2);
	// arithmetics done "in place"
    Complex_& operator+=(const Complex_& C);
    Complex_& operator-=(const Complex_& C);
    Complex_& operator*=(const Complex_& C);
    Complex_& operator/=(const Complex_& C);
    
	// conjugation, fabs (abs value), argument
    Complex_ conjug() const ;
    double fabs() const ;
    double argument() const ;

	// conversions
    operator double() const { return((Im==0.0)? Re: fabs()); }
    
	// powers, sqrt, exp, log
    friend  Complex_ pow(const Complex_& C, double Ex);
    friend  Complex_ pow(const Complex_& C, const Complex_& Ex);
    friend  Complex_ sqrt(const Complex_& C);
    friend  Complex_ exp(const Complex_& C);
    friend  Complex_ log(const Complex_& C);
    
	// trigonometry
    friend  Complex_ sin(const Complex_& C);
    friend  Complex_ cos(const Complex_& C);
    friend  Complex_ tan(const Complex_& C);
    
	// inverse trigs
    friend  Complex_ asin(const Complex_& C);
    friend  Complex_ acos(const Complex_& C);
    friend  Complex_ atan(const Complex_& C);
    
	// input-output (very simple)
    friend istream& operator>>(istream& In, Complex_& C);
    friend ostream& operator<<(ostream& Out, const Complex_& C);
    
};
// END OF CLASS Complex_

// ==== END OF HEADER Cplx.h ====
#endif
