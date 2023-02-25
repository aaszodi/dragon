// ==== FUNCTIONS Cplx.c++ ====

/* Double-precision complex arithmetics class. */

// SGI C++ 4.0, IRIX 5.3, 25. Oct. 1995. Andris */

/* ---- HEADER FILE ---- */

#include "Cplx.h"

/* ==== FUNCTIONS ==== */

    // equality and arithmetics
int operator==(const Complex_& C1, const Complex_& C2)
    { return(C1.Re==C2.Re && C1.Im==C2.Im); }
Complex_& Complex_::operator-()   // unary minus
    { Re*=-1.0; Im*=-1.0; return(*this); }
Complex_ operator+(const Complex_& C1, const Complex_& C2)
    { return Complex_(C1.Re+C2.Re, C1.Im+C2.Im); }
Complex_ operator-(const Complex_& C1, const Complex_& C2)
    { return Complex_(C1.Re-C2.Re, C1.Im-C2.Im); }
Complex_ operator*(const Complex_& C1, const Complex_& C2)
    { return Complex_(C1.Re*C2.Re-C1.Im*C2.Im, C1.Re*C2.Im+C1.Im*C2.Re); }
Complex_ operator/(const Complex_& C1, const Complex_& C2)
{
    register double Abs2=C2.Re*C2.Re+C2.Im*C2.Im;
    
    if (Abs2<=DBL_MIN)
    {
	cerr << "? Complex_:: division by zero, (0, 0) returned\n";
	return(Complex_(0.0, 0.0));
    }
    else return(Complex_((C1.Re*C2.Re+C1.Im*C2.Im)/Abs2, (C1.Im*C2.Re-C1.Re*C2.Im)/Abs2));
}
    // arithmetics done "in place"
Complex_& Complex_::operator+=(const Complex_& C)
    { Re+=C.Re; Im+=C.Im; return(*this); }
Complex_& Complex_::operator-=(const Complex_& C)
    { Re-=C.Re; Im-=C.Im; return(*this); }
Complex_& Complex_::operator*=(const Complex_& C)
    { *this=*this * C; return(*this); }
Complex_& Complex_::operator/=(const Complex_& C)
    { *this=*this / C; return(*this); }

    // conjugation, fabs, argument
Complex_ Complex_::conjug() const { return(Complex_(Re, -Im)); }
double Complex_::fabs() const { return(sqrt(Re*Re+Im*Im)); }
double Complex_::argument() const { return(atan2(Im, Re)); }

/* transcendental functions: these are not members but rather
 * friends of Complex_ since there is no need for calls like
 * "C.exp()";  the "exp(C)" syntax is much better
 */
    // powers, sqrt, exp, log
Complex_ pow(const Complex_& C, double Ex)
{
    register double Mod=C.fabs();
    register double Arg=C.argument();
    
    Mod=pow(Mod, Ex); Arg*=Ex;
    return(Complex_(Mod*cos(Arg), Mod*sin(Arg)));
}
Complex_ sqrt(const Complex_& C)
{
    register double Mod=sqrt(C.fabs());
    register double Arg=C.argument()/2.0;
    return(Complex_(Mod*cos(Arg), Mod*sin(Arg)));
}
Complex_ exp(const Complex_& C)
{
    register double Ex=exp(C.Re);
    return(Complex_(Ex*cos(C.Im), Ex*sin(C.Im)));
}
Complex_ log(const Complex_& C)
    { return(Complex_(log(C.fabs()), C.argument())); }
Complex_ pow(const Complex_& C, const Complex_& Ex)
    { return(exp(Ex*log(C))); }

    // trigonometry
Complex_ sin(const Complex_& C)
    { return(Complex_(sin(C.Re)*cosh(C.Im), cos(C.Re)*sinh(C.Im))); }
Complex_ cos(const Complex_& C)
    { return(Complex_(cos(C.Re)*cosh(C.Im), -sin(C.Re)*sinh(C.Im))); }
Complex_ tan(const Complex_& C)
    { return(sin(C)/cos(C)); }

    // inverse trigs
Complex_ asin(const Complex_& C)
    { return(-IMAG*log(IMAG*C+sqrt(UNITY-C*C))); }
Complex_ acos(const Complex_& C)
    { return(IMAG*log(C+sqrt(C*C-UNITY))); }
Complex_ atan(const Complex_& C)
    { return(Complex_(0.5)*log((UNITY+IMAG*C)/(UNITY-IMAG*C))/IMAG); }

    // input-output (very simple)
istream& operator>>(istream& In, Complex_& C)
    { In >> C.Re >> C.Im; return(In); }
ostream& operator<<(ostream& Out, const Complex_& C)
    { Out << '<' << C.Re << ',' << C.Im << '>'; return(Out); }

// ==== END OF FUNCTIONS Cplx.c++ ====
