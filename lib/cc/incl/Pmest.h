#ifndef __PMEST_HEADER__
#define __PMEST_HEADER__

// ==== HEADER Pmest.h ====

/* Parameter estimation routines: all based on the book:
 * Valko, P., Vajda,  S.: M\H{u}szaki-tudom\'{a}nyos feladatok
 * megold\'{a}sa szem\'{e}lyi sz\'{a}m\'{\i}t\'{o}g\'{e}ppel, 
 * Bp,  1986.
 */

// SGI C++ 4.0, IRIX 5.3, 11. May 1995. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>

// ---- INCLUDE FILES ----

#include "Vector.h"
#include "Matrix.h"
#include "Trimat.h"
#include "Vmutils.h"

// ---- DEFINITIONS ----

#define NLIN_STEPLIM (1e-5)
#define NLIN_SILENT 0
#define NLIN_TALK 1
#define NLIN_CHATTER 2

// ---- TYPEDEFS ----

/* Userfunct_: the type of the function that describes
 * the nonlinear vector:vector mapping for nonlin_reg().
 */
typedef int (*Userfunct_)(const Vector_&, const Vector_&, Vector_&);

/* Userfunct11_ : the type of the function that describes
 * the 1:1 (scalar:scalar) nonlinear mapping for nonlin11_reg().
 */
typedef double (*Userfunct11_)(double, const Vector_&);

// ---- PROTOTYPES ----

/* lin_reg(): multiple weighed linear regression with "ridge" regularisation.
 * Input parameters:-
 * Xmeas: observation matrix of independent variables (Nm x Nx)
 * Ymeas: vector of dependent variable data (Nm)
 * W: weight vector (Nm)
 * Output parameters:-
 * P: vector of estimated parameters (Nx)
 * Sdev: standard deviation of parameters (Nx)
 * Correl: parameter correlation matrix (Nx x Nx)
 * Tstat95: "Student"'s t-statistics at 95 % level
 * Return value: the residual deviation (Q).
 * Nm and Nx are the number of observations and number of variables, 
 * respectively. The input dimensions are checked and if a mismatch
 * is detected then the routine returns with -1.0 and a warning is printed.
 * The output variables' dimensions are adjusted if necessary.
 */
float lin_reg(const Matrix_& Xmeas, const Vector_& Ymeas, const Vector_& W, 
	    Vector_& P, Vector_& Sdev, Trimat_& Correl, float& Tstat95);

/* nonlin_reg(): estimates the parameters in the function Y=Funct(X,P).
 * X is an Nx-long vector, P is the Np-long parameter vector, Y is
 * an Ny-long vector of dependent variables and Nm observations are
 * supplied. 
 * Input parameters:-
 *   Funct: the nonlinear vector-vector mapping between X and Y
 *   Xmeas: Nm x Nx matrix of independent variables
 *   Ymeas: Nm x Ny matrix of dependent variables
 *   W: weight matrix (Nm x Ny)
 *   Steplim: iteration stepsize (default=NLIN_STEPLIM)
 *   Verbose: <=0: no output (default), 1: '....Done', >=2: detailed (itno, Q)
 * Output parameters:-
 *   Sdev: SD vector of parameters (Np)
 *   Correl: correlation matrix (Np x Np)
 *   Tcrit95: "Student"'s t-value at 95 % confidence level
 * Input/output parameters:-
 *   P: initial parameter guess on input, estimated parameters on output
 *   Itmax: max. allowed no. of iterations on input, actual on output
 * Return value:- the residual deviation or -1.0 on error.
 * NOTE: Funct should have the following prototype (typedef Userfunct_):-
 * int Funct(const Vector_& X, const Vector_& P, Vector_& Y)
 * and should return 0 on error, non-zero otherwise. Before the iteration
 * starts, a call is made to Funct() to check possible dim mismatches
 * and if 0 is returned then nonlin_reg() returns with an error.
 * It is the caller's responsibility to ensure that Funct() carries
 * out error checking properly (dim + math errors).
 * During iteration, the original signs of the initial parameter guesses
 * are kept. A parameter set to 0.0 will remain zero. Try sign changes
 * if the fit is bad.
 */
float nonlin_reg(const Matrix_& Xmeas, const Matrix_& Ymeas, const Matrix_& W, 
    Userfunct_ Funct, 
    Vector_& P, Vector_& Sdev, Trimat_& Correl, float& Tcrit95, 
    int& Itmax, float Steplim=NLIN_STEPLIM, int Verbose=NLIN_SILENT);

/* nonlin11_reg(): estimates the parameters in the function Y=Funct(X,P)
 * where X and Y are scalars. P is the Np-long parameter vector.
 * Otherwise the algorithm is the same as in nonlin_reg() above.
 * Input parameters:-
 *   Funct: the nonlinear scalar-scalar mapping between X and Y
 *   Xmeas: Nm vector of independent variables
 *   Ymeas: Nm vector of dependent variables
 *   W: weight vector
 *   Steplim: iteration stepsize (default=NLIN_STEPLIM)
 *   Verbose: <=0: no output (default), 1: '....Done', >=2: detailed (itno, Q)
 * Output parameters:-
 *   Sdev: SD vector of parameters (Np)
 *   Correl: correlation matrix (Np x Np)
 *   Tcrit95: "Student"'s t-value at 95 % confidence level
 * Input/output parameters:-
 *   P: initial parameter guess on input, estimated parameters on output
 *   Itmax: max. allowed no. of iterations on input, actual on output
 * Return value:- the residual deviation or -1.0 on error.
 * NOTE: Funct should have the following prototype (typedef Userfunct11_) :-
 * double Funct(double X, const Vector_& P)
 * It is the caller's responsibility to ensure that Funct() carries
 * out error checking properly (dim + math errors).
 * During iteration, the original signs of the initial parameter guesses
 * are kept. A parameter set to 0.0 will remain zero. Try sign changes
 * if the fit is bad.
 */
float nonlin11_reg(const Vector_& Xmeas, const Vector_& Ymeas, const Vector_& W, 
    Userfunct11_ Funct, 
    Vector_& P, Vector_& Sdev, Trimat_& Correl, float& Tcrit95, 
    int& Itmax, float Steplim, int Verbose);

// ==== END OF HEADER Pmest.h ====

#endif



