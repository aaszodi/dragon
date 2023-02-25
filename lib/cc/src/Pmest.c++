// ==== FUNCTIONS Pmest.c++ ====

/* Parameter estimation routines: all based on the book:
 * Valko, P., Vajda,  S.: M\H{u}szaki-tudom\'{a}nyos feladatok
 * megold\'{a}sa szem\'{e}lyi sz\'{a}m\'{\i}t\'{o}g\'{e}ppel, 
 * Bp,  1986.
 */

// SGI C++ 4.0, IRIX 5.3, 11. May 1995. Andris Aszodi

// ---- HEADER ----

#include "Pmest.h"

// ---- DEFINITIONS ----

#ifndef DBL_MIN
#define DBL_MIN (1e-100)
#endif

// Some math libraries lack the float equivalents of the
// usual double functions (eg. float sqrtf(float) in addition
// to double sqrt(double). Define NO_MATHFLOATFUNC in these cases.
#ifdef NO_MATHFLOATFUNC
#define fabsf fabs
#endif

// ---- PROTOTYPES ----

static float tcrit_95(long Nf);
static int posdef_inv(Trimat_& A);

// ==== FUNCTIONS ====

// ---- Linear regression ----

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
	    Vector_& P, Vector_& Sdev, Trimat_& Correl, float& Tstat95)
{
    // input dim checks (Xmeas dimensions are the "standard")
    unsigned int Nm=Xmeas.rno(), Nx=Xmeas.cno();
    if (Nm<=Nx)
    {
	cerr<<"\n! lin_reg(): Please supply at least "<<Nx+1
	    <<" measurements instead of "<<Nm<<endl;
	return(-1.0);
    }
    
    
    if (Ymeas.dim()!=Nm)
    {
	cerr<<"\n! lin_reg(): Ymeas.dim()="<<Ymeas.dim()<<" instead of "<<Nm<<endl;
	return(-1.0);
    }
    
    if (W.dim()!=Nm)
    {
	cerr<<"\n! lin_reg(): W.dim()="<<W.dim()<<" instead of "<<Nm<<endl;
	return(-1.0);
    }
    
    // output dimension checks
    if (P.dim()!=Nx)
    {
	cerr<<"\n? lin_reg(): P.dim()="<<P.dim()<<" instead of "<<Nx<<", adjusted\n";
	P.dim(Nx);
    }
    
    if (Sdev.dim()!=Nx)
    {
	cerr<<"\n? lin_reg(): Sdev.dim()="<<Sdev.dim()<<" instead of "<<Nx<<", adjusted\n";
	Sdev.dim(Nx);
    }
    
    if (Correl.rno()!=Nx)
    {
	cerr<<"\n? lin_reg(): Correl.rno()="<<Correl.rno()<<" instead of "<<Nx<<", adjusted\n";
	Correl.set_size(Nx);
    }
    
    // get the X'WX matrix and the X'WY vector
    Trimat_ Xtx=trans_mwprod(Xmeas, W);

    double Temp;
    unsigned int im, ip, jp;
    Vector_ Xty(Nx);
    
    for (ip=0; ip<Nx; ip++)
    {
	Temp=0.0;
	for (im=0; im<Nm; im++)
	    Temp+=Xmeas[im][ip]*W[im]*Ymeas[im];
	Xty[ip]=Temp;
    }
    
    // normalise X'WX
    double Trace=Xtx.get_trace();
    Trace/=(Nx*1000.0);	    // no norm with pivots smaller than this
    Vector_ Norm=Xtx.diag();
    for (ip=0; ip<Nx; ip++)
	Norm[ip]=(Norm[ip]>Trace)? sqrt(Norm[ip]): 1.0;

    for (ip=0; ip<Nx; ip++)
	for (jp=0; jp<=ip; jp++)
	    Xtx[ip][jp]/=(Norm[ip]*Norm[jp]);
    Trimat_ Xtxold(Xtx);    // preserve normalised Xtx
    
    // inversion (Rid is increased if Xtx is singular)
    double Rid=0.0;
    while(!posdef_inv(Xtx))
    {
	cerr<<"\n? lin_reg(): Xtx is singular (Rid="<<Rid<<")\n";
	Xtx=Xtxold;	// get back old Xtx
	Rid+=0.01;	// increment ridge parameter
	for (ip=0; ip<Nx; ip++)
	    Xtx[ip][ip]+=Rid;	// and add to diagonal
    }

    // solution
    for (ip=0; ip<Nx; ip++)	// norm back: divide again
	for (jp=0; jp<=ip; jp++)
	    Xtx[ip][jp]/=(Norm[ip]*Norm[jp]);

    P=Xtx*Xty;

    // quality function
    float D, Q=0.0;
    Vector_ Yest=Xmeas*P;   // estimated Y values

    for (im=0; im<Nm; im++)
    {
	D=Yest[im]-Ymeas[im];
	Q+=W[im]*D*D;
    }
    
    // rest of output
    int Nf=Nm-Nx;
    Q=sqrt(Q/Nf);   // residual deviation
    Tstat95=tcrit_95(Nf);   // t-statistics
    
    Correl=Xtx*(Q/Nf);	// covariance matrix
    Sdev=Correl.diag();	// variances
    for (ip=0; ip<Nx; ip++) Sdev[ip]=sqrt(Sdev[ip]);
    
    double Sij;
    for (ip=0; ip<Nx; ip++)	// norm to correlation
    {
	Correl[ip][ip]=1.0;
	for (jp=0; jp<ip; jp++)
	{
	    Sij=Sdev[ip]*Sdev[jp];
	    if (Sij<DBL_EPSILON) Correl[ip][jp]=0.0;
	    else Correl[ip][jp]/=Sij;
	}
    }

    return(Q);
}
// END of lin_reg()

// ---- Nonlinear regression ----

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
 * NOTE: Funct should have the following prototype:-
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
    int& Itmax, float Steplim, int Verbose) 
{
    // dim checks
    int Nm=Xmeas.rno(), Nx=Xmeas.cno(), Ny=Ymeas.cno(), Np=P.dim(), Nf;
    if (Nm!=Ymeas.rno())
    {
	cerr<<"\n! nonlin_reg(): "<<Ymeas.rno()<<" measurements in Ymeas instead of "<<Nm<<endl;
	return(-1.0);
    }
    if (Nm!=W.rno() || Ny!=W.cno())
    {
	cerr<<"\n! nonlin_reg(): W is ("<<W.rno()<<" x "<<W.cno()
	    <<") instead of ("<<Nm<<" x "<<")\n";
	return(-1.0);
    }
    if ((Nf=Nm*Ny-Np)<=0)
    {
	cerr<<"\n! nonlin_reg(): Not enough degrees of freedom ("<<Nf<<")\n";
	return(-1.0);
    }
    
    // less serious trouble (adjust silently)
    if (Np!=Sdev.dim()) Sdev.dim(Np);
    if (Np!=Correl.rno()) Correl.set_size(Np);
    Itmax=abs(Itmax); if (!Itmax) Itmax=100;
    Steplim=fabsf(Steplim);
    
    // adjust verbosity: used to be a param but caused trouble
    if (Verbose<=NLIN_SILENT) Verbose=NLIN_SILENT;	// no chat
    if (Verbose>=NLIN_CHATTER) Verbose=NLIN_CHATTER;	// lots of output

    // some of the local vars
    register int im, iy, ip, jp;
    double Q, Qold;
    
    // get the initial quality 
    Vector_ F(Ny), Dly(Ny);
    Qold=0.0;
    for (im=0; im<Nm; im++)
    {
	if (!Funct(Xmeas.row(im), P, F))
	{
	    cerr<<"\n! nonlin_reg(): Funct() returns error (dim mismatch?)\n";
	    return(-1.0);
	}
	Dly=Ymeas.row(im)-F;    // difference
	for (iy=0; iy<Ny; iy++)
	    Qold+=Dly[iy]*W[im][iy]*Dly[iy];
    }
    if (Verbose==NLIN_TALK) cout<<"\nnonlin_reg():";
    if (Verbose==NLIN_CHATTER) cout<<"\nNonlinear regression:\nItno\tQ\tLm\n";

    // set up matrices and vectors
    Matrix_ Jac(Ny, Np);    // Jacobian: dF/dP
    Vector_ Dy(Ny), Fd(Ny), Pnew(Np), Pd(Np), JtDy(Np), Db(Np), Norm(Np);
    Trimat_ Jtj(Np), Jtjold(Np);
    int Itno=0;
    double Lm=0.01, Stlen, Stfac, Dp, Trace;
    char Grow;
    
    const float DERIV_COEFF=0.001, DERIV_MINSTEP=1E-6;
    
    do	// main iteration cycle
    {
	Jtj.set_values(); JtDy.set_values();	// zeroing
	
	// calculate Jtj=J'WJ and JtDy=J'W(Y-F(X,P)) from observed data
	for (im=0; im<Nm; im++)	    // all measurement points
	{
	    if (!Funct(Xmeas.row(im), P, F))
	    {
		cerr<<"\n? nonlin_reg(): Error in Funct() evaluation\n";
		break;
	    }
	    Dy=Ymeas.row(im)-F; // obs. vs. calc.
	    
	    // numerical derivation for the Jacobian
	    for (ip=0; ip<Np; ip++)
	    {
		Dp = DERIV_COEFF * fabs(P[ip]) + DERIV_MINSTEP;
		Pd=P; Pd[ip]+=Dp;   // shift ip-th parameter
		Funct(Xmeas.row(im), Pd, Fd);	// funct values at shifted pos
		Jac.col((Fd-F)*(P[ip]/Dp), ip);	// normalised Jacobian column
	    }
	    
	    // now J'WJ and the "right-hand side" JtDy
	    for (iy = 0; iy < Ny; iy++)
		for (ip = 0; ip < Np; ip++)
		{
		    for (jp = 0; jp <= ip; jp++)
			Jtj[ip][jp] += Jac[iy][ip] * W[im][iy] * Jac[iy][jp];
		    JtDy[ip] += Jac[iy][ip] * W[im][iy] * Dy[iy];
		}
	}	// for im
	
	// normalise J'WJ
	Trace=Jtj.get_trace();
	Trace/=(Np*1000);   // no norm with pivots smaller than this
	Norm=Jtj.diag();    // copy diagonal
	for (ip=0; ip<Np; ip++)
	    Norm[ip]=(Norm[ip]>Trace)? sqrt(Norm[ip]): 1.0;
	
	for (ip=0; ip<Np; ip++)
	    for (jp=0; jp<=ip; jp++)
		Jtj[ip][jp]/=(Norm[ip]*Norm[jp]);
	Jtjold=Jtj;	// preserve original 
	
	do  // a linreg step, repeat w/ higher Marquardt Lm if Q grows
	{
	    // apply current Marquardt-lambda
	    for (ip=0; ip<Np; ip++) Jtj[ip][ip]+=Lm;
	    
	    // invert Jtj
	    while (!posdef_inv(Jtj))
	    {
		Jtj=Jtjold;
		Lm = (Lm==0.0)? 0.01: 10.0*Lm;
		for (ip=0; ip<Np; ip++) Jtj[ip][ip]+=Lm;
	    }
	    
	    // norm inverse back
	    for (ip=0; ip<Np; ip++)
		for (jp=0; jp<=ip; jp++)
		    Jtj[ip][jp]/=(Norm[ip]*Norm[jp]);
	    Db=Jtj*JtDy;    // linreg step solution
	    
	    // check stepsize
	    Stlen=Db.vec_len();	// overall length
	    Stfac=1.0;
	    for (ip=0; ip<Np; ip++)
		if (Stfac*Db[ip]<-0.95) Stfac=-0.95/Db[ip];
	    Db*=Stfac; Stlen*=Stfac;
	    
	    // get new parameters
	    for (ip = 0; ip < Np; ip++)
		Pnew[ip] = P[ip] * (1.0 + Stfac * Db[ip]);
	    
	    // new quality function
	    Q=0.0;
	    for (im=0; im<Nm; im++)
	    {
		Funct(Xmeas.row(im), Pnew, F);
		Dly=Ymeas.row(im)-F;    // difference
		for (iy=0; iy<Ny; iy++)
		    Q+=Dly[iy]*W[im][iy]*Dly[iy];
	    }
	    if (Verbose==NLIN_CHATTER) cout<<Itno<<'\t'<<Q<<'\t'<<Lm<<'\n';
	    if (Verbose==NLIN_TALK) cout.put('.');

	    // incr Lm if Q grew
	    if (Grow=(Q>=Qold)) Lm=(Lm>0.0)? 10.0*Lm: 0.01;
	    else	// successful step
	    {
		Lm=(Lm > 1e-6)? 0.1*Lm: 0.0;
		P=Pnew; Qold=Q;
	    }
	}
	while (Grow && Stlen >= Steplim / 10);	// linreg step
	Itno++;
    }
    while (Stlen >= Steplim && Itno <= Itmax);	// main cycle end

    // rest of output
    Tcrit95=tcrit_95(Nf);   // t-statistics
    
    Correl=Jtj*(Q/Nf);	// covariance matrix estimated from last inverted J'WJ
    Sdev=Correl.diag();	// relative variances
    for (ip=0; ip<Np; ip++) Sdev[ip]=sqrt(Sdev[ip]);	// relative SD
    
    double Sij;
    for (ip=0; ip<Np; ip++)	// norm to correlation
    {
	Correl[ip][ip]=1.0;
	for (jp=0; jp<ip; jp++)
	{
	    Sij=Sdev[ip]*Sdev[jp];
	    if (Sij<DBL_EPSILON) Correl[ip][jp]=0.0;
	    else Correl[ip][jp]/=Sij;
	}
    }
    for (ip=0; ip<Np; ip++) Sdev[ip]*=fabs(P[ip]);  // SD
    Q=sqrt(Q/Nf);   // residual deviation
    
    if (Verbose==NLIN_TALK) cout<<"Done\n";
    Itmax=Itno; return(Q);
}
// END of nonlin_reg()

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
    int& Itmax, float Steplim, int Verbose) 
{
    // dim checks
    int Nm=Xmeas.dim(), Np=P.dim(), Nf;
    if (Nm!=Ymeas.dim())
    {
	cerr<<"\n! nonlin11_reg(): "<<Ymeas.dim()<<" measurements in Ymeas instead of "<<Nm<<endl;
	return(-1.0);
    }
    if (Nm!=W.dim())
    {
	cerr<<"\n! nonlin11_reg(): W is "<<W.dim()<<" long instead of "<<Nm<<endl;
	return(-1.0);
    }
    if ((Nf=Nm-Np)<=0)
    {
	cerr<<"\n! nonlin11_reg(): Not enough degrees of freedom ("<<Nf<<")\n";
	return(-1.0);
    }
    
    // less serious trouble (adjust silently)
    if (Np!=Sdev.dim()) Sdev.dim(Np);
    if (Np!=Correl.rno()) Correl.set_size(Np);
    Itmax=abs(Itmax); if (!Itmax) Itmax=100;
    Steplim=fabsf(Steplim);
    
    // adjust verbosity: used to be a param but caused trouble
    if (Verbose<=NLIN_SILENT) Verbose=NLIN_SILENT;	// no chat
    if (Verbose>=NLIN_CHATTER) Verbose=NLIN_CHATTER;	// lots of output

    // some of the local vars
    register int im, ip, jp;
    double Q, Qold;
    
    // get the initial quality 
    double F, Dly;
    Qold=0.0;
    for (im=0; im<Nm; im++)
    {
	F=Funct(Xmeas[im], P);
	Dly=Ymeas[im]-F;    // difference
	Qold+=Dly*W[im]*Dly;
    }
    if (Verbose==NLIN_TALK) cout<<"\nnonlin11_reg():";
    if (Verbose==NLIN_CHATTER) cout<<"\nNonlinear regression:\nItno\tQ\tLm\n";

    // set up matrices and vectors
    Vector_ Grad(Np);    // Gradient: dF/dP
    Vector_ Pnew(Np), Pd(Np), JtDy(Np), Db(Np), Norm(Np);
    Trimat_ Jtj(Np), Jtjold(Np);
    int Itno=0;
    double Dy, Fd, Lm=0.01, Stlen, Stfac, Dp, Trace;
    char Grow;
    
    const float DERIV_COEFF=0.001, DERIV_MINSTEP=1E-6;
    
    do	// main iteration cycle
    {
	Jtj.set_values(); JtDy.set_values();	// zeroing
	
	// calculate Jtj=J'WJ and JtDy=J'W(Y-F(X,P)) from observed data
	for (im=0; im<Nm; im++)	    // all measurement points
	{
	    F=Funct(Xmeas[im], P);
	    Dy=Ymeas[im]-F; // obs. vs. calc.
	    
	    // numerical derivation for the Jacobian
	    for (ip=0; ip<Np; ip++)
	    {
		Dp = DERIV_COEFF * fabs(P[ip]) + DERIV_MINSTEP;
		Pd=P; Pd[ip]+=Dp;   // shift ip-th parameter
		Fd=Funct(Xmeas[im], Pd);	// funct values at shifted pos
		Grad[ip]=(Fd-F)*(P[ip]/Dp);	// normalised gradient
	    }
	    
	    // now J'WJ and the "right-hand side" JtDy
	    for (ip = 0; ip < Np; ip++)
	    {
		for (jp = 0; jp <= ip; jp++)
		    Jtj[ip][jp] += Grad[ip] * W[im] * Grad[jp];
		JtDy[ip] += Grad[ip] * W[im] * Dy;
	    }
	}	// for im
	
	// normalise J'WJ
	Trace=Jtj.get_trace();
	Trace/=(Np*1000);   // no norm with pivots smaller than this
	Norm=Jtj.diag();    // copy diagonal
	for (ip=0; ip<Np; ip++)
	    Norm[ip]=(Norm[ip]>Trace)? sqrt(Norm[ip]): 1.0;
	
	for (ip=0; ip<Np; ip++)
	    for (jp=0; jp<=ip; jp++)
		Jtj[ip][jp]/=(Norm[ip]*Norm[jp]);
	Jtjold=Jtj;	// preserve original 
	
	do  // a linreg step, repeat w/ higher Marquardt Lm if Q grows
	{
	    // apply current Marquardt-lambda
	    for (ip=0; ip<Np; ip++) Jtj[ip][ip]+=Lm;
	    
	    // invert Jtj
	    while (!posdef_inv(Jtj))
	    {
		Jtj=Jtjold;
		Lm = (Lm==0.0)? 0.01: 10.0*Lm;
		for (ip=0; ip<Np; ip++) Jtj[ip][ip]+=Lm;
	    }
	    
	    // norm inverse back
	    for (ip=0; ip<Np; ip++)
		for (jp=0; jp<=ip; jp++)
		    Jtj[ip][jp]/=(Norm[ip]*Norm[jp]);
	    Db=Jtj*JtDy;    // linreg step solution
	    
	    // check stepsize
	    Stlen=Db.vec_len();	// overall length
	    Stfac=1.0;
	    for (ip=0; ip<Np; ip++)
		if (Stfac*Db[ip]<-0.95) Stfac=-0.95/Db[ip];
	    Db*=Stfac; Stlen*=Stfac;
	    
	    // get new parameters
	    for (ip = 0; ip < Np; ip++)
		Pnew[ip] = P[ip] * (1.0 + Stfac * Db[ip]);
	    
	    // new quality function
	    Q=0.0;
	    for (im=0; im<Nm; im++)
	    {
		F=Funct(Xmeas[im], Pnew);
		Dly=Ymeas[im]-F;    // difference
		Q+=Dly*W[im]*Dly;
	    }
	    if (Verbose==NLIN_CHATTER) cout<<Itno<<'\t'<<Q<<'\t'<<Lm<<'\n';
	    if (Verbose==NLIN_TALK) cout.put('.');

	    // incr Lm if Q grew
	    if (Grow=(Q>=Qold)) Lm=(Lm>0.0)? 10.0*Lm: 0.01;
	    else	// successful step
	    {
		Lm=(Lm > 1e-6)? 0.1*Lm: 0.0;
		P=Pnew; Qold=Q;
	    }
	}
	while (Grow && Stlen >= Steplim / 10);	// linreg step
	Itno++;
    }
    while (Stlen >= Steplim && Itno <= Itmax);	// main cycle end

    // rest of output
    Tcrit95=tcrit_95(Nf);   // t-statistics
    
    Correl=Jtj*(Q/Nf);	// covariance matrix estimated from last inverted J'WJ
    Sdev=Correl.diag();	// relative variances
    for (ip=0; ip<Np; ip++) Sdev[ip]=sqrt(Sdev[ip]);	// relative SD
    
    double Sij;
    for (ip=0; ip<Np; ip++)	// norm to correlation
    {
	Correl[ip][ip]=1.0;
	for (jp=0; jp<ip; jp++)
	{
	    Sij=Sdev[ip]*Sdev[jp];
	    if (Sij<DBL_EPSILON) Correl[ip][jp]=0.0;
	    else Correl[ip][jp]/=Sij;
	}
    }
    for (ip=0; ip<Np; ip++) Sdev[ip]*=fabs(P[ip]);  // SD
    Q=sqrt(Q/Nf);   // residual deviation
    
    if (Verbose==NLIN_TALK) cout<<"Done\n";
    Itmax=Itno; return(Q);
}
// END of nonlin11_reg()

// ---- Auxiliaries ----

/* tcrit_95() returns the value of the Nf degrees-of-freedom t-distribution
 *  at 95% significance level. Good for confidence intervals. */
static float tcrit_95(long Nf)
{
  double t, TEMP;

  if (Nf <= 5)
    t = 6.415 + (5 - Nf) *
	  (0.289 + (4 - Nf) * (0.0575 + (3 - Nf) * (0.022 + (2 - Nf) * 0.020552)));
  else if (Nf <= 30)
    t = 7.6278 - 0.2316 * Nf + 0.00421 * Nf * Nf - 0.186 * sin(0.2214 * Nf) -
	0.0116 * sin(0.4428 * Nf) + 0.0186 * cos(0.4428 * Nf);
  else if (Nf < 40)
    t = 4.4067 + 0.0296 * (30 - Nf);
  else if (Nf < 60)
    t = 4.1108 + 0.021 * (40 - Nf);
  else
    t = 3.6888 + 0.0116 * (60 - Nf);
  t = exp(t);
  if (Nf == 7)
    t += 2;
  if ((unsigned long)Nf < 32 && ((1L << Nf) & 0x100a100L) != 0)
    t++;
  modf(t, &TEMP);
  return ((float)(TEMP / 1000.0));   /* round to x.xxx */
}
// END of tcrit_95()

/* posdef_inv: inverts the positive definite symmetric matrix given
 * in A. Since the inverse is symmetric, too, it is returned in A
 * (original values are overwritten). 
 * Return value: 1 if OK, 0 if singular.
 */
static int posdef_inv(Trimat_& A)
{
  int i, j, k, N=A.rno();
  double *H=new double [N];
  double At;

  /* main elimination cycle */
  for (k = N; k >= 1; k--)
  {  
    if (fabs(A[0][0]) < DBL_MIN) { delete [] H; return(0);}    /* singularity */
    H[N - 1] = 1 / A[0][0];
    for (i = 2; i <= N; i++)
    {
      At = A[i - 1][0] * H[N - 1];
	H[i - 2] = (i>k)? At: -At;
      for (j = 2; j <= i; j++)
	A[i - 2][j - 2] = A[i - 1][j - 1] + A[i - 1][0] * H[j - 2];
    }
    for (i = 0; i < N; i++)
      A[N - 1][i] = H[i];
  }
  delete [] H; return(1);
}
// END of posdef_inv()

// ==== END OF FUNCTIONS Pmest.c++ ====
