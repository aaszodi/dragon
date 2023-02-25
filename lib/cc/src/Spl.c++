// ==== FUNCTIONS Spl.c++ ====

/* Classic third-order splines. Based on the "spl" module
 * written in C which is based on a Pascal routine collection
 * written in Hungary ages ago...
 */

// SGI C++ 3.2.1, IRIX 5.2, 10. Jan. 1995. Andris 

// ---- HEADER ----

#include "Spl.h"

// ==== Spl_ METHODS ====

// ---- FULL ARRAY ACCESS ----

/* x_arr(), y_arr(): copy the contents of the argument array
 * into the corresponding X or Y array. A warning is printed
 * if the length of the argument is different from the pre-set
 * spline length (cf. the len() method) and no action is taken
 * in order to preserve the lengths of the spline member arrays.
 */
void Spl_::x_arr(const Array_<double>& Xa)
{
    if (X.len()!=Xa.len())
    {
	cerr<<"? Spl_::x_arr(): Array length mismatch\n";
	return;
    }
    X=Xa; Ev=0;
}

void Spl_::y_arr(const Array_<double>& Ya)
{
    if (Y.len()!=Ya.len())
    {
	cerr<<"? Spl_::y_arr(): Array length mismatch\n";
	return;
    }
    Y=Ya; Ev=0;
}
// END of x_arr(), y_arr()

// ---- FITTING AND EVALUATION ----

/* fit_spl(): fits a cubic spline to a series of data points already
 * stored in the calling object. Does not do anything if Ev=1.
 *   yp1,ypn: user-supplied first derivatives at the 1st and last
 *            points. If they are >=SPL_MAX1DER then the second derivatives are
 *            assumed to be 0 (natural spline, default)
 *   Return value: 0 if OK, <0 if something went wrong.
 */
int Spl_::fit_spl(double yp1, double ypn)
{
    if (Ev) return(0);	// fitting already done
    
    int i, k, n=len();
    double p, Qn, Sig, Un, dx;
    double *u;

    // init temp array (conventional)
    u=new double [n];
    if (u==NULL) return(-2);	// out of memory
    
    if (yp1 >=SPL_MAX1DER) Y2[0] = u[0] = 0.0;
    else
    {
	Y2[0] = -0.5;
	u[0] = 3.0 / (X[1] - X[0]) * ((Y[1] - Y[0]) / (X[1] - X[0]) - yp1);
    }
    
    for (i = 2; i < n; i++)
    {
	Sig = (X[i - 1] - X[i - 2]) / (X[i] - X[i - 2]);
	p = Sig * Y2[i - 2] + 2.0;
	Y2[i - 1] = (Sig - 1.0) / p;
	u[i - 1] = (Y[i] - Y[i - 1]) / (X[i] - X[i - 1]) -
		   (Y[i - 1] - Y[i - 2]) / (X[i - 1] - X[i - 2]);
	u[i - 1] = (6.0 * u[i - 1] / (X[i] - X[i - 2]) - Sig * u[i - 2]) / p;
    }
    
    if (ypn >=SPL_MAX1DER) Qn = Un = 0.0;
    else
    {
	Qn = 0.5;
	Un = 3.0 / (X[n - 1] - X[n - 2]) *
	     (ypn - (Y[n - 1] - Y[n - 2]) / (X[n - 1] - X[n - 2]));
    }
    
    Y2[n - 1] = (Un - Qn * u[n - 2]) / (Qn * Y2[n - 2] + 1.0);
    for (k = n - 1; k >= 1; k--)
    Y2[k - 1] = Y2[k - 1] * Y2[k] + u[k - 1];
    
    // calculating the integral
    Yin[0] = 0.0;
    for (k = 2; k <= n; k++)
    {
	dx = X[k - 1] - X[k - 2];
	Yin[k - 1] = dx * (Y[k - 1] + Y[k - 2]) / 2 -
		     dx * dx * dx * (Y2[k - 1] + Y2[k - 2]) / 24 + Yin[k - 2];
    }
    
    delete [] u ; Ev=1;
    return(0);	/* OK */
}
// END of fit_spl()

/* eval_spl: evaluates a spline at a point xi. If xi is out of
 * range then 0.0 is returned. Otherwise the return value is y(xi). 
 * Adapted from Numerical Recipes.
 * The derivatives and the integral are returned only if
 * they are not NULL: this way,  the user can calc any combination
 * of the first three derivatives and the integral.
 *   xi: interpolation point
 *   Der1..3: 1st..3rd derivative of y at xi
 *   Integ: definite integral of y(x) between X[0] and xi.
 * The derivative and integral ptrs are NULL by default.
 * Prints a warning and returns 0.0 if the values were modified (Ev=0).
 */
double Spl_::eval_spl(double xi, double *Der1, double *Der2,
		double *Der3, double *Integ) const 
{
    if (!Ev)	// re-fitting necessary
    {
	cerr<<"? Spl_::eval_spl(): Data modified, please call fit_spl()\n";
	return(0.0);
    }
    
    int kl, kh, k, n=len();
    double a, b, h, yh, yl, y2h, y2l, dx, yi;

    /* locate the position of xi: this is also coded in a separate
     * routine pos_finder() below. 0.0 is returned if xi was out of
     * the range or some other error occurred.
     */
    if (xi > X[n - 1] || xi < X[0])
    {
	if (Der1!=NULL) (*Der1)=0.0;  // out of range
	if (Der2!=NULL) (*Der2)=0.0;
	if (Der3!=NULL) (*Der3)=0.0;
	if (Integ!=NULL) (*Integ)=0.0;
	cerr<<"? Spl_::eval_spl(): X="<<xi<<" is out of range\n";
	return(0.0);
    }

    // find the smallest interval containing xi
    kl = 1;
    kh = n;
    while (kh - kl > 1)
    {
	k = (kh + kl) / 2;
	if (X[k - 1] > xi) kh = k;
	else kl = k;
    }
    h = X[kh - 1] - X[kl - 1];   // length of the interval
    if (h == 0.0)
    {
	if (Der1!=NULL) (*Der1)=0.0;  // zero interval length
	if (Der2!=NULL) (*Der2)=0.0;
	if (Der3!=NULL) (*Der3)=0.0;
	if (Integ!=NULL) (*Integ)=0.0;
	cerr<<"? Spl_::eval_spl(): Zero interval length\n";
	return(0.0);
    }
    
    a = (X[kh - 1] - xi) / h; b = 1 - a;    // useful temps
    yh = Y[kh - 1]; yl = Y[kl - 1];
    y2h = Y2[kh - 1]; y2l = Y2[kl - 1];
    
    // value 
    yi = a * yl + b * yh +
	((a * a * a - a) * y2l + (b * b * b - b) * y2h) * h * h / 6.0;
	
    // derivatives
    if (Der1!=NULL)
	*Der1 = (yh - yl) / h + h * ((1 - 3 * a * a) * y2l + (3 * b * b - 1) * y2h) / 6;
    if (Der2!=NULL)
	*Der2 = a * y2l + b * y2h;
    if (Der3!=NULL)
	*Der3 = (y2h - y2l) / h;
    
    // integral 
    if (Integ!=NULL)
    {
	dx = xi - X[kl - 1];
	*Integ = dx * (yi + yl) / 2 - dx * dx * dx * 
	    (a * y2l + b * y2h + y2l) / 24 + Yin[kl - 1];
    }
    
    return(yi);
}
// END of eval_spl()

/* integ_spl(): returns the definite integral of the spline Spl
 * between Low and Up or 0.0 if range error.
 */
double Spl_::integ_spl(double Low, double Up) const
{
    if (!Ev)	// re-fitting necessary
    {
	cerr<<"? Spl_::integ_spl(): Data modified, please call fit_spl()\n";
	return(0.0);
    }
    
    int n=len();
    if (Low < X[0] || Low > X[n - 1] || Up < X[0] || Up > X[n - 1])
    {
	cerr<<"? Spl_::integ_spl(): limit(s) out of range\n";
	return(0.0);
    }
    
    if (Low <= Up)
    return (int_0x(Up) - int_0x(Low));
    else
    return (int_0x(Low) - int_0x(Up));
}
// END of integ_spl()

/* int_0x(): returns the definite integral of the spline between X[0] and xi.
 * xi is assumed to be within the legal range (this method is called
 * by integ_spl() only and is private).
 */
double Spl_::int_0x(double xi) const
{
    int kl, kh, k, n=len();
    double a, b, h, yi, yh, yl, Der2, y2h, y2l, dx;
    
    // find the smallest interval containing xi
    kl = 1; kh = n;
    while (kh - kl > 1)
    {
	k = (kh + kl) / 2;
	if (X[k - 1] > xi) kh = k; else kl = k;
    }
    
    h = X[kh - 1] - X[kl - 1];   // length of the interval
    if (h == 0.0) return(0.0);

    a = (X[kh - 1] - xi) / h; b = 1 - a;
    yh = Y[kh - 1]; yl = Y[kl - 1];
    y2h = Y2[kh - 1]; y2l = Y2[kl - 1];
    
    yi = a * yl + b * yh +
       ((a * a * a - a) * y2l + (b * b * b - b) * y2h) * h * h / 6.0;
    // value
    Der2 = a * y2l + b * y2h;
    dx = xi - X[kl - 1];
    return (dx * (yi + yl) / 2 - dx * dx * dx * (Der2 + y2l) / 24 + Yin[kl - 1]);
}
// END of int_0x()

// ==== END OF Spl_ METHODS ====

// ==== END OF FUNCTIONS Spl.c++ ====
