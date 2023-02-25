// ==== METHODS Cdf.c++ ====

/* Estimation of cumulative distribution functions. */

// SGI C++ 4.0, IRIX 5.3, 11. May 1995. Andris Aszodi 

// ---- CLASS HEADER ----

#include "Cdf.h"

// ==== Cdf_ METHODS ====

// ---- Constructor ----

/* Inits to have Binno>=2 bins equally spaced between Low and Up.
 * If Low>Up then they're swapped silently. If Binno<=1 then Binno==2 is used.
 */
Cdf_::Cdf_(unsigned int Binno, double Low, double Up):
    X(Binno), Y(Binno), Yint(Binno), N(0), Ev(1)
{
    if (Binno<=1)
    {
	Binno=2; X.len(2); Y.len(2); Yint.len(2);
    }
    
    double Step;
    if (Low>Up) { Step=Low; Low=Up; Up=Step; } // swap
    
    Step=(Up-Low)/(Binno-1);
    for (unsigned int d=0; d<Binno; d++)
    {
	X[d]=Low+d*Step; Y[d]=0.0; Yint[d]=0;
    }
}

// ---- Reset ----

/* reset(): resets the calling object to have Binno>=2 bins equally spaced
 * between Low and Up. Zeroes everything. Very similar to what the ctor does.
 */
void Cdf_::reset(unsigned int Binno, double Low, double Up)
{
    if (Binno<=1) Binno=2;
    X.len(Binno); Y.len(Binno); Yint.len(Binno);
    
    double Step;
    if (Low>Up) { Step=Low; Low=Up; Up=Step; } // swap
    
    Step=(Up-Low)/(Binno-1);
    for (unsigned int d=0; d<Binno; d++)
    {
	X[d]=Low+d*Step; Y[d]=0.0; Yint[d]=0;
    }
    N=0; Ev=1;
}
// END of reset()

// ---- Data update ----

/* Cdf+=V, Cdf-=V: add or remove a value V to/from the Cdf object. */
Cdf_& Cdf_::operator+=(double V)
{
    if (V>=X[X.len()-1]) return(*this);	    // too large
    
    register unsigned int Dup=get_index(V);
    Yint[Dup]++; N++; Ev=0;
    return(*this);
}
// END of +=

Cdf_& Cdf_::operator-=(double V)
{
    if (V>=X[X.len()-1]) return(*this);	    // too large
    
    register unsigned int Dup=get_index(V);
    if (Yint[Dup]>=1)
    {
	Yint[Dup]--; N--; Ev=0;  // cannot subtract if already 0
    }
    return(*this);
}
// END of +=

/* get_index(): finds the index kh which satisfies X[kh-1]<V<=X[kh]
 * (X[] is assumed to hold values in ascending order).
 * Return value: kh or 0 if V<X[0] or Len-1 if V>X[Len-1]. Private
 */
unsigned int Cdf_::get_index(double V) const
{
    
    // out of range: return first or last index
    unsigned int Len=X.len();
    if (V<X[0]) return(0);
    if (V>X[Len-1]) return(Len-1);
    
    /* search */
    register unsigned int kl, k, kh;
    kl=0; kh=Len;
    while (kl<kh-1)
    {
	k=(kl+kh)/2;
	if (V==X[k]) { kh=k; break; }
	if (V<X[k]) kh=k; else kl=k;
    }
    return(kh);
}
// END of get_index()

// ---- Access ----

const Array_<double>& Cdf_::y_arr()
{
    if (!Ev)
    {
	eval_cdf(); Ev=1;
    }
    return(Y); 
}
// END of y_arr()

Vector_ Cdf_::x_vec() const
{
    Vector_ Xv(X.len());
    for (register unsigned int i=0; i<X.len(); i++) Xv[i]=X[i];
    return(Xv); 
}
// END of x_vec()

Vector_ Cdf_::y_vec()
{
    if (!Ev) { eval_cdf(); Ev=1; }
    Vector_ Yv(Y.len());
    for (register unsigned int i=0; i<Y.len(); i++) Yv[i]=Y[i];
    return(Yv); 
}
// END of y_vec()

/* eval(): calculates the CDF and puts the estimated values 
 * in Y (a normalisation). Leaves the original Yint[]
 * intact. Private
 * Return value: 1 if OK, 0 if no action was taken, -1 if no data were in Cdf.
 */
int Cdf_::eval_cdf()
{
    if (!N) return(-1);	    // no data
    if (Ev) return(0);	    // has already been evaluated ("clean")
    
    register unsigned int d, Cum;
    register double Total=double(N);   /* convert to real */
    
    Cum=0;
    for (d=0; d<X.len(); d++)
    {
	Cum+=Yint[d];
	Y[d]=Cum/Total;
    }
    return(1);
}
// END of eval_cdf()

// ==== END OF METHODS Cdf.c++ ====
