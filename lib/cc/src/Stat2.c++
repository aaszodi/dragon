// ==== FUNCTIONS Stat2.c++ ====

/* Very simple 1 and 2-variable statistics (avg,SD,corr). */

// SGI C++ 3.2.1, IRIX 5.2, 17. Jan. 1995. Andris

// ---- HEADER ----

#include "Stat2.h"

// ==== FUNCTIONS ====

// ---- Stat_ METHODS ----

/* Stat+=Val: adds a new data point Val to Stat. Returns calling object */
Stat_& Stat_::operator+=(double Val)
{
    Sx+=Val; N++;
    if (Val<Min) Min=Val;
    if (Val>Max) Max=Val;
    Sx2+=(Val*Val);
    return(*this);
}
// END of +=

/* min(), max(): return the minimal/maximal value or 0.0 if
 * there were no data items (w/ a warning).
 */
double Stat_::min() const
{
    if (!N)
    {
	cerr<<"? Stat_::min(): no data items\n";
	return(0.0);
    }
    else return(Min);
}

double Stat_::max() const
{
    if (!N)
    {
	cerr<<"? Stat_::max(): no data items\n";
	return(0.0);
    }
    else return(Max);
}
// END of min(), max()

/* avg(): returns the average. Prints a warning and returns 0.0 if
 * there were no data items at all.
 */
double Stat_::avg() const
{
    if (!N)
    {
	cerr<<"? Stat_::avg(): no data items\n";
	return(0.0);
    }
    else return(Sx/N);
}
// END of avg()

/* sd(): returns the standard deviation. Prints a warning and returns 0.0
 * if there were no data items at all. If N=1, 0.0 is returned automatically.
 */
double Stat_::sd() const
{
    if (!N)
    {
	cerr<<"? Stat_::sd(): no data items\n";
	return(0.0);
    }
    if (N==1) return(0.0);
    double Xavg=Sx/N;
    double Numer=Sx2-N*Xavg*Xavg;
    return((fabs(Numer)<DBL_EPSILON)? 0.0: sqrt(Numer/N));
}
// END of sd()

// ---- Stat2_ METHODS ----

/* corr(): returns the correlation coefficient between the X and Y
 * points. Prints a warning and returns 0.0 if there were no data points.
 */
double Stat2_::corr() const
{
    unsigned int N;
    if ((N=data_no())<2)
    {
	cerr<<"? Stat2_::corr(): not enough data items\n";
	return(0.0);
    }
    double Sd=Xs.sd()*Ys.sd();
    if (Sd<DBL_EPSILON) return(0.0);	// no corr if no SD
    return((Sxy-N*Xs.avg()*Ys.avg())/(N*Sd));
}
// END of corr()

// ==== END OF FUNCTIONS Stat2.c++ ====
