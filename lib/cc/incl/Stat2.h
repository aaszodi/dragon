#ifndef STAT2_CLASSES
#define STAT2_CLASSES

// ==== HEADER Stat2.h ====

/* Very simple 1 and 2-variable statistics (avg,SD,corr). */

// SGI C++ 3.2.1, IRIX 5.2, 17. Jan. 1995. Andris

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <math.h>
#include <limits.h>

// ---- DEFINITIONS ----

#define BIG_NUMBER 1.7e308  /* this should do */
#ifndef DBL_EPSILON
#define DBL_EPSILON 1e-15
#endif

// ==== CLASSES ====

/* Class Stat_: one-variable statistics. A Stat_ object works
 * much in the same way as the stat operations on a programmable
 * calculator. Double data are accumulated in "stat sums" and
 * the average and SD can be obtained from these.
 */
class Stat_
{
    // data
    double Sx, Sx2, Min, Max;	// sum,sum of squares,minimal and maximal value
    unsigned int N;	// number of data
    
    // methods
    public:
    
	// constructor
    Stat_():
      Sx(0.0),
      Sx2(0.0),
      Min(BIG_NUMBER),
      Max(-BIG_NUMBER),
      N(0)
      {}
    
	// data
    /* clear(): resets the calling object. */
    void clear() { Sx=Sx2=0.0; Min=BIG_NUMBER; Max=-Min; N=0; }
    
    /* Stat+=Val: adds a new data point Val to Stat. Returns calling object */
    Stat_& operator+=(double Val);
    
	// results
    unsigned int data_no() const { return(N); }
    
    /* min(), max(): return the minimal/maximal value or 0.0 if
     * there were no data items (w/ a warning).
     */
    double min() const;
    double max() const;
        
    /* avg(): returns the average. Prints a warning and returns 0.0 if
     * there were no data items at all.
     */
    double avg() const;
    
    /* sd(): returns the standard deviation. Prints a warning and returns 0.0
     * if there were no data items at all. If N=1, 0.0 is returned automatically.
     */
    double sd() const;
};
// END OF CLASS Stat_

/* Class Stat2_: Two-variable simple statistics class. A Stat2_ object
 * contains two Stat_ objects (Xs and Ys) plus a mixed-sum member
 * Sxy for correlation. Data point addition is done by pairs. The
 * one-var stats (min/max/avg, SD) are done via returning Xs or Ys
 * and then using the appropriate Stat_ methods.
 */
class Stat2_
{
    // data
    Stat_ Xs, Ys;   // X and Y fields
    double Sxy;	    // mixed product sum
    
    // methods
    public:
	// constructor
    Stat2_(): Xs(), Ys(), Sxy(0.0) {}
    
	// data
    /* clear(): clear the calling object. */
    void clear() { Xs.clear(); Ys.clear(); Sxy=0.0; }
    
    /* add(): adds a pair of values X,Y to the calling object. */
    void add(double X, double Y) { Xs+=X; Ys+=Y; Sxy+=X*Y; }
    
	// results
    /* data_no(): returns the no. of data points (from Xs, but
     * Xs.N==Ys.N all the time).
     */
    unsigned int data_no() const { return(Xs.data_no()); }
    
    /* xs(), ys(): return the one-stat objects for X and Y. Use
     * the Stat_ methods on these to obtain avg, sd, min, max.
     */
    Stat_ xs() const { return(Xs); }
    Stat_ ys() const { return(Ys); }
    
    /* corr(): returns the correlation coefficient between the X and Y
     * points. Prints a warning and returns 0.0 if there were no data points.
     */
    double corr() const;
};
// END OF CLASS Stat2_

// ==== END OF HEADER Stat2.h ====

#endif	/* STAT2_CLASSES */
