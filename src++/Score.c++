// ==== PROJECT DRAGON: METHODS Score.c++ ====

/* Keeps track of scores, relative changes and exit criteria. */

// SGI C++ 7.1, IRIX 6.2, 4. Apr. 1997. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <iostream.h>
#include <iomanip.h>

// ---- MODULE HEADER ----

#include "Score.h"

// ---- DEFINITIONS ----

#ifndef DBL_EPSILON
#define DBL_EPSILON (1.0e-15)
#endif

// ==== Sco_ METHODS ====

// ---- Access ----

/* score(S): replaces the current score with S and
 * moves it to Previous. If a new score was built up
 * beforehand (using operator+=), it is discarded
 * and the summation facility reset. Returns previous score.
 */
double Sco_::score(double S)
{
    Previous=Current; Current=S;
    Sum=0.0; No=0;
    return(Previous);
}
// END of score()

/* set_noexit(): sets the current and previous scores so that
 * both are higher than Minscore and their relative difference is
 * higher than Minchange. Use on entering loops to avoid early exits.
 */
void Sco_::set_noexit()
{
    Current=((Minscore>=0.0)? 2.0: -2.0)*Minscore;
    Previous=10.0*Minchange*Current;
    Sum=0.0; No=0;
}
// END of set_noexit()

/* is_exit(): returns 1 if Current<Minscore,  
 * returns 2 if the relative change
 * between Current and previous is smaller than Minchange.
 * Returns 0 otherwise.
 */
int Sco_::is_exit() const
{
    if (Current<Minscore) return(1);
    else if (rel_change()<Minchange) return(2);
    else return(0);
}
// END of is_exit()

/* rel_change(): returns the relative change of the current score
 * with respect to the previous score. Always >=0.0
 */
double Sco_::rel_change() const
{
    if (fabs(Previous)>=DBL_EPSILON)
	return(fabs((Current-Previous)/Previous));
    else
	return(fabs(Current-Previous));
}
// END of rel_change()

/* change(): returns -1 if the score went down, +1 if it went up,
 * 0 if the change was less than Minchange.
 */
int Sco_::change() const
{
    double Rch=rel_change();
    if (Rch==0.0 || Rch<Minchange) return(0);
    return(Current<Previous? -1: 1);
}
// END of change()

/* update(): updates the score. The current value becomes the
 * previous, and the new current value is the result of the
 * latest summation divided by the normalisation factor.
 * No update is done and a warning is written to stderr if
 * the norm factor was 0.0.
 * Return value: the new current score.
 */
double Sco_::update()
{
    // ### if (!No) return(Current);	// no += since last update
    if (Norm==0.0)
    {
	cerr<<"\n? Sco_::update(): Norm factor is 0.0\n";
	return(0.0);
    }
    Sum/=Norm; No=0;
    Previous=Current; Current=Sum; Sum=0.0;
    return(Current);
}
// END of update()

// ==== Scores_ STATIC INIT ====

const int Scores_::SCO_NO=5;		// 5 different scores
const double Scores_::MAX_RELINCR=0.1;	// 10 % increase allowed

// ==== Scores_ METHODS ====

/* accept_new(): this method compares the current scores in Snew
 * to those in the calling object and returns 1 if Snew represents
 * a score set "more acceptable" than the calling object's set, 
 * 0 otherwise.
 * The acceptance criteria: the BOND and RESTR scores must not grow, 
 * only one of the NONBD|SECSTR score may grow, and if one of
 * them grows then the relative change must be less than the
 * static const MAX_RELINCR. ACCESS may go up by twice of that amount.
 */
int Scores_::accept_new(const Scores_& Snew) const
{
    // these scores must go down
    if (Scos[BOND].score()<Snew[BOND].score())
	return(0);	// bond score would grow, no good
    if (Scos[RESTR].score()<Snew[RESTR].score())
	return(0);	// restraint score would grow, no good
    
    // test the growth of the access score
    Sco_ Sacc=Scos[ACCESS];
    Sacc.score(Snew[ACCESS].score());
    if (Sacc.change()>0 && Sacc.rel_change()>2.0*MAX_RELINCR)
	return(0);  // access score would grow too much
    
    Sco_ Snonbd=Scos[NONBD], Ssecstr=Scos[SECSTR];
    Snonbd.score(Snew[NONBD].score());
    Ssecstr.score(Snew[SECSTR].score());
    
    int Growno=0;   // counts how many scores have grown: only one may grow

    // if one would grow too much, no matter what happens to the others
    if (Snonbd.change()>0)
    {
	if (Snonbd.rel_change()>MAX_RELINCR) return(0);
	Growno++;
    }
    
    #if 0
    if (Srestr.change()>0)
    {
	if (Growno || Srestr.rel_change()>MAX_RELINCR)
	    return(0);
	Growno++;
    }
    #endif
    
    if (Ssecstr.change()>0)
    {
	if (Growno || Ssecstr.rel_change()>MAX_RELINCR)
	    return(0);
	Growno++;
    }
    
    // must be OK
    return(1);
}
// END of accept_new()

/* <<: nice output, one line only (no endl). */
ostream& operator<<(ostream& Out, const Scores_& S)
{
    Out<<"BD="<<S[Scores_::BOND].score()
	<<", NB="<<S[Scores_::NONBD].score()
	<<", RS="<<S[Scores_::RESTR].score()
	<<", SC="<<S[Scores_::SECSTR].score()
	<<", AC="<<S[Scores_::ACCESS].score();
    return(Out);
}
// END of <<

// ==== END OF METHODS Score.c++ ====
