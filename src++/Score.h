#ifndef SCORE_CLASSES
#define SCORE_CLASSES

// ==== PROJECT DRAGON: HEADER Score.c++ ====

/* Keeps track of scores, relative changes and exit criteria. */

// SGI C++ 7.1, IRIX 6.2, 2. Apr. 1997. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <math.h>

// ==== CLASSES ==== 

/* Class Sco_ : holds current and previous score, an absolute
 * and a relative limit. The idea is that in an iteration the
 * exit criterion could be either that the current score should be
 * lower than a preset limit or that the relative change should
 * be smaller than another limit.
 */
class Sco_
{
    // data
    private:
    
    double Sum, Norm;		// for summing scores and normalisation factor
    double Current, Previous;	// current and previous score
    double Minscore, Minchange;	// the absolute and relative limits
    int No;	// number of additions to Sum since last reset
    
    // methods
    public:
    
	// constructors
    /* Inits to hold a minimal limit Minlim and a minimal change Minchg.
     * The defaults for both are 0.0. Sets the object to "no-exit" status.
     */
    Sco_(double Minlim=0.0, double Minchg=0.0):
	    Sum(0.0), Norm(1.0), No(0), 
	    Minscore(Minlim), Minchange(fabs(Minchg)) { set_noexit(); }
    
	// access
    /* score(): returns the current score.
     * score(S): replaces the current score with S and
     * moves it to Previous. If a new score was built up
     * beforehand (using operator+=), it is discarded
     * and the summation facility reset. Returns previous score.
     */
    double score() const { return(Current); }
    double score(double S);

    // sum_reset(): zeroes the summation facility.
    void sum_reset() { Sum=0.0; No=0; }
    
    // +=: updates the sum of scores with V.
    Sco_& operator+=(double V) { Sum+=V; No++; return(*this); }
    
    // norm(): adjust the norm for the score (default 1.0). No query
    void norm(double N=1.0) { Norm=(N<=0.0)? 1.0: N; }
    
    /* min_score(), min_change(): adjust the absolute and relative
     * limits. There is no way to query these.
     */
    void min_score(double Minsco) { Minscore=Minsco; }
    void min_change(double Minchg) { Minchange=fabs(Minchg); }
    
    /* set_noexit(): sets the current and previous scores so that
     * both are higher than Minscore and their relative difference is
     * higher than Minchange. Use on entering loops to avoid early exits.
     */
    void set_noexit();
    
    /* is_exit(): returns 1 if Current<Minscore,  
     * returns 2 if the relative change
     * between Current and previous is smaller than Minchange.
     * Returns 0 otherwise.
     */
    int is_exit() const;

    /* rel_change(): returns the relative change of the current score
     * with respect to the previous score. Always >=0.0
     */
    double rel_change() const;

    /* change(): returns -1 if the score went down, +1 if it went up,
     * 0 if the change was less than Minchange.
     */
    int change() const;

    /* update(): updates the score. The current value becomes the
     * previous, and the new current value is the result of the
     * latest summation divided by the normalisation factor.
     * No update is done and a warning is written to stderr if
     * the norm factor was 0.0.
     * Return value: the new current score.
     */
    double update();
};
// END OF CLASS Sco_

/* Scores_: this object holds a set of various Sco_ subobjects
 * each representing a different kind of score. Score types are defined
 * in the Scotype_ enum. The sub-scores may be accessed separately, 
 * and there are some operations which can be carried out on
 * all of them at the same time. 
 */
class Scores_
{
    public:
    
    /* Type of available scores (cf. Viol.h) */
    enum Scotype_ {BOND=0, NONBD, RESTR, ACCESS, SECSTR};
    
    // data
    private:
    
    static const int SCO_NO;	// number of sub-scores
    static const double MAX_RELINCR;	// maximal relative increase
    
    Sco_ *Scos;	    // array of sub-scores
    
    // methods
    public:
    
	// constructors and destructor, assignment
    /* Init all sub-scores to hold the minimal score limit Minsco
     * and the minimal relative change Minchg (default 0.0 for both).
     */
    Scores_(double Minsco=0.0, double Minchg=0.0): Scos(new Sco_ [SCO_NO])
    {
	for (register int i=0; i<SCO_NO; i++)
	{
	    Scos[i].min_score(Minsco);
	    Scos[i].min_change(Minchg);
	}
    }
    
    // copy constructor
    Scores_(const Scores_& S): Scos(new Sco_ [SCO_NO])
    { for (register int i=0; i<SCO_NO; i++) Scos[i]=S.Scos[i]; }
    
    // destructor
    ~Scores_() { delete [] Scos; }

    // assignment
    Scores_& operator=(const Scores_& S)
    {
	if (this!=&S) 
	    for (register int i=0; i<SCO_NO; i++) Scos[i]=S.Scos[i];
	return(*this);
    }
    
	// access
    const Sco_& operator[](Scotype_ s) const { return(Scos[s]); }
    Sco_& operator[](Scotype_ s) { return(Scos[s]); }
    
	// global operations
    /* The following methods execute the corresponding Sco_
     * methods with the same name for all sub-score members.
     * See comments in Sco_.
     */
    void min_score(double Minsco) 
	{ for (int i=0; i<SCO_NO; i++) Scos[i].min_score(Minsco); }
    void min_change(double Minchg)
	{ for (int i=0; i<SCO_NO; i++) Scos[i].min_change(Minchg); }
    void set_noexit()
	{ for (int i=0; i<SCO_NO; i++) Scos[i].set_noexit(); }
    int is_exit()
    {
	int Isexit=1;
	for (int i=0; i<SCO_NO; i++) Isexit=Isexit && Scos[i].is_exit();
	return(Isexit);
    }
    
    /* update(): does Sco_::update() for all member scores.
     * update(S): moves the current scores from S to the calling object.
     */
    void update()
	{ for (int i=0; i<SCO_NO; i++) Scos[i].update(); }
    void update(const Scores_& S)
	{ for (int i=0; i<SCO_NO; i++) Scos[i].score(S.Scos[i].score()); }
    
    /* accept_new(): this method compares the current scores in Snew
     * to those in the calling object and returns 1 if Snew represents
     * a score set "more acceptable" than the calling object's set, 
     * 0 otherwise.
     * The acceptance criteria: the BOND and RESTR scores must not grow, 
     * only one of the NONBD|SECSTR score may grow, and if one of
     * them grows then the relative change must be less than the
     * static const MAX_RELINCR. ACCESS may go up by twice of that amount.
     */
    int accept_new(const Scores_& Snew) const;

	// output
    /* <<: nice output, one line only (no endl). */
    friend ostream& operator<<(ostream& Out, const Scores_& S);
};
// END OF CLASS Scores_

// ==== END OF HEADER Score.h ====
#endif	/* SCORE_CLASSES */
