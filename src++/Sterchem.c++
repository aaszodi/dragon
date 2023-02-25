// ==== PROJECT DRAGON: FUNCTIONS Sterchem.c++ ====

/* General stereochemical adjustment routines. */

// SGI C++ 4.0, IRIX 5.3, 22. May 1996. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <iostream.h>
#include <iomanip.h>

// ---- MODULE HEADERS ----

#include "Sterchem.h"

// ==== FUNCTIONS ====

// ---- Ideal secondary structure ----

/* apply_secstruct(): RMS fits all secondary structure elements in
 * Pieces onto Model, if it is 3-dimensional. 
 * Returns maximal RMS value if the average is closer to the maximum
 * than to the minimum; returns the average otherwise.
 */
double apply_secstruct(const Pieces_& Pieces, Points_& Model)
{
    double Rms, Maxrms=0.0, Avgrms=0.0, Minrms=HUGE_VAL;
    
    Clist1_<Sstr_> Slist=Pieces.secs();    // secondary struct list iterator
    if (Slist.len())
    {
	for (Slist.begin(); Slist!=NULL; Slist++)
	{
	    Rms=(*Slist)->ideal_struct(Model);
	    if (Rms<0.0) continue;
	    
	    if (Rms>Maxrms) Maxrms=Rms;
	    if (Rms<Minrms) Minrms=Rms;
	    Avgrms+=Rms;
	}
	Avgrms/=Slist.len();	// cheat if there was Rms<0
    }
    return((Maxrms-Avgrms<Avgrms-Minrms)? Maxrms: Avgrms);
}
// END of apply_secstruct()

// ---- Overall handedness check ----

/* hand_check(): Works for 3D molecules with secondary structure only.
 * Checks the torsion angles in the secstr regions of the Model and
 * flips the structure if there were more bad than good angles.
 * Return value: 1 if the original model was kept, -1 if the mirror image.
 */
int hand_check(const Pieces_& Pieces, Points_& Model)
{
    // do nothing in high dims or no-secstr cases
    if (Model.dim()!=3 || !Pieces.hbond_bits().on_no())
	return(0);
    
    register unsigned int Goodstr=0, Badstr=0, 
	    Goodsum=0, Badsum=0, Good, Bad;
    register int Flip;
    
    // test secstr regions
    Clist1_<Sstr_> Slist=Pieces.secs();
    for (Slist.begin(); Slist!=NULL; Slist++)
    {
	Flip=(*Slist)->check_torsion(Model, Good, Bad);
	if (!Flip) continue;
	if (Flip>0) ++Goodstr; else ++Badstr;
	Goodsum+=Good; Badsum+=Bad;
    }
    
    cout<<"HAND: (secstr) Good:Bad="<<Goodstr<<":"<<Badstr<<" ("<<Goodsum<<":"<<Badsum<<")\n";
    if (Goodsum<Badsum || Goodsum==Badsum && Goodstr<Badstr)
    {
	Model.mask(true);
	for (unsigned int i=0; i<Model.len(); i++)
	    Model[i][0]*=-1.0;
	return(-1);
    }
    else return(1);
}
// END of hand_check()

// ==== END OF FUNCTIONS Sterchem.c++ ====
