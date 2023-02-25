// ==== PROJECT DRAGON: METHODS Viol.c++ ====

/* For keeping track of restraint violations. */

// SGI C++ 7.1, IRIX 6.2, 2. Apr. 1997. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <fstream.h>

// ---- MODULES ----

#include "Viol.h"

// ==== METHODS ====

// ==== Viol_ METHODS ====

ostream& operator<<(ostream& Out, const Viol_& V)
{
    int Oprec=Out.precision(2);	// store old precision
    long Oformat=Out.flags();	// and format flags
    
    Out.setf(ios::fixed, ios::floatfield);
    Out<<setw(3)<<V.Atom1<<"["<<setw(4)<<V.Res1<<"]:"
	<<setw(3)<<V.Atom2<<"["<<setw(4)<<V.Res2<<"] ";
    switch(V.Vtype)
    {
	case Viol_::BOND: Out<<" BOND"; break;
	case Viol_::NONBD: Out<<"NONBD"; break;
	case Viol_::RESTR: Out<<"RESTR"; break;
	case Viol_::HELIX: Out<<"HELIX"; break;
	case Viol_::SHEET: Out<<"SHEET"; break;
	case Viol_::UNDEF:
	default: Out<<"UNDEF"; break;
    }
    Out<<" "<<setw(5)<<V.Actual
	<<((V.Actual<V.Ideal)? " < ":" > ")
	<<setw(5)<<V.Ideal<<" ("<<setw(4)<<V.Strict<<") "
	<<setw(5)<<V.Viol<<" "<<setw(5)<<setprecision(1)
	<<(100.0*V.rel_error())<<" %\n";
    Out.precision(Oprec);   // reset precision
    Out.flags(Oformat);	    // and format
    return(Out);
}

// ==== Viollist_ METHODS ====

/* add_viol(): adds the violation object V to the violation
 * list if the relative error is larger than Minrelv.
 * Return value: 0 if no insertion was performed, non-0 otherwise.
 */
int Viollist_::add_viol(const Viol_& V, float Minrelv)
{
    /* don't do anything if V's relative error is smaller
     * than Minrelv
     */
    if (Minrelv<0.0 || V.rel_viol()<Minrelv) return(0);
    
    if (!Vl.len()) { Vl+=V; return(1); }    // list was empty: add first
    
    // find first item with a rel. violation smaller than the current one
    for (Vl.begin(); Vl!=NULL && Vl->rel_viol()>V.rel_viol(); Vl++);
    
    Vl.insert(V);   // add new item before first smaller
    return(Vl.len());
}
// END of add_viol()

/* write_file(): writes the contents of the calling object to
 * a file Outfile or to cout if Outfile==NULL or cannot be opened.
 * Return value: 1 if OK, 0 if written to cout.
 */
int Viollist_::write_file(const char *Outfile) const
{
    if (Outfile==NULL) { cout<<(*this); return(0); }
    ofstream Outf(Outfile);
    if (!Outf) { cout<<(*this); return(0); }
    Outf<<(*this); Outf.close(); return(1);
}
// END of write_file()

ostream& operator<<(ostream& Out, const Viollist_& Vl)
{
    Clist1_<Viol_> Cvl=Vl.Vl;
    Out<<"# Restraint violations: "<<Cvl.len()<<endl;
    Out<<"#     Atom pair     Type  Actual Ideal (Strict) Rel.viol Error\n";
    for (Cvl.begin(); Cvl!=NULL; Cvl++)
	Out<<(*Cvl);
    return(Out);
}

// ==== END OF METHODS Viol.c++ ====
