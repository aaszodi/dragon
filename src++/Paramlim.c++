#ifndef PARAMLIM_TMPL_DEFS
#define PARAMLIM_TMPL_DEFS

// ==== PROJECT DRAGON: TEMPLATE METHODS Paramlim.c++ ====

/* Stores global numeric parameters. */

// SGI C++ 7.1, IRIX 6.2, 7. Mar. 1997. Andris Aszodi

// ==== Paramlim_ METHODS ====

// ---- Constructor ----

/* Inits to hold the default value Defval which falls between
 * the lower and upper bounds L and U. The name and description
 * strings are specified in Nm and Ds (both defaults to NULL).
 * If L>U then the values are silently swapped. If Defval is 
 * outside the range defined by [L..U], then it is appropriately
 * modified w/o warning.
 */
template <class T_>
Paramlim_<T_>::Paramlim_(const T_& Defval, const T_& L, const T_& U, 
	const char *Nm, const char *Ds):
	    Default(Defval), Parambase_(Nm, Ds)
{
    if (L>U) { Low=U; Up=L; } else { Low=L; Up=U; } // save limits
    if (Default<Low) Default=Low;
    if (Default>Up) Default=Up;
    Value=Default;
}

// ---- Access ----

/* set_deflims(): resets the default value and the limits
 * (similar to the ctor). The value will be reset to the new
 * default value.
 */
template <class T_>
void Paramlim_<T_>::set_deflims(const T_& Defval, const T_& L, const T_& U)
{
    if (L>U) { Low=U; Up=L; } else { Low=L; Up=U; } // save limits
    Default=Defval;
    if (Default<Low) Default=Low;
    if (Default>Up) Default=Up;
    Value=Default; Changed=true;
}
// END of set_deflims()

// ---- Input/output ----

/* read_from(): reads the following line from an input stream In:
 * "NAME value"
 * where NAME should match the Name member of the object and value
 * should conform to the actual type of Value (T_). If "value" is
 * outside the limits set by Low and Up and it will be modified
 * silently. 
 */
template <class T_>
int Paramlim_<T_>::read_from(istream& In)
{
    if (!In.good()) return(-1);
    
    String_ Namebuf(Name);
    Namebuf[0]='\0'; // make sure sizes match
    streampos Pos=In.tellg();	// save initial position
    
    In>>Namebuf;
    if (!In.good() || Namebuf!=Name)
    {
	In.seekg(Pos);	// name did not match, reset and return
	return(0);
    }
    
    In.flags(ios::skipws);
    In>>Value;
    if (!In.good())
    {
	cerr<<"\n? Param_::read_from(): Cannot read value for variable "<<Name<<endl;
	Value=Default;
	In.clear(ios::failbit|In.rdstate());
	return(-1);
    }
    if (Value<Low) Value=Low;
    if (Value>Up) Value=Up;
    Changed=true;	// value has changed
    return(1);
}
// END of read_from()

/* write_to(): writes the calling object to the output stream Out
 * in the following format:-
 * "# DESCRIPTION_STRING (default DEFAULT, limits LOW..UP)\n"
 * "NAME value\n"
 * which can be parsed by the input routine. The # line is
 * written only if Comments==true (the default).
 */
template <class T_>
ostream& Paramlim_<T_>::write_to(ostream& Out, bool Comments) const
{
    if (Comments)
	Out<<"\n# "<<Descr<<" (default="<<Default<<", limits: ["<<Low<<" .. "<<Up<<"] )\n";
    Out<<Name<<" "<<Value<<"\n";
    return(Out);
}
// END of write_to()

// ==== END OF TEMPLATE METHODS Paramlim.c++ ====

#endif	/* PARAMLIM_TMPL_DEFS */
