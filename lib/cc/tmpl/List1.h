#ifndef LIST1_TMPL_DECLS
#define LIST1_TMPL_DECLS

// ==== TEMPLATE HEADER List1.h ====

/* Unidirectional linked list templates. */

// SGI C++ 4.0, IRIX 5.3, 7. Mar. 1996. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>

// ==== CLASSES ====

/* This is a pretty complex bunch of classes. The list items are
 * kept in Item1_ objects linked together and the handler class
 * Listhnd1_ keeps track of the head, tail and length of them.
 * The handler objects are completely private. Two public wrapper interfaces
 * are provided: Clist1_, which accesses the list in a non-destructive
 * way (read item at current position, move current position pointer Cur), 
 * and List1_ (derived from Clist1_), which supports destructive
 * manipulations (insert/delete, write item at current position)
 * as well. More than one Clist1_ or List1_ wrappers can access
 * the same Listhnd1_ through Lptr pointers; the handler counts 
 * the wrappers which point to it and dies only if all wrappers are dead. 
 * Here is an illustration of the situation, with just one 
 * Clist1_ accessing Listhnd1_ (but there could be several, 
 * and List1_ objects as well):-
 * 
 *   Item1_ --> Item1_ --> ... ---> Item1_ ---> Item1_
 *     A                              A           A
 *     |  Head  +-----------+         |    Tail   |
 *     +------- | Listhnd1_ | --------------------+
 *              +-----------+         |
 *                    A               |
 *                    | Lptr          |
 *              +-----------+   Cur   |
 *              |  Clist1_  |---------+
 *              +-----------+
 */

// NOTE: for GNU C++'s sake, the Item_ struct which originally
// was a member class of Listhnd1_ is moved outside. 6-Jul-1997

// forward declarations
template <class T_> class Listhnd1_;
template <class T_> class Clist1_;
template <class T_> class List1_;

    /* The items are stored in a container class template Item1_
     * which provides direct access to its members. The Next pointer
     * is NULL if the item is at the end of the list.
     */
    template <class T_>
    class Item1_
    {
        friend class Listhnd1_<T_>;
        friend class Clist1_<T_>;
        friend class List1_<T_>;

        private:
	T_ Value;   // the item
	Item1_ *Next;	// ptr to the next item
	
	// constructor
	Item1_(const T_& Val): Value(Val), Next(NULL) {}
    };
    // END OF CLASSLET Item_

/* Listhnd1_: the list handler class. Holds the ptr-linked
 * list of Item1_ objects between two pointers, Head and Tail
 * and keeps track of the length of the list in Len.
 * This class is inaccessible to the outside world, only
 * Clist1_ and List1_ (its friends) can use it. Both classes
 * access the list handler via a pointer. Listhnd1_ keeps track
 * of the number of Clist1_ and List1_ objects pointing to it
 * and removes the list only if the last of these was removed.
 */
template <class T_>
class Listhnd1_
{
    friend class Clist1_<T_>;
    friend class List1_<T_>;
    
    private:    
    Item1_<T_> *Head, *Tail;    // ptrs to beginning and end of list
    unsigned int Len, Ptrno;	// list length and no. of objects pointing to it
    
    // methods (all private)
    
	// constructors
    /* NOTE: since all ctors are usable by the Clist1_/List1_
     * classes only, we assume that at least 1 ptr will point to
     * the new Listhnd1_ object. 
     */
    /* Init to empty */
    Listhnd1_(): Head(NULL), Tail(NULL), Len(0), Ptrno(1) {}
    
    /* Init to a 1-member list holding the Value */
    Listhnd1_(const T_& Item): Len(1), Ptrno(1)
    {
	Head=Tail=new Item1_<T_>(Item);
    }
    
    /* Copy constructor. */
    Listhnd1_(const Listhnd1_<T_>& L): Len(L.make_copy(Head, Tail)), Ptrno(1) {}
    
	// destructor
    /* The destructor is tricky. It is called whenever a Clist1_/List1_
     * object dies that pointed to the current object. However, there
     * might be others pointing to it: so only the Ptrno is decremented
     * and the actual destruction happens only when Ptrno==1.
     */
    ~Listhnd1_() { if (!(--Ptrno)) remove_all(); }
    
    /* = : Creates an almost exact copy of the rhs in the lhs. Note that the
     * Ptrno is unchanged! This assignment is used in the "deep copy"
     * methods of the List1_ class.
     */
    Listhnd1_<T_>& operator=(const Listhnd1_<T_>& Rhs);
    
	// memory management
    /* make_copy(): constructs a copy of the calling list. The beginning
     * and end of the new list will be returned in Beg and End, the current
     * position on the new list (corresponding to the original Cur) 
     * will be returned in *Cp if it isn't NULL (the default is NULL).
     * Return value: the length of the new list.
     */
    unsigned int make_copy(Item1_<T_>*& Beg, Item1_<T_>*& End) const;

    /* remove_all(): removes all items from the calling object. */
    void remove_all();
    
};
// END OF CLASS Listhnd1_

/* Clist1_: a unidirectional linked list template. Items are
 * accessed via an internal "current position" pointer which
 * can be moved from the Head of the list towards the Tail.
 * Items at the current position can be read but not modified.
 * Use this class for "constant lists".
 */
template <class T_>
class Clist1_
{
    // data
    protected:
    
    Listhnd1_<T_> *Lptr;	// ptr to the list handler
    Item1_<T_> *Cur;  // ptr to current position
    
    // methods
    public:
    
	// constructors
    /* Copy constructor ("shallow"). NOTE: this is the only public ctor!
     * It is assumed that a "non-const" List1_ object has already created
     * a list and the Clist1_ object just accesses that. Use either to
     * shallow-copy another Clist1_ object (the two would scan the same
     * underlying Listhnd1_) or to init with a List1_ object which is
     * derived from Clist1_ so will be acceptable here.
     */
    Clist1_(const Clist1_<T_>& L): Lptr(L.Lptr), Cur(L.Cur) { (Lptr->Ptrno)++; }
    
    protected:
    /* This constructor is to be used by List1_ only!
     * Lptr and Cur are assigned appropriately as determined
     * by the List1_ ctors.
     */
    Clist1_(Listhnd1_<T_> *Lh=NULL, Item1_<T_> *Cp=NULL):
	    Lptr(Lh), Cur(Cp) {}
    public:
    
	// Destructor
    /* Calls the destructor of the underlying list object
     * explicitly which deletes it only if the current Clist_
     * object was the last that pointed to it.
     */
    virtual ~Clist1_()
    {
	Lptr->Listhnd1_<T_>::~Listhnd1_();
	if (!Lptr->Ptrno) delete Lptr;	// remove handler if last dies
    }
        
	// assignment
    /* = : "shallow" assignment. The Listhnd1_ pointed to by the lhs
     * is "released" and the lhs then will point to Rhs's list.
     */
    Clist1_<T_>& operator=(const Clist1_<T_>& Rhs);
    
	// Access
    /* len(): returns the current length. */
    unsigned int len() const { return(Lptr->Len); }
    
    /* !: returns 1 if the list is empty, 0 otherwise. */
    int operator!() const { return(Lptr->Len? 0: 1); }
    
    /* Conversion to const void* : the status of the Cur pointer can be
     * checked via this: NULL is returned when Cur==NULL, indicating
     * that we are at the end of the list. The pointer value obtained
     * by this conversion must not be used for anything else.
     */
    operator const void* () const { return((const void*)Cur); }
    
    /* *(): returns a ref. to the value of the item at the current position or
     * aborts with an error message if Cur==NULL.
     */
    const T_& operator*() const
    {
	if (Cur==NULL)
	{
	    cerr<<"\n? *Clist1_: Illegal const access attempted\n";
	    abort();
	}
	return(Cur->Value);
    }
    
    /* ->: returns a const ptr to the value of the item at the current
     * position or aborts with an error message if Cur==NULL.
     */
    const T_* operator->() const
    {
	if (Cur==NULL)
	{
	    cerr<<"\n? Clist1_->: Illegal const access attempted\n";
	    abort();
	}
	return((const T_*)&(Cur->Value)); // GCC: 7-Jul-1997.
    }
    
	// Move
    /* begin(), end(): moves the current pointer to the head/tail
     * of the list.
     */
    void begin() { Cur=Lptr->Head; }
    void end() { Cur=Lptr->Tail; }
    
    /* Postfix ++: moves the Cur ptr to the next position and returns 1
     * if it was not already at the end in which case it returns 0.
     * Also returns 0 if Cur==NULL. Note that Cur is allowed to walk
     * off the list: if Cur==Tail and ++ is invoked, then 0 is returned
     * but Cur will become NULL.
     */
    int operator++(int)
    {
	return((Cur!=NULL)? Cur=Cur->Next, Cur!=NULL: 0);
    }
    
    /* forward(): moves the current ptr forward by N steps (default 1)
     * if possible. May walk off the list at the end.
     * Returns the actual number of steps taken.
     */
    unsigned int forward(unsigned int N=1)
    {
	register unsigned int Stepno;
	for (Stepno=0; Stepno<N && Cur!=NULL; Stepno++, Cur=Cur->Next);
	return(Stepno);
    }
    
    #ifdef __LIST1_DEBUG__
	friend ostream& operator<<(ostream& Out, const Clist1_<T_>& List);
    #endif

};
// END OF CLASS Clist1_

/* List1_ : the general unidirectional linked list template class.
 * Derived from Clist1_, this class adds "non-constant" operations
 * (insertions, deletions, item modification) to the basic set of
 * access methods.
 */
template <class T_>
class List1_: public Clist1_<T_>
{
    // data
    protected:
    
    // methods
    public:
    
	// constructors
    /* Starts an empty list. */
    List1_(): Clist1_<T_>(new Listhnd1_<T_>()) {}
    
    /* Starts a list with one item in it. */
    List1_(const T_& Item): Clist1_<T_>(new Listhnd1_<T_>(Item)) { Cur=Lptr->Head; }
    
    /* Copy constructor. This is a "deep" copy: the underlying
     * Listhnd1_ of L is cloned.
     */
    List1_(const List1_<T_>& L): Clist1_<T_>(new Listhnd1_<T_>(*(L.Lptr))) { copy_curpos(L); }
    
	// Access
    /* operator *(): non-const access to the current item.
     * Aborts with a warning if Cur==NULL.
     */
    T_& operator*()
    {
	if (Cur==NULL)
	{
	    cerr<<"\n? *Clist1_: Illegal non-const access attempt\n";
	    abort();
	}
	return(Cur->Value);
    }

    /* ->: returns a ptr to the value of the item at the current
     * position or aborts with an error message if Cur==NULL.
     */
    T_* operator->()
    {
	if (Cur==NULL)
	{
	    cerr<<"\n? Clist1_->: Illegal non-const access attempted\n";
	    abort();
	}
	return((T_*)&(Cur->Value)); // GCC 7-Jul-1997
    }
    
	// Assignment
    /* = : Creates an exact copy of the rhs in the lhs (even the relative
     * position of the current pointer is preserved). This is a "deep" copy.
     */
    List1_<T_>& operator=(const List1_<T_>& Rhs)
    {
	if (this==&Rhs) return(*this);	// x=x
	
	*(this->Lptr)=*(Rhs.Lptr);  // deep copy of underlying list
	copy_curpos(Rhs);  // adjust cursor position
	return(*this);
    }
    
	// Insertions
    /* insert(Val): inserts a new item into the calling object so that
     * the current item pointed to by Cur will follow it, i. e.
     * ...->Current->... will become ...->New->Current->... and Cur
     * will point to New. If Cur==NULL, then the new item will be
     * appended to the end of the list (cf. operator+= ).
     * Returns the calling object.
     */
    List1_<T_>& insert(const T_& Val);
    
    /* insert(List): inserts a whole List into the calling object
     * at the current position so that ...->Current->... (where Current
     * is pointed to by Cur) becomes ...->[List]->Current->...
     * and Cur will point to the first item in List. If Cur==NULL, 
     * then List will be appended to the end of the calling object.
     * Returns calling object.
     */
    List1_<T_>& insert(const List1_<T_>& List);
    
    /* List1+=Val : appends a new item with the value Val at the end of the list.
     * List1+=List: appends a copy of List at the end of the calling object.
     * Both return the calling object.
     */
    List1_<T_>& operator+=(const T_& Val);
    List1_<T_>& operator+=(const List1_<T_>& List);
    
    /* List1^=Val: appends a new item with value Val before the beginning of List1.
     * List1^=List: appends List before the first item of List1.
     * The operator ^= was chosen because '^' means "the beginning of a string"
     * in some regular expressions (reminiscent of "beginning of a list).
     * Both return the calling object.
     */
    List1_<T_>& operator^=(const T_& Val);
    List1_<T_>& operator^=(const List1_<T_>& List);
    
	// Deletions
    /* del(): deletes N items (default 1) from the list, beginning with the
     * current item pointed to by Cur. If Cur==NULL or N==0, nothing happens. If N is
     * larger than the length of the tail of the list from Cur, then the
     * whole tail is deleted. The new current position will be the first
     * item after the deleted region.
     * Return value: the number of items actually deleted.
     */
    unsigned int del(unsigned int N=1);
    
    /* clear(): removes all items from the calling object and returns it. */
    List1_<T_>& clear() { Lptr->remove_all(); Cur=NULL; return(*this); }
    
	// Auxiliaries
    private:
    
    /* copy_curpos(): sets the Cur current pointer of the calling object to
     * the position where L's Cur is: e.g. if L.Cur points to the fifth
     * item in its list, then Cur will point also to the fifth item.
     * It is assumed that the two objects have same-length lists inside
     * (used for deep copies only). Private
     */
    void copy_curpos(const List1_<T_>& L)
    {
	Item1_<T_> *Ocp;
	for (Ocp=L.Lptr->Head, Cur=Lptr->Head; 
	    Ocp!=NULL && Cur!=NULL && Ocp!=L.Cur;
	    Ocp=Ocp->Next, Cur=Cur->Next);
    }

    #ifdef __LIST1_DEBUG__
	friend ostream& operator<<(ostream& Out, const Clist1_<T_>& List);
    #endif
};
// END OF CLASS List1_

#ifdef INCLUDE_TMPL_DEFS
#include "List1.c++"
#endif

// ==== END OF TEMPLATE HEADER List1.h ====

#endif	/* LIST1_TMPL_DECLS */
