#ifndef LIST1_TMPL_DEFS
#define LIST1_TMPL_DEFS

// ==== TEMPLATE DEFINITIONS List1.c++ ====

/* Unidirectional linked list templates. */

// SGI C++ 4.0, IRIX 5.3, 7. Mar. 1996. Andris Aszodi

// ==== Listhnd1_ METHODS ====

// ---- Assignment ----

/* = : Creates an almost exact copy of the rhs in the lhs. Note that the
 * Ptrno is unchanged! This assignment is used in the "deep copy"
 * methods of the List1_ class.
 */
template <class T_>
Listhnd1_<T_>& Listhnd1_<T_>::operator=(const Listhnd1_<T_>& Rhs)
{
    if (this==&Rhs) return(*this);  // x=x
    
    Item1_<T_> *Ocp, *Cp, *Prev;
    
    // copy items into existing storage in lhs
    for (Cp=Head, Ocp=Rhs.Head; Cp!=NULL && Ocp!=NULL; Prev=Cp, Cp=Cp->Next, Ocp=Ocp->Next)
	Cp->Value=Ocp->Value;
    
    if (Cp==NULL && Ocp!=NULL)	// rhs is longer, add its tail to lhs
	while (Ocp!=NULL)
	{
	    if (Tail==NULL)	// was empty
	    {
		Head=Tail=new Item1_<T_>(Ocp->Value);
	    }
	    else
	    {
		Tail->Next=new Item1_<T_>(Ocp->Value);
		Tail=Tail->Next;
	    }
	    Ocp=Ocp->Next;
	}
    else if (Cp!=NULL && Ocp==NULL)	// lhs is longer, truncate
    {
	for (Item1_<T_> *Del=Cp; Del!=NULL; Del=Cp)
	{
	    Cp=Del->Next; delete Del;
	}
	Tail=Prev; Tail->Next=NULL;
    }
    Len=Rhs.Len;
    return(*this);
}
// END of =

// ---- Memory management ----

/* make_copy(): constructs a copy of the calling list. The beginning
 * and end of the new list will be returned in Beg and End.
 * Return value: the length of the new list. Private
 */
template <class T_>
unsigned int Listhnd1_<T_>::make_copy(Item1_<T_>*& Beg, Item1_<T_>*& End) const
{
    if (!Len)	// empty list
    {
	Beg=End=NULL;
	return(0);
    }
    
    Item1_<T_> *Newitem, *Ocp;
    
    // process the old list
    Beg=NULL;
    Ocp=Head;	// point to first in original list
    do
    {
	Newitem=new Item1_<T_>(Ocp->Value); // create a new copy of old item
	if (Beg==NULL)	// add first item
	    Beg=End=Newitem;
	else
	    End->Next=Newitem; End=Newitem;	    // add to end of new list
	Ocp=Ocp->Next;  // move to next position
    }
    while (Ocp!=NULL);
    
    return(Len);
}
// END of make_copy()

/* remove_all(): removes all items from the calling object. Private */
template <class T_>
void Listhnd1_<T_>::remove_all()
{
    for (Item1_<T_> *Cp=Head; Cp!=NULL; Cp=Head)
    {
	Head=Cp->Next; delete Cp;
    }
    Len=0; Tail=NULL;
}
// END of remove_all()

// ==== Clist1_ METHODS ====

// ---- Assignment ----

/* = : "shallow" assignment. The Listhnd1_ pointed to by the lhs
 * is "released" and the lhs then will point to Rhs's list.
 */
template <class T_>
Clist1_<T_>& Clist1_<T_>::operator=(const Clist1_<T_>& Rhs)
{
    if (this==&Rhs) return(*this);	// x=x
    
    Lptr->Listhnd1_<T_>::~Listhnd1_();	// notify the list handler that one less will point to it
    Lptr=Rhs.Lptr;
    Lptr->Ptrno++;	// tell the Rhs's list that one more will point to it
    Cur=Rhs.Cur;
    return(*this);
}
// END of =

// ==== List1_ METHODS ====

// ---- Insertions ----

/* insert(Val): insert a new item into the calling object so that
 * the current item pointed to by Cur will follow it, i. e.
 * ...->Current->... will become ...->New->Current->... and Cur
 * will point to New. If Cur==NULL, then the new item will be
 * appended to the end of the list (cf. operator+= ).
 * Returns the calling object.
 */
template <class T_>
List1_<T_>& List1_<T_>::insert(const T_& Val)
{
    if (Cur==NULL)  // simply append to end
    {
	(*this)+=Val; Cur=Lptr->Tail;	// point to new item
	return(*this);
    }
    
    // create new container and move the current value into it
    Item1_<T_> *Newitem=new Item1_<T_>(Cur->Value);
    
    // include the new item into the list
    Newitem->Next=Cur->Next;
    Cur->Next=Newitem;
    
    // was it appended to the tail?
    if (Cur==Lptr->Tail) Lptr->Tail=Newitem;
    
    Cur->Value=Val; // new value goes into old current container
    (Lptr->Len)++;  // grow length
    return(*this);
}
// END of insert(Val)

/* insert(List): insert a whole List into the calling object
 * at the current position so that ...->Current->... (where Current
 * is pointed to by Cur) becomes ...->[List]->Current->...
 * and Cur will point to the first item in List. If Cur==NULL, 
 * then List will be appended to the end of the calling object.
 * Returns calling object.
 */
template <class T_>
List1_<T_>& List1_<T_>::insert(const List1_<T_>& List)
{
    if (!List) return(*this);	// empty List: don't bother
    
    Item1_<T_> *Newhead=NULL, *Newtail;
    unsigned int Llen=List.len();
    
    if (!len() || Cur==NULL)  // empty or current position is "outside"
    {
	List.Lptr->make_copy(Newhead, Newtail);
	if (!len()) Lptr->Head=Newhead; else Lptr->Tail->Next=Newhead;
	Cur=Newhead; Lptr->Tail=Newtail;
    }
    else    // general insertion at a valid Cur position
    {
	Item1_<T_> *Newitem, *Ocp=List.Lptr->Head;
	
	/* The first item from List goes into the container
	 * pointed to by Cur. Its previous value is saved
	 */
	Item1_<T_> Curval=Cur->Value;
	Cur->Value=Ocp->Value;
	Ocp=Ocp->Next;
	
	/* Create a "shifted-copy" list between Newhead and Newtail,
	 * the i:th item of List becomes the i-1:th
	 */
	while (Ocp!=NULL)
	{
	    Newitem=new Item1_<T_>(Ocp->Value);
	    if (Newhead==NULL) Newhead=Newitem;
	    else Newtail->Next=Newitem;
	    Newtail=Newitem;
	    Ocp=Ocp->Next;
	}
	
	/* Create the very last item which will contain Curval
	 * and append it to the new list
	 */
	Newtail->Next=new Item1_<T_>(Curval);
	Newtail=Newtail->Next;
	
	// build the new list into the old after Cur
	Newtail->Next=Cur->Next;
	Cur->Next=Newhead;
	
	// update Lptr->Tail if necessary
	if (Cur==Lptr->Tail) Lptr->Tail=Newtail;
    }
    Lptr->Len+=Llen;
    return(*this);
}
// END of insert(List)

/* List1+=Val : append a new item with the value Val at the end of the list.
 * List1+=List: append a copy of List at the end of the calling object.
 * Both return the calling object.
 */
template <class T_>
List1_<T_>& List1_<T_>::operator+=(const T_& Val)
{
    if (!len())	// was empty
	Lptr->Head=Lptr->Tail=Cur=new Item1_<T_>(Val);
    else
    {
	Lptr->Tail->Next=new Item1_<T_>(Val);
	Lptr->Tail=Lptr->Tail->Next;
    }
    (Lptr->Len)++;
    return(*this);
}

template <class T_>
List1_<T_>& List1_<T_>::operator+=(const List1_<T_>& List)
{
    if (!List) return(*this);	// empty, don't bother
    
    // create a copy of List
    Item1_<T_> *Newbeg, *Newend;
    unsigned int Llen=List.Lptr->make_copy(Newbeg, Newend);
    
    if (!len())	// was empty
	Lptr->Head=Cur=Newbeg;
    else
	Lptr->Tail->Next=Newbeg;
    Lptr->Tail=Newend;
    Lptr->Len+=Llen;
    return(*this);
}
// END of operator +=

/* List1^=Val: append a new item with value Val before the beginning of List1.
 * List1^=List: append List before the first item of List1.
 * The operator ^= was chosen because '^' means "the beginning of a string"
 * in some regular expressions (reminiscent of "beginning of a list).
 * Both return the calling object.
 */
template <class T_>
List1_<T_>& List1_<T_>::operator^=(const T_& Val)
{
    if (!len())	// was empty
	Lptr->Head=Lptr->Tail=Cur=new Item1_<T_>(Val);
    else
    {
	Item1_<T_> *Newitem=new Item1_<T_>(Val);
	Newitem->Next=Lptr->Head;
	Lptr->Head=Newitem;
    }
    (Lptr->Len)++;
    return(*this);
}

template <class T_>
List1_<T_>& List1_<T_>::operator^=(const List1_<T_>& List)
{
    if (!List) return(*this);	// empty, don't bother
    
    // create a copy of List
    Item1_<T_> *Newbeg, *Newend;
    unsigned int Llen=List.Lptr->make_copy(Newbeg, Newend);
    
    if (!len())	// was empty
    {
	Lptr->Tail=Newend; Cur=Newbeg;
    }
    else
	Newend->Next=Lptr->Head;
    Lptr->Head=Newbeg;
    Lptr->Len+=Llen;
    return(*this);
}
// END of operator ^=

// ---- Deletions ----

/* del(): delete N items (default 1) from the list, beginning with the
 * current item pointed to by Cur. If Cur==NULL or N==0, nothing happens. If N is
 * larger than the length of the tail of the list from Cur, then the
 * whole tail is deleted. The new current position will be the first
 * item after the deleted region.
 * Return value: the number of items actually deleted.
 */
template <class T_>
unsigned int List1_<T_>::del(unsigned int N)
{
    if (!N || Cur==NULL) return(0);	// do nothing
    
    Item1_<T_> *Prev, *Del;
    register unsigned int Dno;
    
    // set Prev to the item before Cur: clumsy
    if (Cur==Lptr->Head) Prev=NULL;
    else for (Prev=Lptr->Head; Prev!=NULL && Prev->Next!=Cur; Prev=Prev->Next);
    
    // delete items until end or N
    for (Dno=0, Del=Cur; Dno<N && Del!=NULL; Dno++, Del=Cur)
    {
	Cur=Del->Next; delete Del;
    }
    
    // reconnect
    if (Prev!=NULL) Prev->Next=Cur; else Lptr->Head=Cur;
    if (Cur==NULL) Lptr->Tail=Prev;
    (Lptr->Len)-=Dno;
    
    return(Dno);
}
// END of del()

// ==== END OF TEMPLATE DEFINITIONS List1.c++ ====

#endif	/* LIST1_TMPL_DEFS */
