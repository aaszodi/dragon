// ==== PROJECT DRAGON: METHODS Sigproc.c++ ====

/* Signal trapping and multi-process management class. */

// SGI C++ 7.1, IRIX 6.5, 11-Nov-1999. Andris Aszodi

// ---- MODULE HEADER ----

#include "Sigproc.h"

// ---- STANDARD HEADERS ----

#include <string.h>
#include <sys/wait.h>

// ==== Sigproc_ METHODS ====

// ---- Signal trapping ----

/* set_signal(): sets up the signal trap. Most signals which do not
 * cause a core dump will be trapped by signal_handler(). Some fatal
 * signals will also be caught.
 */
void Sigproc_::set_signal(SIG_PF handler_func)
{
    // System V / BSD traps
    #if defined(_SYSV_SIGNALS) || defined(_BSD_SIGNALS)
	// non-fatals
	signal(SIGHUP, (SIG_PF)handler_func);	// hangup
	signal(SIGQUIT, (SIG_PF)handler_func);	// quit
	signal(SIGPIPE, (SIG_PF)handler_func);	// Broken pipe
	signal(SIGALRM, (SIG_PF)handler_func);	// software alarm
	signal(SIGTERM, (SIG_PF)handler_func);	// termination
	//fatals
	signal(SIGILL, (SIG_PF)handler_func);	// illegal instruction
	signal(SIGFPE, (SIG_PF)handler_func);	// floating-point exception
	signal(SIGBUS, (SIG_PF)handler_func);	// bus error
	signal(SIGSEGV, (SIG_PF)handler_func);	// segmentation fault
    #endif
    
    // POSIX traps
    #ifdef _POSIX_SIGNALS
	struct sigaction Sigact;
	
	Sigact.sa_handler=(SIG_PF)handler_func;	// install handler
	sigemptyset(&Sigact.sa_mask);	// clear signal block mask
	Sigact.sa_flags=0;		// clear signal flags
	// non-fatal
	sigaddset(&Sigact.sa_mask, SIGHUP);	// construct the mask
	sigaddset(&Sigact.sa_mask, SIGQUIT);	// of all signals caught
	sigaddset(&Sigact.sa_mask, SIGPIPE);	// by handler_func()
	sigaddset(&Sigact.sa_mask, SIGALRM);
	sigaddset(&Sigact.sa_mask, SIGTERM);
	// fatal
	sigaddset(&Sigact.sa_mask, SIGILL);
	sigaddset(&Sigact.sa_mask, SIGFPE);
	sigaddset(&Sigact.sa_mask, SIGBUS);
	sigaddset(&Sigact.sa_mask, SIGSEGV);
	
	// non-fatal
	sigaction(SIGHUP, &Sigact, NULL);   // specify action for each
	sigaction(SIGQUIT, &Sigact, NULL);
	sigaction(SIGPIPE, &Sigact, NULL);
	sigaction(SIGALRM, &Sigact, NULL);
	sigaction(SIGTERM, &Sigact, NULL);
	// fatal
	sigaction(SIGILL, &Sigact, NULL);
	sigaction(SIGFPE, &Sigact, NULL);
	sigaction(SIGBUS, &Sigact, NULL);
	sigaction(SIGSEGV, &Sigact, NULL);
    #endif
    
    // Ctrl-C is for adults only
    if (!is_child())
    {
	#if defined(_SYSV_SIGNALS) || defined(_BSD_SIGNALS)
	    signal(SIGINT, (SIG_PF)handler_func);
	#endif
	#ifdef _POSIX_SIGNALS
	    sigaddset(&Sigact.sa_mask, SIGINT);	// holds handler_func() already
	    sigaction(SIGINT, &Sigact, NULL);
	#endif
    }
    else    // children ignore Ctrl-C
    {
	#if defined(_SYSV_SIGNALS) || defined(_BSD_SIGNALS)
	    signal(SIGINT, SIG_IGN);
	#endif
	#ifdef _POSIX_SIGNALS
	    struct sigaction Ctrlcact;
	    
	    Ctrlcact.sa_handler=SIG_IGN;
	    sigemptyset(&Ctrlcact.sa_mask);
	    sigaddset(&Ctrlcact.sa_mask, SIGINT);
	    sigaction(SIGINT, &Ctrlcact, NULL);
	#endif
    }
    
    /* From Version 4.10 on, DRAGON can spawn separate processes
     * for parallel runs. Trapping SIGCHLD will take care of the
     * children; handler_func() implements the waiting cycle
     * for them. The logic follows the examples given in the 
     * IRIX wait(2) manpage.
     */
    if (is_parent())
    {
	#ifdef _SYSV_SIGNALS
	    sigset(SIGCHLD, (SIG_PF)handler_func);	// child status changed 
	    sighold(SIGCHLD);    // add to signal mask: block delivery 
	#endif
	#ifdef _POSIX_SIGNALS
	    sigaddset(&Sigact.sa_mask, SIGCHLD);	// holds handler_func() already
	    sigaction(SIGCHLD, &Sigact, NULL);
	#endif
	#ifdef _BSD_SIGNALS
	    struct sigvec Sigvec;
	    
	    Sigvec.sv_handler=(SIG_PF)handler_func;  // install handler 
	    Sigvec.sv_mask=sigmask(SIGCHLD);    // add to mask 
	    Sigvec.sv_flags=0;
	    sigvec(SIGCHLD, &Sigvec, NULL);
	    sigsetmask(sigmask(SIGCHLD));
	#endif
    }
}
// END of set_signal() 

/* signal_handler(): catches signals, prints a message to cerr,
 * and throws a Sigexcept_ exception with the signal
 * as the return value. In multiprocess runs, SIGCHLD (child status
 * changed) is also caught and the child status is updated in
 * the external Sigproc object.
 */
void signal_handler(int Sigtype)
{
    extern Sigproc_ Sigproc;
    
    cerr<<flush;

    /* If this is a parent process, then the handler waits for a child
     * process to return and when this happens, a notification
     * is written to cerr. No exception is thrown.
     */
    if (Sigproc.is_parent() && Sigtype==SIGCHLD)
    {
	int Pid, Status;
	
	#ifdef _SYSV_SIGNALS
	    Pid=wait(&Status);
	    Sigproc.Childno--;
	    cerr<<"\nCHILD PROCESS "<<Pid<<" exited with code "
		<<WEXITSTATUS(Status)<<endl;
	    sigset(SIGCHLD, (SIG_PF)signal_handler);	// reinstall signal
	#endif
	#ifdef _POSIX_SIGNALS
	    while((Pid=waitpid(-1, &Status, WNOHANG))>0)
	    {
		Sigproc.Childno--;
		cerr<<"\nCHILD PROCESS "<<Pid<<" exited with code "
		    <<WEXITSTATUS(Status)<<endl;
	    }
	#endif
	#ifdef _BSD_SIGNALS
	    // NOTE: this won't compile under IRIX 6.2
	    // because the prototype is int wait3()
	    // and typechecking is suspended
	    // which is OK in C but illegal in C++
	    while((Pid=wait3(&Status, WNOHANG, NULL))>0)
	    {
		Sigproc.Childno--;
		cerr<<"\nCHILD PROCESS "<<Pid<<" exited with code "
		    <<WEXITSTATUS(Status)<<endl;
	    }
	#endif
	return;
    }
    
    // print a message and reset default traps
    signal_message(Sigtype);
    Sigproc.set_signal(SIG_DFL);
    
    // Under Linux, compiled with GCC/G++ 2.8 and up to 2.95 at least,
    // the exception mechanism is not entirely perfect: 
    // on receiving Ctrl-C, the program aborts.
    // Let's ignore SIGINT-s and emit a warning
#if defined (__linux) && defined (__GNUC__) && (__GNUC__ <=2) && (__GNUC_MINOR__ <=95)
    if (Sigtype==SIGINT)
    {
	cerr<<"\n! signal_handler(): Cannot throw exception for SIGINT, continue\n";
	cerr<<"(OS is Linux, GCC version "<<__GNUC__<<"."
		<<__GNUC_MINOR__<<".x)\n";
	return;
    }
#endif
    throw(Sigexcept_(Sigtype));
}
// END of signal_handler()

/* signal_message(): prints an explanatory message to "cerr"
 * for the signal Sigtype.
 */
void signal_message(int Sigtype)
{
    cerr<<flush;
    if (Sigtype==SIGINT)
	cerr<<"User interrupt requested (Ctrl-C)\n";
    else
    {
	cerr<<"\nWARNING: Exiting on signal "<<Sigtype;
	switch(Sigtype)
	{
	    case SIGILL:	// fatals
	    case SIGFPE:
	    case SIGBUS:
	    case SIGSEGV:
		cerr<<" -- HORRIBLE ERROR CAUSED BY BUGGY CODE\n";
		abort();
	    break;
	    
	    // non-fatals
	    case SIGHUP:
	    cerr<<" -- hangup\n"; break;
	    case SIGQUIT:
	    cerr<<" -- quit\n"; break;
	    case SIGPIPE:
	    cerr<<" -- broken pipe\n"; break;
	    case SIGALRM:
	    cerr<<" -- software alarm\n"; break;
	    case SIGTERM:
	    cerr<<" -- termination\n"; break;
	    default:
	    cerr<<" -- unknown signal\n"; break;
	}
    }
    cerr<<flush;
}
// END of signal_message() 

// ---- Multiple process management ----

/* spawn_children(): attempts to spawn child processes
 * to run the Runno simulations among each other. The
 * calling object will spawn Mproc processes (should be set
 * beforehand) but only if Runno>=2. The actual number of
 * processes may be less (depending on system limits etc.).
 * The child processes will divide the runs among themselves
 * (almost) equally, the parent process will just monitor them
 * and won't carry out any simulations itself.
 * Return value: the actual number of children spawned in the
 * parent process, 0 otherwise.
 */
int Sigproc_::spawn_children(int Runno)
{
    // check if a previous run has switched a parent to single
    if (Maxprocno && is_single()) Stat=PARENT;
    
    if (!is_parent()) return(0);    // children and singles don't breed
    if (Runno<=1)
    {
	cerr<<"\n? Sigproc_::spawn_children("<<Runno<<"): Too few runs\n";
	Stat=SINGLE;	// switch back to serial mode
	return(0);
    }
    
    // prepare for launch
    set_signal((SIG_PF)signal_handler);	// set up signal trap
    
    // get no. of serial runs per process: Maxprocno is the maximum
    if (Runno<Maxprocno)
    { Maxchildno=Runno; Runpp=1; }  // one run per process
    else
    { Runpp=Runno/Maxprocno; Maxchildno=Maxprocno; }
    Childno=0;	// actual number of living children
    
    // storage for children's PIDs
    delete [] Children;
    Children=new pid_t [Maxchildno];
    
    // spawn child processes
    pid_t Forkval;
    errno=0;
    for (Procno=0; Procno<Maxchildno; Procno++)
    {
	Forkval=fork();	// spawn
	if (errno)	// could not fork
	{
	    cerr<<"\n! Sigproc_::spawn_children(): "<<strerror(errno)<<endl;
	    Childno=Maxchildno=Procno;
	    errno=0; break;
	}
	
	if (Forkval)    // I am the parent
	{
	    Stat=PARENT;
	    Children[Childno]=Forkval;  // store child PID
	    Childno++; 
	    cout<<"PROCESS #"<<(Procno+1)<<" (PID="<<Forkval<<") started\n"<<flush;
	}
	else	// I am the child
	{
	    Stat=CHILD;
	    break; // children mustn't fork
	}
    }
    return(Stat==PARENT? Childno: 0);
}
// END of spawn_children()

/* wait_4children(): when invoked in a parent process, this
 * method monitors the children. (When one of them finishes, 
 * signal_handler() takes a note.) If a signal other than
 * SIGCHLD is caught and if there are any active
 * children left, they are terminated. The method does
 * nothing if invoked in single or child processes.
 * Return value: the signal caught in the parent or 0.
 */
int Sigproc_::wait_4children()
{
    if (!is_parent()) return(0);    // don't do anything

    // wait for children to finish
    try
    {
	#ifdef _POSIX_SIGNALS
	    sigset_t Emptyset;  // needed for sigsuspend()
	    sigemptyset(&Emptyset);
	#endif
	
	while(Childno>0)
	{
	    #ifdef _SYSV_SIGNALS
		sigpause(SIGCHLD);
		sighold(SIGCHLD);
	    #endif
	    #ifdef _POSIX_SIGNALS
		sigsuspend(&Emptyset);
	    #endif
	    #ifdef _BSD_SIGNALS
		sigpause(0);
	    #endif
	}
	set_signal(SIG_DFL);
    }
    catch(Sigexcept_ Sigexc)	    // parent interrupted or killed: kill the family
    {
	if (Childno && Children!=NULL)
	{
	    // terminate child processes
	    for (unsigned int Ch=0; Ch<Maxchildno; Ch++)
	    {
		// don't even attempt to kill processes 0 or 1
		if (Children[Ch]>=2 && !kill(Children[Ch], SIGKILL))
		    cout<<"PROCESS (PID="<<Children[Ch]<<") killed\n";
	    }
	    delete [] Children;
	    Children=NULL; Childno=0;
	}
	set_signal(SIG_DFL);	// don't catch anything any more
	return(Sigexc.sigval());
    }
    return(0);	// normal termination
}
// END of wait_4children()

/* get_runlimits(): when multiple processes are spawned, the
 * Runno parallel simulation runs are divided among them
 * more or less equally. For each child process, the number
 * of the first and last run to be done is returned in Rcyclo
 * and Rcychi, respectively. (Eg for 3 children and Runno==15, 
 * the first child will do runs #1..#5, the second #6..#10, the
 * third #11..#15) In this way the results don't overwrite each
 * other. For single processes, Rcyclo==1, Rcychi==Runno.
 * If invoked in a parent process, the method does nothing.
 */
void Sigproc_::get_runlimits(int Runno, int& Rcyclo, int& Rcychi) const
{
    if (is_parent()) return;
    if (is_child())
    {
	Rcyclo=Procno*Runpp+1;
	Rcychi=(Procno==Maxchildno-1)? Runno: (Procno+1)*Runpp;
    }
    else { Rcyclo=1; Rcychi=Runno; }	// single process
}
// END of get_runlimits()

// ==== MATH LIBRARY EXCEPTION HANDLING ====

// This is unfortunately not implemented on every architecture.
// Works with IRIX 6.5 and the MipsPro 7.3 compilers,
// but not under Linux or Solaris.
// Define HAS_MATHERR for those which can do this trick.
#ifdef HAS_MATHERR

#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <float.h>

// SGI-specific definitions. <math.h> must come after _SGIAPI is defined.
#ifdef __sgi
    #include <ieeefp.h>
    #ifndef _SGIAPI
	#define _SGIAPI
    #endif
    
#endif

#include <math.h>

/* matherr(): under System V, this function, when defined,
 * is invoked at certain floating-point exceptions. There is
 * no way at present to intercept divisions-by-zero, but the
 * result shows up sooner or later as a NaN.
 * DOMAIN errors are trapped and the program is aborted unless
 * the error occurred because of a negative argument to sqrt()
 * or sqrtf() in which case the absolute value of the argument
 * is used. UNDERFLOW and OVERFLOW are also trapped and the return
 * values are set to 0.0 or DBL_MAX, respectively. All other errors
 * are treated by the default mechanism in the library.
 */
int matherr(register struct math_exception *Ex)
{
    switch(Ex->type)
    {
	case DOMAIN:
	if (!finite(Ex->arg1) || !finite(Ex->arg2))	/* abort on NaN-s */
	{
	    fprintf(stderr, "\n! DOMAIN NaN/inf argument(s) for %s, dying\n", Ex->name);
	    kill(getpid(), SIGFPE); /* kill myself */
	}
	else if (!strncmp(Ex->name, "sqrt", 4))  /* replace arg with abs val */
	{
	    fprintf(stderr, "\n? DOMAIN fp exception: %s(%f), abs val used\n", Ex->name, Ex->arg1);
	    Ex->retval=sqrt(-(Ex->arg1));
	}
	else
	{
	    fprintf(stderr, "\n! DOMAIN fp exception: op=%s, arg1=%f, arg2=%f, dying\n", 
		Ex->name, Ex->arg1, Ex->arg2);
	    kill(getpid(), SIGFPE);
	}
	return(1);  /* no default handling */
	
	case UNDERFLOW:
	fprintf(stderr, "\n? UNDERFLOW fp exception: op=%s, arg1=%f, arg2=%f, result 0.0\n", 
	    Ex->name, Ex->arg1, Ex->arg2);
	Ex->retval=0.0;
	return(1);

	case OVERFLOW:
	fprintf(stderr, "\n? OVERFLOW fp exception: op=%s, arg1=%f, arg2=%f, result DBL_MAX\n", 
	    Ex->name, Ex->arg1, Ex->arg2);
	Ex->retval=DBL_MAX;
	return(1);
	
	default: return(0); /* handle as usual */
    }
}
/* END of matherr() */
#endif	/* HAS_MATHERR */

// ==== END OF METHODS Sigproc.c++ ====
