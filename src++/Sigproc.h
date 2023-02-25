#ifndef SIGPROC_CLASS
#define SIGPROC_CLASS

// ==== PROJECT DRAGON: HEADER Sigproc.h ====

/* Signal trapping and multi-process management class. */

// SGI C++ 7.1, IRIX 6.2, 17. Oct. 1997. Andris Aszodi

/* NOTES TO DEVELOPERS:-
 * 
 * 1) Signal trapping is implemented using exception.
 * This will NOT work with older compilers. For SGI C++ 7.1, 
 * use the "-exceptions" flag when compiling for the O32 ABI.
 * 
 * 2) Use the System V or POSIX signalling mechanism for SGI.
 * Due to a bug in the <sys/wait.h> header, the BSD signalling
 * assumes that no typechecking is done which is not accepted
 * by the C++ compiler.
 */

// ---- DEFINITIONS ----

/* Three UNIX signalling mechanisms are supported: BSD, System V
 * and POSIX. System V is the default. Define the macro _BSD_SIGNALS
 * for BSD behaviour, and _POSIX_SIGNALS for POSIX behaviour.
 * When more macros are defined by mistake, the order of
 * precedence is SysV > POSIX > BSD.
 * WARNING: mixing these various flavours of signal handling
 * can lead to "unpredictable results"!
 */

#if !defined(_SYSV_SIGNALS) && !defined(_POSIX_SIGNALS) && !defined(_BSD_SIGNALS)
    #define _SYSV_SIGNALS   /* set default signal flavour */
#endif

// Linux signal handling: we use POSIX by all means
// SIG_PF must be defined
#ifdef __linux__
    #define _POSIX_SIGNALS  /* overrides BSD */
    #undef _SYSV_SIGNALS    /* this must be explicitly switched off */
    #define SIG_PF __sighandler_t
#endif

// SGI signal handling: we may use SYSV or POSIX but not BSD
#ifdef __sgi
#undef _BSD_SIGNALS
    #define _POSIX_SIGNALS  /* Sys V overrides if defined earlier */
#endif

#ifdef _SYSV_SIGNALS
    #ifdef _POSIX_SIGNALS
	#undef _POSIX_SIGNALS
	#undef _POSIX_C_SOURCE
    #endif
    #ifdef _BSD_SIGNALS
	#undef _BSD_SIGNALS
	#undef _BSD_COMPAT
    #endif
    #ifndef _SVR4_SOURCE
        #define _SVR4_SOURCE
    #endif
#else
    #undef _SVR4_SOURCE
#endif	/* _SYSV_SIGNALS */

#ifdef _POSIX_SIGNALS
    #ifdef _BSD_SIGNALS
	#undef _BSD_SIGNALS
	#undef _BSD_COMPAT
    #endif
    #ifndef _POSIX_C_SOURCE
        #define _POSIX_C_SOURCE 1
    #endif
#else
    #undef _POSIX_C_SOURCE
#endif  /* _POSIX_SIGNALS */

#ifdef _BSD_SIGNALS
    #define _BSD_COMPAT
#endif	/* _BSD_SIGNALS */

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <errno.h>
#include <iostream.h>
#include <iomanip.h>
#include <signal.h>

// 27-Mar-1997: SGI still has no std c++ library
#ifdef __sgi
    #include <exception.h>
#else
    #include <exception>
#endif
#include <unistd.h>

/* SuSE Linux 6.3, g++ 2.95.2:
 * pid_t seems to be undefined, very mysterious
 * 13-Feb-2000
 */
#ifdef __linux
    #ifndef pid_t
    typedef __pid_t pid_t;
    #define pid_t pid_t
    #endif
#endif

// ==== CLASSES ====

/* Sigexcept_: an exception object which stores
 * the signal number or 0 if no signal was caught.
 */
class Sigexcept_
{
    private:
    volatile int Sigval;	// SIGxxx value
    
    public:
    Sigexcept_(int Sig=0): Sigval(Sig) {} 
    int sigval() const { return(Sigval); }
};
// END OF CLASS Sigexcept_

/* Sigproc_: signal trapping and multiple process management class.
 * The signal handler function is a friend of this class. DRAGON
 * can spawn multiple copies of itself for parallel processing and
 * the details of the child processes are kept here.
 */
class Sigproc_
{
    // data
    private:
    
    // Status_: SINGLE (one process), PARENT, CHILD
    enum Status_ {SINGLE, PARENT, CHILD};
     
    pid_t Pid;	// process ID
    pid_t *Children;   // array of child process IDs (parent only)
    Status_ Stat;	    // what kind of process
    int Maxprocno, Maxchildno, 	// max. and actual child counts
	Runpp, Procno;	// no. of runs per process, current process no.
    volatile int Childno;   // can be changed by signal handler!
    
    // methods
    public:
    
	// constructor
    /* Sets up the object to manage Mprocno child processes. If Mprocno==0
     * (the default), then there will be just a SINGLE process; otherwise,
     * abs(Mprocno) processes will be spawned when required.
     */
    Sigproc_(int Mprocno=0):
	Pid(0), Children(NULL), Childno(0),
	Maxchildno(0), Runpp(0), Procno(0)
    { set_maxprocno(Mprocno); }

	// destructor
    ~Sigproc_() {delete [] Children; }
    
	// access
    /* set_maxprocno(): sets the no. of maximal allowed child processes
     * to Mprocno. If Mprocno==0, then the process status will be SINGLE, 
     * otherwise PARENT.
     * Return value: the maximal number of child processes.
     */
    int set_maxprocno(int Mprocno)
    {
	Maxprocno=abs(Mprocno); if (Maxprocno==1) Maxprocno=2;
	Stat=Maxprocno? PARENT: SINGLE;
	return(Maxprocno);
    }
    
    /* process type */
    int is_single() const { return(Stat==SINGLE); }
    int is_parent() const { return(Stat==PARENT); }
    int is_child() const { return(Stat==CHILD); }

	// signal trapping
    /* set_signal(): sets up the signal trap. Most signals which do not
     * cause a core dump will be trapped by signal_handler(). The fatal
     * signals will be caught and dealt with by the system.
     */
    void set_signal(SIG_PF handler_func);
    
    /* signal_handler(): catches signals, prints a message to cerr,
     * and throws a Sigexcept_ exception with the signal
     * as the return value. In multiprocess runs, SIGCHLD (child status
     * changed) is also caught and the child status is updated in
     * the external Sigproc object.
     */
    friend void signal_handler(int Sigtype);
    
	// multiple process management
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
    int spawn_children(int Runno);
    
    /* wait_4children(): when invoked in a parent process, this
     * method monitors the children. (When one of them finishes, 
     * signal_handler() takes a note.) If a signal other than
     * SIGCHLD is caught and if there are any active
     * children left, they are terminated. The method does
     * nothing if invoked in single or child processes.
     * Return value: the signal caught in the parent or 0.
     */
    int wait_4children();

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
    void get_runlimits(int Runno, int& Rcyclo, int& Rcychi) const;

    // forbidden methods
    private:
    
    Sigproc_(const Sigproc_&);
    Sigproc_& operator=(const Sigproc_&);
};
// END OF CLASS Sigproc_

// ---- GLOBALS ----

/* signal_message(): prints an explanatory message to "cerr"
 * for the signal Sigtype.
 */
void signal_message(int Sigtype);

// ==== END OF HEADER Sigproc.h ====
#endif	/* SIGPROC_CLASS */
