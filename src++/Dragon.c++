// ==== PROJECT DRAGON: Dragon.c++ ====

/* Fold prediction using hierarchic distance matrix projection. */

// 5-Aug-2000. Andras Aszodi

/*******************************************/
/*                     ,     ,             */
/*     Software by Andras Aszodi, PhD      */
/*          All rights reserved.           */
/*        Alle Rechte vorbehalten.         */
/*         Minden jog fenntartva.          */
/*                                         */
/*           Previous address:             */
/*    Division of Mathematical Biology     */
/* National Institute for Medical Research */
/* The Ridgeway, Mill Hill, London NW7 1AA */
/*            UNITED KINGDOM               */
/*  e-mail: aaszodi@nimr.mrc.ac.uk         */
/*                                         */
/*            Present address:             */
/*    Novartis Forschungsinstitut GmbH     */
/*    Brunnerstrasse 59, A-1235 Vienna     */
/*                AUSTRIA                  */
/*                                         */
/* Tel: (++43 1) 866 34, ext. 9052         */
/* Fax: (++43 1) 866 34 727                */
/* Andras.Aszodi@pharma.novartis.com       */
/*                                         */
/*******************************************/

// ---- VERSION HISTORY ----

// 4.0: 6-Jun-95
// 4.1: 20-Jun-95
// 4.2: 23-Jun-95
// 4.3: 3-Jul-95
// 4.4: 4-Jul-95
// 4.5: 7-Jul-95
// 4.6: 11-Aug-95
// 4.7: 20-Oct-95
// 4.8: 27-Nov-95
// 4.9: 18-Dec-95
// 4.10: 4-Jan-96
// 4.11: 15-Jan-96
// 4.12: 26-Feb-96
// 4.13: 7-May-96
// 4.14: 15-May-96
// 4.15: 21-May-96
// 4.16: 10-July-96
// 4.17: 16-Jan-97
// 4.18: March 1998 (?)

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <math.h>

// ---- C++ MODULE HEADERS ----

#include "Access.h"
#include "Clip.h"
#include "Density.h"
#include "Restr.h"
#include "Hmom.h"
#include "Homodel.h"
#include "Iproj.h"
#include "Params.h"
#ifdef USE_PVM
    #include "Pvmtask.h"
#endif 
#include "Pieces.h"
#include "Polymer.h"
#include "Output.h"
#include "Score.h"
#include "Sigproc.h"
#include "Score.h"
#include "Steric.h"
#include "Sterchem.h"
#include "Tangles.h"
#include "Viol.h"

#ifdef USE_OPENGL_GRAPHICS
    #include "Graphics.h"
#endif

// ---- C MODULE HEADERS ----

#include "version.h"

// ---- C++ UTILITY HEADERS ----

#include "String.h"

// ---- C UTILITY HEADERS ----

#include "cmdopt.h"
#include "tstamp.h"

// ---- GLOBAL VARIABLES ----

Sigproc_ Sigproc;   // signal trapping and multiprocess management

#ifdef USE_PVM
Pvmtask_ Pvmtask;   // Parallel Virtual Machine support
#endif

// ---- STATIC VARIABLES ----

static Params_ Params;
static Polymer_ Polymer;
static Restraints_ Restraints;
static Homodel_ Homodel(Polymer);	// bind to Polymer
static Access_ Access;
static Steric_ Steric;

static unsigned int Rno=10;	    // just a non-0 value
static Pieces_ Pieces(Rno);	    // must have a ctor

// ---- PROTOTYPES ----

#ifdef USE_PVM
unsigned int master_pvmrun(unsigned int Jobno=1);
static void init_pvmslave();
#endif

unsigned int dragon_run(unsigned int Runno=1);
static void init_dragon();
static void merge_distmat(const Trimat_& Bestdist, Trimat_& Dist);

// ==== MAIN ====

int main(int argc,  char *argv[])
{
    /* The program understands the following switches:-
     * -p param_file: reads the parameters from param_file (3.x behaviour)
     * and runs only once
     * -p param_file -r run_no: like above but runs run_no times
     * -c command_file: interprets commands from command_file
     * Neither -p nor -c: interactive mode (cf. "Clip" module)
     * -h: prints a short help
     * -m procno: spawns procno processes for parallel runs (min. 2)
     * -M: spawns a slave task on every node in the PVM if available
     * -A: give The Answer and exit
     * The options are processed by the "cmdopt" module.
     */
    parse_optstr("hA c%s<command_file> m%d<process_no> M p%s<param_file> r%d<run_no>");
    if (get_options(argc, argv)<0 || optval_bool('h'))
    {
	char *Help=opt_helpstr();   // generate help string
	cerr<<flush<<version_string()<<endl;
	cerr<<"\nUsage: "<<argv[0]<<" "<<Help<<endl;
	cerr<<"Options:-\n";
	cerr<<"No options: run in interactive mode (press \'h\' for help)\n";
	cerr<<"-c <command_file>: execute commands from <command_file>\n";
	cerr<<"-h: print this help and exit\n";
	cerr<<"-m <process_no>: spawn <process_no> processes (>=2) for parallel runs\n";
#ifdef USE_PVM
	cerr<<"-M: spawn a slave on every node in the PVM\n";
#endif
	cerr<<"-p <param_file>: perform one run with parameters in <param_file>\n";
	cerr<<"-p <param_file> -r <run_no>: perform <run_no> runs with parameters in <param_file>\n";
	cerr<<"-A: give The Answer and exit\n";
	free(Help);
	return(EXIT_FAILURE);
    }
    
    // greet the user
    cout<<endl<<greeting()<<"!\n";
    cout<<"Welcome to "<<version_string()<<endl;
    cout<<"                                      ,     ,\n";
    cout<<"Algorithms by William R. Taylor & Andras Aszodi\n";
    cout<<"                      ,     ,\n";
    cout<<"Implementation by Andras Aszodi\n";
    cout<<"(C) 1993-2000. All rights reserved.\n\n";
#if defined(__sun)
    cout<<"SUN Solaris port by Nigel W. Douglas\n";
#endif
#if defined(USE_PVM) && (defined(__sgi) || defined(__sun) || defined(__linux))
    cout<<"MP support under PVM by J. Hungershoefer\n";
#endif
    
    // if you ported DRAGON to another architecture,
    // you may add your name here
#if defined(__some_architecture__)
    cout<<"SOME MACHINE port by YOUR NAME HERE\n";
#endif
        
    // list optional components
    cout<<"PVM: ";
    #ifdef USE_PVM
	cout<<"supported\n";
    #else
	cout<<"not supported\n";
    #endif
    cout<<"OpenGL graphics: ";
    #ifdef USE_OPENGL_GRAPHICS
	cout<<"supported\n";
    #else
	cout<<"not supported\n";
    #endif
    
    char *Parfnm=NULL;	// parameter file name
    int Dretval=0;    // DRAGON return value (0 is OK)
    
    // this option is not entirely serious :-)
    if (optval_bool('A'))
    {
      cerr<<"The Answer is 42.\n";
      exit(42);
    }

    /* The -p <param_file> option just reads the <param_file> and
     * then performs a single run or multiple runs if -r <run_no> was specified.
     * Overrides -c, -m, -M
     */
    if (optval_str('p', &Parfnm))
    {
	cout<<"# Reading parameters from file \""<<Parfnm<<"\"\n";
	if (!Params.read_file(Parfnm))
	    cerr<<"\n? Using default parameters\n";
    }
    
    // check if runs were requested on the command line
    int Runno=0;
    if (optval_int('r', &Runno))
    {
	Runno=abs(Runno);	// silently convert negative values
	if (!Runno)
	{
	    cerr<<"\n? Zero runs requested, exiting...\n";
	    return(0);
	}
    }	    // here Runno==0 means that -r was not seen
    
    // Otherwise, the command-line interpreter runs the show.
    Clip_ Clip(Params, "DRAGON");	// command-line interpreter
    char *Cmdfnm=NULL;	// command script file name (from -c option)
    
    // if a command script was specified with the -c option, then use it
    // provided no -r option was specified
    if (!Runno) optval_str('c', &Cmdfnm);
    
    // enable PVM? If yes, overrides -m
#ifdef USE_PVM
    if (optval_bool('M'))
    {
	Pvmtask.enrol_pvm("dragon");
	
	// PVM master
	if (Pvmtask.is_master())
	{
	    cout<<"PVM enabled.\n";
	    if (Runno)	// -r was specified
	    {
		Dretval=master_pvmrun(Runno);	// run immediately
		return Dretval;     	    	// then exit
	    }
	    else    // run under the command interpreter's control
	    	Dretval=Clip.get_command(Cmdfnm, master_pvmrun);
	}
	
	// PVM slave
	if (Pvmtask.is_slave())
	{
	    int Result=0, Tag, Plen;
            // determine number of cpus on this node
            // and send information to master 
            (void)Pvmtask.send_ncpus();
            // make process a nice one (ignore failure)
            // (fortunately the call is the same on Solaris, Linux and Irix)
            Result=nice(10);
        
	    while (1)
	    {
		// wait for a master message to arrive
		Tag=Pvmtask_::ANY;
		Result=Pvmtask.wait_master(Tag);
		if (Result<0) break;
		
		// parameters have been sent
		if (Tag==Pvmtask_::PARAMS)
		{
		    Plen=Pvmtask.recv_params(Params);
		    cout<<"Parameters received: Plen="<<Plen<<endl;
		    if (Plen<0) break;
		    
		    // some of the parameters have been changed, update
		    if (Plen>0) init_pvmslave();
		    
		    // tell master that I'm ready to run
		    Pvmtask.slave_ready();
		    continue;
		}
		
		// a run has been requested
		if (Tag==Pvmtask_::RUN)
		{
		    int Jobno=0;
		    Jobno=Pvmtask.recv_job();	// Jobno>0 is the job number
		    if (Jobno<0) break;	// error
		    else cout<<"Received job #"<<Jobno<<endl;
		    
		    /* Do the Jobno-th simulation. dragon_run()
		     * will know from Pvmtask that Jobno should be
		     * interpreted as the no. of the job (not the
		     * max. number of simulations)
		     */
		    Dretval=dragon_run(Jobno);
		    
		    // exit if a signal was caught during calc
		    if (Dretval)
		    {
			cout<<"Signal "<<Dretval<<" caught after dragon_run():"<<Pvmtask.id_str()<<endl;
			break;
		    }
		    
		    // report that the job is done
		    Pvmtask.job_status(Pvmtask_::SLAVE_DONE, Jobno);
		    Pvmtask.slave_ready();
		}
	    }	    // while(1)
	    // errors, signals and master's death land us here
	}	// if (PVM slave)
    }
    #endif
    
#ifdef USE_PVM
    // non-PVM run
    if (Pvmtask.no_pvm())
    {
#endif
	// set the multiple process management object
	int Mproc=0;
	optval_int('m', &Mproc);
	if (Mproc=Sigproc.set_maxprocno(Mproc))	// = intended
	    cout<<Mproc<<" parallel processes enabled.\n";

	if (Runno)
	    Dretval=dragon_run(Runno); // run Runno times
	else
	    Dretval=Clip.get_command(Cmdfnm, dragon_run);
#ifdef USE_PVM
    }
#endif
    cout<<"\nThank you for using DRAGON. Goodbye.\n";
    
    return((Dretval>0 && Dretval!=SIGINT)? Dretval: EXIT_SUCCESS);
}

// ==== FUNCTIONS ====

// ---- PVM support ----

#ifdef USE_PVM
/* master_pvmrun(): sends the parameter values to the PVM slaves
 * and requests Jobno simulations to be done.
 * Catches signals and passes them to the slaves via
 * the Sigproc mechanism. Used by Clip_ in place of dragon_run().
 * Return value: 0 upon normal termination.
 */
unsigned int master_pvmrun(unsigned int Jobno)
{
    int Signal=0;
    
    Pvmtask.spawn_slaves();	// start up slaves
    Sigproc.set_signal((SIG_PF)signal_pvm);    // set the PVM master signal trap

    // send parameters and do the runs
    if (Pvmtask.send_params(Params)>=0)
	Signal=Pvmtask.send_jobs(Params, Jobno);
    else Jobno=0;
    
    Sigproc.set_signal(SIG_DFL);    // unset the trap
    cout<<Jobno<<" job"<<(Jobno==1? "":"s")<<" done."<<endl;
    return(Signal);
}
// END of master_pvmrun()

/* init_pvmslave(): initialises the static global variables
 * before a simulation in a slave running under PVM. The
 * main difference is that the contents of the data files arrive
 * as strings.
 */
static void init_pvmslave()
{
    static bool Chainchg=true;    // has the chain changed?
    char *Fstr=NULL;	// points to static data buffer
    
    // initialise the polymer
    if (Params.changed("Alnfnm") || Params.changed("Masterno"))
    {
	Fstr=Pvmtask.recv_filestr(Pvmtask_::ALN);
	Polymer.str_aln(Fstr, Params.i_value("Masterno"));
	Chainchg=true;
	Params.reset_changed("Alnfnm");	// touch this as well
    }
    if (!Polymer.len())
    {
	cerr<<"\n! No valid polymer chain, exiting...\n";
	exit(0);
    }
    Rno=Polymer.len();
     
    if (Params.changed("Phobfnm"))
    {
	Fstr=Pvmtask.recv_filestr(Pvmtask_::PHO);
	Polymer.str_phob(Fstr);
	Params.reset_changed("Phobfnm");
    }
    if (Params.changed("Volfnm"))
    {
	Fstr=Pvmtask.recv_filestr(Pvmtask_::VOL);
	Polymer.str_vol(Fstr);
	Params.reset_changed("Volfnm");
    }
    if (Params.changed("Simfnm"))
    {
	Fstr=Pvmtask.recv_filestr(Pvmtask_::SIM);
	Polymer.str_simil(Fstr);
	Params.reset_changed("Simfnm");
    }
    if (Params.changed("Adistfnm"))
    {
	Fstr=Pvmtask.recv_filestr(Pvmtask_::ACD);
	Polymer.str_acdist(Fstr);
	Params.reset_changed("Adistfnm");
    }
    
    cout<<"\n=== THE MODEL CHAIN ===\n\n"<<Polymer;
    
    // initialise the distance limits
    Restraints.set_size(Rno);
    if (Chainchg || Params.changed("Restrfnm") || 
	    Params.changed("Homfnm") || Params.changed("Maxdist") ||
	    Params.changed("Minsepar"))
    {
	Fstr=Pvmtask.recv_filestr(Pvmtask_::RESTR);
	istrstream Ifs(Fstr); Ifs>>Restraints;
	Restraints.convert_restraints(Polymer);
	Fstr=Pvmtask.recv_filestr(Pvmtask_::HOM);
	int Kstr=Homodel.str_readknown(Fstr);
	if (Kstr>0)
	    Restraints.add_restrs(Homodel.make_restrs(Params.f_value("Maxdist"), Params.i_value("Minsepar")));
	else
	    cout<<"<No homology-derived distance restraints>\n";
	Params.reset_changed("Restrfnm");
	Params.reset_changed("Homfnm");
	Params.reset_changed("Maxdist");
	Params.reset_changed("Minsepar");
    }
    cout<<"\n=== DISTANCE LIMITS ===\n\n"<<Restraints
	<<"Total number of restraints: "<<Restraints.restr_no()<<endl;

    // initialise accessibility
    Access.set_size(Rno);
    if (Chainchg || Params.changed("Accfnm"))
    {
	Fstr=Pvmtask.recv_filestr(Pvmtask_::ACC);
	istrstream Ifs(Fstr); Ifs>>Access;
	Params.reset_changed("Accfnm");
    }
    cout<<"\n=== KNOWN ACCESSIBILITIES ===\n\n"<<Access;
    
    // get the secondary structure pieces (if any)
    if (Chainchg)	    // chain size changed
    {
	Pieces.res_no(Rno);	// resize and always read the secstr specs
	Fstr=Pvmtask.recv_filestr(Pvmtask_::SSTR);
	istrstream Ifs(Fstr); Ifs>>Pieces;
	Params.reset_changed("Sstrfnm");
	Chainchg=false;
    }
    else if (Params.changed("Sstrfnm"))
    {
	Fstr=Pvmtask.recv_filestr(Pvmtask_::SSTR);
	istrstream Ifs(Fstr); Ifs>>Pieces;
	Params.reset_changed("Sstrfnm");
    }
    Restraints.setup_restr(Pieces, Polymer);
    Steric.setup(Rno);
    
    cout<<"\n=== SECONDARY STRUCTURE ===\n\n"<<Pieces;
}
// END of init_pvmslave()
#endif	/* USE_PVM */

// ---- Simulation ----

/* init_dragon(): initialises the static global variables
 * before a simulation. The objects are updated only if a
 * relevant parameter has changed.
 * Note: see init_pvmslave() which does the same thing for
 * PVM runs.
 */
static void init_dragon()
{
    static bool Chainchg=true;	// has the chain changed?
    
    // initialise the polymer
    if (Params.changed("Alnfnm") || Params.changed("Masterno"))
    {
	Polymer.read_aln(Params.s_value("Alnfnm"), Params.i_value("Masterno"));
	Chainchg=true;
    }
    if (!Polymer.len())
    {
	cerr<<"\n! No valid polymer chain, exiting...\n";
	exit(0);
    }
    Rno=Polymer.len();
    
    if (Params.changed("Phobfnm"))
	Polymer.read_phob(Params.s_value("Phobfnm"));
    if (Params.changed("Volfnm"))
	Polymer.read_vol(Params.s_value("Volfnm"));
    if (Params.changed("Simfnm"))
	Polymer.read_simil(Params.s_value("Simfnm"));
    if (Params.changed("Adistfnm"))
	Polymer.read_acdist(Params.s_value("Adistfnm"));
    
    cout<<"\n=== THE MODEL CHAIN ===\n\n"<<Polymer;
    
    // initialise the distance limits
    Restraints.set_size(Rno);
    if (Chainchg || Params.changed("Restrfnm") || 
	    Params.changed("Homfnm") || Params.changed("Maxdist") ||
	    Params.changed("Minsepar"))
    {
	Restraints.read_restrs(Params.s_value("Restrfnm"), Polymer);
	int Kstr=Homodel.read_knownstr(Params.s_value("Homfnm"));
	if (Kstr>0)
	    Restraints.add_restrs(Homodel.make_restrs(Params.f_value("Maxdist"), Params.i_value("Minsepar")));
	else
	    cout<<"<No homology-derived distance restraints>\n";
    }
    
    // do not list more than 200 restraints to a terminal
    cout<<"\n=== DISTANCE LIMITS ===\n\n";
    if (Restraints.restr_no()>200 && isatty(STDOUT_FILENO))
	cout<<"More than 200 restraints, listing suppressed\n";
    else cout<<Restraints;
    cout<<"Total number of restraints: "<<Restraints.restr_no()<<endl;

    // initialise accessibility
    Access.set_size(Rno);
    if (Chainchg || Params.changed("Accfnm"))
	Access.read_file(Params.s_value("Accfnm"));
    cout<<"\n=== KNOWN ACCESSIBILITIES ===\n\n"<<Access;
    
    // get the secondary structure pieces (if any)
    if (Chainchg)	    // chain size changed
    {
	Pieces.res_no(Rno);	// resize and always read the secstr specs
	Pieces.read_secstr(Params.s_value("Sstrfnm"));
	Chainchg=false;
    }
    else
    {
	if (Params.changed("Sstrfnm"))
	    Pieces.read_secstr(Params.s_value("Sstrfnm"));
    }
    Restraints.setup_restr(Pieces, Polymer);
    Steric.setup(Rno);
    
    cout<<"\n=== SECONDARY STRUCTURE ===\n\n"<<Pieces;
}
// END of init_dragon()

/* dragon_run(): performs a full DRAGON simulation run Runno (default 1)
 * times using the parameters in Params. The objects inside are
 * updated only at the entry and only if the corresponding parameters
 * have changed. If the global variable Mproc>=2, then Mproc separate
 * processes are spawned which execute in parallel. Each of these
 * processes will run Runno/Mproc simulations sequentially.
 * If Mproc==1, then a separate process will be launched for each
 * of the Runno simulations. No process will be spawned if Runno==1.
 * Return value: 0 if OK, otherwise the value of a signal caught inside.
 * From Version 4.11 on, PVM support is also built in. The PVM status
 * is read from the static global Pvmtask object.
 */
unsigned int dragon_run(unsigned int Runno)
{
    // ---- Initialisation ----
    
#ifdef USE_PVM
    if(!Pvmtask.is_slave())  // task not under PVM
#endif
	init_dragon();  // non-PVM cases

    // set up projection
    Iproj_ Iproj(Rno+2);
    Iproj.set_size(Rno+2);
    Iproj.make_clusters();

    // set up detangling
    Tangles_ Tangles(Pieces);
    static const double TADJ=0.5;   // tangle adjustment scaling
    unsigned int Tangviol=0, Tangiter=Params.i_value("Tangiter");
    
    // set up the distance matrix and the coordinates: allow extra 2 points for N/C term
    Trimat_ Dista(Rno+2), Distbest(Rno+2);
    Fakebeta_ Fakebeta(Rno);	// puts extra 2 points there automagically
    Points_ Model(Rno+2, Rno), Best(Rno+2);
    
    // set up scores
    Scores_ Distsco(Params.f_value("Minscore"), Params.f_value("Minchange")), 
	Euclsco(Params.f_value("Minscore"), Params.f_value("Minchange")), 
	Bestsco(Params.f_value("Minscore"), Params.f_value("Minchange"));

    // init graphics if enabled
    int Graph=Params.i_value("Graph");
    #ifdef USE_OPENGL_GRAPHICS
	Graphics_ Draw;
	if (Graph) Draw.update_polymer(Polymer);
    #endif
    
    // set the floating-point output for cout and cerr
    long Coutf=cout.flags(), Cerrf=cerr.flags();
    int Oldoutprec=cout.precision(3), Olderrprec=cerr.precision(3);
    cout.precision(3); cout.setf(ios::scientific, ios::floatfield);
    cerr.precision(3); cerr.setf(ios::scientific, ios::floatfield);
    
    // ---- Main iteration cycle ----
    
    // codes for various exit types
    typedef enum {NOEXIT=0, EXIT_SIGNAL, EXIT_CTRLC, EXIT_SCOREOK,
	EXIT_MAXITER, EXIT_REPROJ} Exitreason_ ;
    Exitreason_ Exreason=NOEXIT;
    
    unsigned int Itno=0, It3dno=0, Dim=Rno+2, Oldim=Rno+2, 
	Bestfound=0, Reprojmax, Repriter, Reprojno, 
	Speciter=Params.i_value("Speciter");
    float Stress=0.0, Rmss=0.0, Densfact=0.0, Speceps=Params.f_value("Speceps");
    int Signal=0, Handflip=1, Logfd, Workdone=0, Noconv=0;
    String_ Outname, Logname;
    
    // set 3D reprojections
    Reprojmax=Params.i_value("Maxiter")/10+1;
    if (Reprojmax<3) Reprojmax=3;
    
    // init the scoring system
    Steric.reset_viol(Restraints, Rno+2, Distsco);
    Steric.reset_viol(Restraints, Rno+2, Euclsco);
    Steric.reset_viol(Restraints, Rno+2, Bestsco);
    
    /* Set up multiple process spawns. DRAGON can run in parallel
     * either if a -m flag requested that several copies be spawned
     * using fork() (cf. "Sigproc" module) on the same machine, 
     * or if a -M flag requested PVM support (no forking, copies
     * run on different machines). Under PVM, Sigproc.is_single()
     * is "true". If PVM wasn't linked in or no PVM run was requested, 
     * then Pvmtask.no_pvm() is "true".
     */
    
#ifdef USE_PVM
    if (Pvmtask.is_slave()) Graph=0;	// no graphics in PVM runs
    else
    {
#endif
	Sigproc.spawn_children(Runno);
	if (!Sigproc.is_single()) Graph=0;	// no graphics in multiprocess runs
#ifdef USE_PVM
    }
#endif

    // will be true only if no PVM run was requested
    if (Sigproc.is_parent())
	Signal=Sigproc.wait_4children();   // sit here until children finish
    else    // not parent
    {
	/* Run the simulations: either as children of a parent,
	 * or a simple serial run (w/o the "-m" option) or
	 * a multiple run which could not spawn children
	 * or a PVM slave
	 */
	int Rcyc, Rcyclo, Rcychi;
	
#ifdef USE_PVM
	/* for PVM runs, interpret Runno as the simulation limits */
	if (Pvmtask.is_slave())
	    Rcyclo=Rcychi=Runno;    
	else
#endif
	    Sigproc.get_runlimits(Runno, Rcyclo, Rcychi);

	cout<<"RUN from "<<Rcyclo<<" to "<<Rcychi<<endl;
	
	/* The main simulation cycle */
	for (Rcyc=Rcyclo; !Signal && Rcyc<=Rcychi; Rcyc++)
	{
	
	    /* Output redirection to a logfile is done
	     * in multiprocess runs only: for PVM runs, 
	     * each slave task maintains its own logfile
	     */
	    if (Sigproc.is_child())	// I am a child
	    {	
		// redirect cout and cerr to a logfile
		Logname=Params.s_value("Outfnm");
		make_outname(Logname, Rcyc, "log");
		errno=0;
		Logfd=open(Logname, O_CREAT|O_WRONLY, 0644);
		if (Logfd<0)
		{
		    cerr<<"Child logfile:"<<strerror(errno)<<endl;
		    errno=0;
		}
		else
		{
		    cout<<flush; cerr<<flush;
		    dup2(Logfd, STDOUT_FILENO);
		    dup2(Logfd, STDERR_FILENO);
		    cout<<"CHILD PROCESS ID="<<getpid()<<endl;
		}
		
	    }
	    
	    cout<<"\nRUN "<<Rcyc<<" STARTED: "<<time_stamp()<<endl;
	    start_timer();
	    
	    /* Initialise the distance matrix to random values
	     * within the pre-calculated bounds, modified by the
	     * hydrophobic distances for "soft" restraints.
	     * If the "Randseed" parameter is 0, then the time(NULL)
	     * value is used as the seed (this is fairly random).
	     * For parallel runs, this value is "spiced"
	     * with either the task ID (under PVM) or the child
	     * process ID (multiprocess run) to avoid identical
	     * seeds on time-sychronised machines.
	     */
	    long Randseed=Params.i_value("Randseed");
	    if (!Randseed || Runno>1)  // "random" start (always w/ multiples)
	    {
		Randseed=time(NULL);
	    #ifdef USE_PVM
		if (Pvmtask.is_slave())
		    Randseed+=Pvmtask.tid();
	    #endif
		if (Sigproc.is_child())
		    Randseed+=1024*getpid();	// give it a "large" perturbation
	    }
	    cout<<"# Randseed="<<Randseed<<endl;
	    Restraints.init_distmat(Dista, Polymer, Randseed);
    
	    Itno=It3dno=Repriter=Reprojno=0;
	    Oldim=Dim=Rno+2; Bestfound=0;
	    Rmss=0.0; Densfact=0.0;
	    Distsco.set_noexit();	// "prime" the scores
	    Euclsco.set_noexit();
	    Bestsco.set_noexit();
	    Exreason=NOEXIT;
	    
	    /* Signal traps: most non-fatal signals are trapped
	     * and then the handler throws a Sigexcept_
	     * exception.
	     */
	    Sigproc.set_signal((SIG_PF)signal_handler);	    // set up signal trap
	    do	// <--cycle until finished, interrupted or crashed
	    {
		try	// look for signal exceptions
		{
		    // print how much we have done
		    Workdone=int(100.0*((Dim==3)? Rno-1+It3dno: 3.0-Dim+Rno-1.0)/(Rno-1.0+Params.i_value("Maxiter")));
	    	    stop_timer();
		    cout<<"CYCLE: "<<(Itno+1)<<" ("
			<<Workdone<<"%, "
			<<time_string(timer_results(TS_UTIME|TS_STIME))
			<<")\n";
		#ifdef USE_PVM
		    // tell the master which cycle is being done
		    if (Pvmtask.is_slave())
			Pvmtask.job_status(Pvmtask_::SLAVE_RUNNING, Itno+1);
		#endif
		    // distance "space" adjustments in hyperspace:
		    if (Dim>3 || Repriter==Reprojmax)
		    {
			// in 3D, count how many reprojections were done
			if (Dim==3) Reprojno++;
			
			// make distance matrix from previous coords
			if (Itno) Model.dist_mat2(Dista);
    
			// adjust density
			Densfact=scale_distdens(Dista,
			    Restraints.exp_rad(Rno, Params.f_value("Density")));
    
			// "blend in" previous best distmat
			if (Dim==3 && Bestfound)
			    merge_distmat(Distbest, Dista);
			
			// ideal distances
			Fakebeta.update(Dista, Polymer);    // get C:beta-related distances
			Steric.ideal_dist(Dista, Fakebeta, Restraints, Polymer, 
				Pieces, Steric_::ALL | Steric_::RESTR | Steric_::SPECGRAD);
			Steric.adjust_dist(Dista, Pieces, Steric_::ALL);
			
			// display new distance matrix
			#ifdef USE_OPENGL_GRAPHICS
			    if (Graph) Draw.display_dist(Dista);
			#endif
			
			/* perform full projection: in 3.x compatibility mode
			 * the overall projection will automatically be used
			 * (cf. above the Iproj.set_size() call)
			 */
			Dim=Iproj.full_project(Dista, Params.f_value("Evfract"), Oldim, Model);
	    
			// post-projection refinement
			Densfact=proj_dens(Dista, Pieces, Model);
			Noconv=0;
			Stress=Steric.adjust_xyz(Model, Speciter, Speceps, Noconv);	// uses the pre-projection dists
			if (Noconv || Stress<0.0)	// on error or no convergence
			{
			    Model.dist_mat2(Dista);
			    Steric.adjust_xyz(Dista, Model, Pieces, Steric_::ALL);
			}
			
			/* See if the BOND distance scores improved compared to
			 * the previous embedding: if yes, then employ
			 * a bolder dim reduction strategy 
			 */
			Model.dist_mat2(Dista);
			Fakebeta.update(Dista, Polymer);    // get C:beta-related distances
			Steric.ideal_dist(Dista, Fakebeta, Restraints, Polymer, 
				Pieces, Steric_::ALL | Steric_::RESTR | Steric_::SCORE, &Distsco);
			if (Distsco[Scores_::BOND].change()<0)
			    Oldim=(Dim>4)? (2*(Dim-3))/3+3: 4;
			else
			    Oldim=(Dim>3)? Dim: 4;
			
			/* 3D-specific things: overall handedness adjustment,
			 * reset 3D reproj counter
			 */
			if (Dim==3)
			{
			    // get correct enantiomer from comparison to known homol. struct
			    if (Homodel.known_no()>0)
				Handflip=Homodel.hand_check(Model);
			    else	// we are on our own (check is based on secstr hands)
				Handflip=hand_check(Pieces, Model);
			    Repriter=0;
			}
			
			// get accessibility score and list
			Distsco[Scores_::ACCESS].score(Access.score_dist(Polymer, Dista));
			cout<<"DIST: "<<Distsco<<endl;
			cout<<"PROJ: Dim="<<Dim<<", Df="<<Densfact
			    <<",  STR="<<Stress<<" ";
			if (Dim==3 && Handflip==-1)
			    cout<<", flip";
			cout<<endl;
			
		    }	// if projection
	
		    // ---- Euclidean space adjustments ----
		    
		    #ifdef USE_OPENGL_GRAPHICS
			if (Graph)
			{
			    Draw.display_eucl(Dista);
			    Draw.display_coords(Model);
			}
		    #endif
		    
		    // detangling and RBA
		    if (Pieces.clu_no()>1)
		    {
			Tangiter=Params.i_value("Tangiter");	// reset
			Tangviol=Tangles.tangle_elim(Pieces, Model, TADJ, Tangiter);
			cout<<"TNGL: "<<Tangviol<<" (cyc="<<Tangiter<<")"<<endl;
			
			if (Tangiter)   // had to do detangling
			{
			    Model.dist_mat2(Dista);
			    Fakebeta.update(Dista, Polymer);
			    Steric.ideal_dist(Dista, Fakebeta, Restraints, Polymer, 
				    Pieces, Steric_::BETWEEN | Steric_::RESTR);
			    Steric.adjust_xyz(Dista, Model, Pieces, Steric_::BETWEEN);
			}
		    }
		    // accessibility
		    Access.solvent_xyz(Polymer, Pieces.hbond_bits(), Model);
		    
		    // 3D isotropic ellipsoidal density adjustment
		    if (Dim==3)
			Densfact=ellips_dens(Params.f_value("Density"), Pieces, Model);
		    
		    /* If the model is composed of several clusters, then first
		     * the clusters are adjusted ("WITHIN"), then the relative
		     * positions of the clusters ("BETWEEN"), finally all
		     * atoms move together ("ALL"). Sub-adjustments are for
		     * external restraints ("NMR" and 2o str), 
		     * then all.
		     */
		    cout<<"EUCL: ";
		    if (Pieces.clu_no()>1)
		    {

			// WITHIN-external
			Model.dist_mat2(Dista);
			Fakebeta.update(Dista, Polymer);
			Steric.ideal_dist(Dista, Fakebeta, Restraints, Polymer, 
				Pieces, Steric_::WITHIN | Steric_::REXT);
			Steric.adjust_xyz(Dista, Model, Pieces, Steric_::WITHIN);
			
			Model.dist_mat2(Dista);
			Fakebeta.update(Dista, Polymer);
			Steric.ideal_dist(Dista, Fakebeta, Restraints, Polymer, 
				Pieces, Steric_::WITHIN | Steric_::REXT | Steric_::SPECGRAD);
			Stress=Steric.adjust_xyz(Model, Speciter, Speceps, Noconv);
			if (Noconv || Stress<0.0)
			    Steric.adjust_xyz(Dista, Model, Pieces, Steric_::WITHIN);
			
			// WITHIN-all
			Model.dist_mat2(Dista);
			Fakebeta.update(Dista, Polymer);
			Steric.ideal_dist(Dista, Fakebeta, Restraints, Polymer, 
				Pieces, Steric_::WITHIN | Steric_::RESTR);
			Steric.adjust_xyz(Dista, Model, Pieces, Steric_::WITHIN);
			
			Model.dist_mat2(Dista);
			Fakebeta.update(Dista, Polymer);
			Steric.ideal_dist(Dista, Fakebeta, Restraints, Polymer, 
				Pieces, Steric_::WITHIN | Steric_::RESTR | Steric_::SPECGRAD);
			Stress=Steric.adjust_xyz(Model, Speciter, Speceps, Noconv);
			if (Noconv || Stress<0.0)	// on error or no convergence
			{
			    Steric.adjust_xyz(Dista, Model, Pieces, Steric_::WITHIN);
			    cout<<"IN=???";
			}
			else cout<<"IN="<<Stress;
			
			#ifdef USE_OPENGL_GRAPHICS
			    if (Graph) Draw.display_coords(Model);
			#endif
			
			// cluster hydrophobic moment rotation
    	    	    	hmom_clurot(Pieces, Polymer, Model);
		    
			// BETWEEN-external
			Model.dist_mat2(Dista);
			Fakebeta.update(Dista, Polymer);
			Steric.ideal_dist(Dista, Fakebeta, Restraints, Polymer, 
				Pieces, Steric_::BETWEEN | Steric_::REXT);
			Steric.adjust_xyz(Dista, Model, Pieces, Steric_::BETWEEN);
			
			// BETWEEN-all (RBA)
			Model.dist_mat2(Dista);
			Fakebeta.update(Dista, Polymer);
			Steric.ideal_dist(Dista, Fakebeta, Restraints, Polymer, 
				Pieces, Steric_::BETWEEN | Steric_::RESTR);
			Steric.adjust_xyz(Dista, Model, Pieces, Steric_::BETWEEN);
			
			#ifdef USE_OPENGL_GRAPHICS
			    if (Graph) Draw.display_coords(Model);
			#endif
			
		    }
		    
		    /* Adjusting all atoms together. If the model
		     * is just one piece, then this is carried out 3 times
		     * to compensate for the lost WITHIN/BETWEEN adjustments.
		     */
		    for (int i=0; i<(Pieces.clu_no()>1? 1: 3); i++)
		    {
			// ALL-external
			Model.dist_mat2(Dista);
			Fakebeta.update(Dista, Polymer);
			Steric.ideal_dist(Dista, Fakebeta, Restraints, Polymer, 
			    Pieces, Steric_::ALL | Steric_::REXT);
			Steric.adjust_xyz(Dista, Model, Pieces, Steric_::ALL);
		    
			Model.dist_mat2(Dista);
			Fakebeta.update(Dista, Polymer);
			Steric.ideal_dist(Dista, Fakebeta, Restraints, Polymer, 
			    Pieces, Steric_::ALL | Steric_::REXT | Steric_::SPECGRAD);
			Stress=Steric.adjust_xyz(Model, Speciter, Speceps, Noconv);
			if (Noconv || Stress<0.0)
			    Steric.adjust_xyz(Dista, Model, Pieces, Steric_::ALL);
			    
			// secondary structure adjustment (in 3D only!)
			if (Dim==3)
			{
			    Rmss=apply_secstruct(Pieces, Model);    // fit ideal secstr in 3D
			    cout<<" 2oSTR="<<Rmss;
			}
			// ALL-all
			Model.dist_mat2(Dista);
			Fakebeta.update(Dista, Polymer);
			Steric.ideal_dist(Dista, Fakebeta, Restraints, Polymer, 
			    Pieces, Steric_::ALL | Steric_::RESTR);
			Steric.adjust_xyz(Dista, Model, Pieces, Steric_::ALL);
		    
			Model.dist_mat2(Dista);
			Fakebeta.update(Dista, Polymer);
			Steric.ideal_dist(Dista, Fakebeta, Restraints, Polymer, 
			    Pieces, Steric_::ALL | Steric_::RESTR | Steric_::SPECGRAD);
			Stress=Steric.adjust_xyz(Model, Speciter, Speceps, Noconv);
			if (Noconv || Stress<0.0)	// on error or no convergence
			{
			    Steric.adjust_xyz(Dista, Model, Pieces, Steric_::ALL);
			    cout<<" ALL=???";
			}
			else cout<<" ALL="<<Stress;
		    }
		    cout<<endl;
		    
		    #ifdef USE_OPENGL_GRAPHICS
			if (Graph)
			{
			    Draw.display_eucl(Dista);
			    Draw.display_coords(Model);
			}
		    #endif
		    
		    // adjust CA:CA bonds and CA(i):CA(i+2) only
		    Model.dist_mat2(Dista);	// Fakebeta ignored, no update: only CAs
		    Steric.ideal_dist(Dista, Fakebeta, Restraints, Polymer, 
			Pieces, Steric_::ALL|Steric_::BOND);
		    Steric.adjust_xyz(Dista, Model, Pieces, Steric_::ALL);
		    
		    // same with Specgrad as well: usually converges after Willie's adjustment
		    Model.dist_mat2(Dista);
		    Steric.ideal_dist(Dista, Fakebeta, Restraints, Polymer, 
			Pieces, Steric_::ALL|Steric_::BOND|Steric_::SPECGRAD);
		    Steric.adjust_xyz(Model, Speciter, Speceps, Noconv);
    
		    // generate violation score (different from Stress)
		    Model.dist_mat2(Dista);
		    Fakebeta.update(Dista, Polymer);
		    Steric.ideal_dist(Dista, Fakebeta, Restraints, Polymer, 
			    Pieces, Steric_::ALL | Steric_::RESTR | Steric_::SCORE, &Euclsco);
		    Euclsco[Scores_::ACCESS].score(Access.score_xyz(Polymer, Model));
    
		    // save best if in 3D and score was acceptable and wasn't tangled
		    if (Dim==3)
		    {
			bool Tangled=Tangles.tangle_detect(Pieces, Model);
			if (!Tangled && (!Bestfound || Bestsco.accept_new(Euclsco)))
			{
			    Best=Model; Bestsco.update(Euclsco);
			    Distbest=Dista;
			    Bestfound++; 
			    Repriter=Reprojno=0;
			    cout<<"** BEST: "<<Bestsco<<endl;
			}
			else    // count iterations in 3D
			{
			    It3dno++; 
			    if (!Tangled && Bestfound)
				Repriter++;   // don't shortcut reproj if tangled
			    if (It3dno>=Params.i_value("Maxiter")) Exreason=EXIT_MAXITER;
			    if (Reprojno==2) Exreason=EXIT_REPROJ;
				cout<<"EUCL: "<<Euclsco<<endl;
			}
		    }
		    else cout<<"EUCL: "<<Euclsco<<endl;
	    
		    // leave early if in 3D and score was good
		    if (Dim==3 && Bestfound && Bestsco.is_exit())
			Exreason=EXIT_SCOREOK;
		    
		    // count overall iterations
		    Itno++;
		}	// end of try-block
		catch(Sigexcept_ Sigexc)	// signals land us here
		{
		    Signal=Sigexc.sigval();
		    Exreason=(Signal==SIGINT)? EXIT_CTRLC: EXIT_SIGNAL;
		}
	    }
	    while (Exreason==NOEXIT);   // end of big do-cycle
	    Sigproc.set_signal(SIG_DFL);    // don't catch signals any more
	
	    cout<<"EXIT: ";
	    switch (Exreason)
	    {
		case EXIT_SIGNAL: cerr<<"on signal "<<Signal<<endl; break;
		case EXIT_CTRLC: cerr<<"user interrupt requested\n"; break;
		case EXIT_SCOREOK: cerr<<"score convergence criterion satisfied\n"; break;
		case EXIT_MAXITER: cerr<<"maximal number of iterations reached\n"; break;
		case EXIT_REPROJ: cerr<<"no further improvement on 3D reprojection\n"; break;
		default: cerr<<"reason unknown (not implemented)\n"; break;
	    }
	    
	    // output
	    stop_timer();
	    cout<<"TIME: "<<time_string(timer_results(TS_UTIME|TS_STIME))<<endl;
	    if (Bestfound)
	    {
		cout<<"END: "<<Bestsco<<", Itno:"<<Itno<<"="<<(Itno-It3dno)<<"+"<<It3dno<<endl;
		
		#ifdef USE_OPENGL_GRAPHICS
		    if (Graph)
		    {
			Draw.display_eucl(Distbest);
			Draw.display_coords(Best);
		    }
		#endif

		// get output file
		Outname=Params.s_value("Outfnm");
		make_outname(Outname, Rcyc, "pdb");
		cout<<"SAVE: "<<Outname<<endl;
		pdb_result(Outname, Best, Polymer, Pieces, Bestsco);
		
		// write violation file
		Viollist_ Viollist;
		
		Outname=Params.s_value("Outfnm");
		make_outname(Outname, Rcyc, "viol");
		Fakebeta.update(Distbest, Polymer);
		Steric.ideal_dist(Distbest, Fakebeta, Restraints, Polymer, Pieces, 
		    Steric_::ALL | Steric_::RESTR | Steric_::SCORE, &Euclsco, &Viollist);
		Viollist.write_file(Outname);
		cout<<"VIOLS: "<<Outname<<endl;
		cout<<"\nRun "<<Rcyc<<" finished: "<<time_stamp()<<endl;
	    }
	    else if (Dim==3)   // no result and no signals: repeat run
	    {
		// save last conformation anyway, but no violation file is made
		Outname=Params.s_value("Outfnm");
		Outname+="_TEMPORARY";
		make_outname(Outname, Rcyc, "pdb");
		pdb_result(Outname, Model, Polymer, Pieces, Euclsco);
		cout<<"END: Temporary result, possibly tangled! Repeating run "<<Rcyc
		    <<endl<<"SAVE: "<<Outname<<endl;
		
		// attempt detangling
		if (Pieces.clu_no()>1)
		{
		    Tangiter=2*Params.i_value("Tangiter");	// reset to generous value
		    Tangviol=Tangles.tangle_elim(Pieces, Model, TADJ, Tangiter);
		    cout<<"TNGL: "<<Tangviol<<" (cyc="<<Tangiter<<")"<<endl;
		    Outname=Params.s_value("Outfnm");
		    Outname+="_DETANGLED";
		    make_outname(Outname, Rcyc, "pdb");
		    pdb_result(Outname, Model, Polymer, Pieces, Euclsco);
		    cout<<"SAVE: "<<Outname<<endl;
		}
		
		if (!Params.i_value("Randseed") && !Signal)
		{
		    // start from different random matrix
		    Rcyc--;
		}
		else cout<<endl;
	    }
	    cerr<<flush; cout<<flush;
	    
	}	    // for Rcyc (all simulations)
	
	// close windows if graphics was on
	#ifdef USE_OPENGL_GRAPHICS
	    if (Graph) Draw.close_window();
	#endif

    	if (Sigproc.is_child())
	{
	    close(Logfd);
	    exit(Signal);
	}
    }	    // else: child or single simulation process

    // reset the floating-point output for cout and cerr
    cout.flags(Coutf); cerr.flags(Cerrf);
    cout.precision(Oldoutprec); cerr.precision(Olderrprec);
    
    return(Signal);
}
// END of dragon_run()

// ---- Auxiliaries ----

/* merge_distmat: mixes the distances in Bestdist into Dist so
 * that local distances (diagonals close to the main diag) will be
 * more or less the same, whereas global distances will come from
 * Bestdist. Based on Willie's idea: the weighting function
 * is an inverted Gaussian.
 */
static void merge_distmat(const Trimat_& Bestdist, Trimat_& Dist)
{
    static const unsigned int GSTEP=10;
    static double Mix[GSTEP], Mixm1[GSTEP];
    static char Init=1;
    
    register double Nij, Dij, Ndij;
    register unsigned int d, i, j, Ptno=Bestdist.rno();
    
    // fill up the blend coefficient arrays on first call
    if (Init)
    {
	for (d=0; d<GSTEP; d++)
	{
	    Mix[d]=exp(-(d*d/double(GSTEP)));
	    Mixm1[d]=1.0-Mix[d];
	}
	Init=0;
    }
    
    // adjust from second off-diagonal
    for (d=2; d<Ptno; d++)
    {
	if (d>GSTEP)	// "global" distances
	{
	    for (i=d; i<Ptno; i++)
	    {
		j=i-d;
		Dist[i][j]=Bestdist[i][j];
	    }
	    continue;
	}
	
	// "local" distances
	for (i=d; i<Ptno; i++)
	{
	    j=i-d;
	    Dij=sqrt(Dist[i][j]); Nij=sqrt(Bestdist[i][j]);

	    Ndij=Mix[d-1]*Dij+Mixm1[d-1]*Nij;
	    Dist[i][j]=Ndij*Ndij;
	}
    }
}
// END of merge_distmat()

// ==== END OF PROGRAM Dragon.c++ ====
