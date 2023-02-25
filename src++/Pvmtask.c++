#ifdef USE_PVM
// ==== PROJECT DRAGON: METHODS Pvmtask.c++ ====

/* Task management on the Parallel Virtual Machine (PVM).
 * PVM is free software originally developed at the
 * Oak Ridge National Laboratory, Oak Ridge, Tennessee, USA.
 */

// 5-Aug-2000. Andras Aszodi
// Multiple slaves on MP nodes: support by Jan Hungershoefer

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <strstream.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <signal.h>
#include <fcntl.h>

/* NOTE: when compiled with the "-ansi" flag on SGIs, gethostname()
 * will not be recognised unless _SGI_SOURCE is defined.
 */
#if defined (__sgi) && !defined (_SGI_SOURCE)
    #define _SGI_SOURCE
#endif
#include <unistd.h>
#include <stdio.h> // for popen, pclose
#if defined(__sgi) && defined(_SGI_SOURCE)
    #undef _SGI_SOURCE
#endif

#ifdef __sun
    #include <kstat.h>
    #include <sys/cpuvar.h>
    #include <sys/param.h>
    #define loaddouble(la) ((double)(la) / FSCALE)
#endif

// ---- MODULE HEADERS ----

#include "Pvmtask.h"
#include "String.h"

// ---- DEFINITIONS ----

// Define this if you want average load detection on multiprocessor nodes
#define DETERMINE_AVGLOAD

// Size of hostname string
#ifndef MAXHOSTNAMELEN
#define MAXHOSTNAMELEN 65
#endif

// ==== Pvmtask_ METHODS ====

// ---- Static initialisation ----

unsigned int Pvmtask_::Objno=0;	// only 1 object per pgm is allowed

// ---- Constructor and destructor ----

// Inits the object to empty.
Pvmtask_::Pvmtask_(): 
    Pvmstat(NO_PVM), Slaves(NULL), Slaveno(0), 
    Idstr(NULL), Slavexec(NULL), Mastertid(-1),
    Hosts(NULL), Nhosts(0)
{
    if (Objno)
    {
	cerr<<"\n! Pvmtask_(): Only one object per program is allowed!\n";
	Pvmstat=NO_PVM;
    }
    else Objno=1;
}

// enrol_pvm(Slexenm):
// enrols the process in PVM and finds out if it runs as a master or slave.
// The optional Slexenm argument specifies the slave executable name.
void Pvmtask_::enrol_pvm(const char *Slexenm)
{
    // turn off automatic error reporting
    pvm_setopt(PvmAutoErr, 0);
    
    if (!no_pvm()) return;
    
    // enrol into PVM
    Tid=pvm_mytid();
    if (PvmSysErr==Tid)
    {
	cerr<<"\n? Pvmtask_(): PVM not available\n";
	Pvmstat=NO_PVM; // PVM is unavailable
	return;
    }
    
    // generate ID string "[Tid]@hostname"
    Idstr=new char [MAXHOSTNAMELEN+20];	// make it large enough
    memset(Idstr, 0, MAXHOSTNAMELEN+20);    // just in case...
    ostrstream(Idstr, 20)<<hex<<Tid<<dec<<"@";
    if (gethostname(Idstr+strlen(Idstr), MAXHOSTNAMELEN)<0)
	strcpy(Idstr+strlen(Idstr), "unknown");
    
    // decide if master or slave
    Mastertid=pvm_parent();
    if (Mastertid==PvmNoParent)
    {
	// I am the master
	Mastertid=-1;
	Pvmstat=MASTER;
	
	if (Slexenm!=NULL && strlen(Slexenm))
	{
	    Slavexec=new char [strlen(Slexenm)+1];
	    strcpy(Slavexec, Slexenm);
	}
	else
	{
	    cerr<<"\n? Pvmtask_::enrol_pvm(MASTER): Unspecified slave executable name\n";
	}
	cout<<"Master "<<id_str()<<" ready.\n";
    }
    else
    {
	// I am a slave
	Pvmstat=SLAVE;
	
	// open personal logfile: the name is the Idstr
	int Logfd=open(id_str(), O_CREAT|O_WRONLY, 0644);
	if (Logfd<0)
	    prt_error("enrol_pvm()", "Cannot open slave log");
	else
	{
	    cout<<flush; cerr<<flush;
	    dup2(Logfd, STDOUT_FILENO);
	    dup2(Logfd, STDERR_FILENO);
	}

	cout<<"Slave "<<id_str()
	    <<" , Master=["<<hex<<Mastertid<<dec<<"] ready"<<endl;
    }
}
// END of enrol_pvm()

/* When the master exits, it attempts to kill all remaining slaves
 * before dying. This assumes that the calling object was global
 * and dies only with the program.
 */
Pvmtask_::~Pvmtask_()
{
    if (is_master())
    {
	check_slaves();
	for (unsigned int i=0; i<Slaveno; i++)
	{
	    pvm_kill(Slaves[i].Tid);
	    cout<<"Killed slave ["<<hex<<Slaves[i].Tid
		<<dec<<"]: Jobno="
		<<Slaves[i].Jobs<<endl;
	}
	delete [] Slaves; Slaves=NULL;
	delete [] Slavexec; Slavexec=NULL;
    delete [] Hosts; Hosts=NULL;
    }
    if (!no_pvm())
    {
	cout<<"Leaving PVM: "<<id_str()<<endl;
	pvm_exit();
    }
    delete [] Idstr;
    Objno--;
}
    
// ---- Task management ----

/* spawn_slaves(): attempts to launch one slave on each node
 * of the virtual machine. Checks old slaves and virtual machine
 * status, and spawns slaves on the nodes which have none.
 * The slave executable's name is Slavenm, no command-line
 * arguments are passed: if left unspecified, then the
 * internal Slavexec string is used. Can be invoked by the master only.
 * Return value: if PVM is not available, then <0 is returned, 
 * otherwise the no. of newly created slaves is returned.
 */
int Pvmtask_::spawn_slaves(const char *Slavenm)
{
    if (!is_master()) return(0);
    
    check_slaves();	// update slave info
    
    // update machine configuration info
    int Nhost=0, Narch=0, Result=0;
    struct pvmhostinfo *Pvmhosts=NULL;
    
    Result=pvm_config(&Nhost, &Narch, &Pvmhosts);
    if (Result<0)
    {
	prt_error("spawn_slaves(): [pvm_config]", Result);
	return(Result);
    }
    
    // prepare argument array: pass -M flag to slaves, too! (5-Aug-2000)
    char *Slaveargs[2]={"-M", NULL};
    
    /* Check each machine in turn: if a slave already runs on one,
     * then its tid and host daemon tid are just copied into the
     * new slave array. Otherwise, a new slave is spawned
     */
    Slave_ *Newslaves=NULL;
    int h, s, Newtid, Newslaveno, Fresh;

    int Soh;  // slaves on host
    int Newslavesize = 2*Nhost; // assume 2 processors per host as default
    Newslaves=new Slave_[ Newslavesize ];

    for (Fresh=Newslaveno=h=0; h<Nhost; h++)
    {
        // count the number of slaves on host
        for (s=Soh=0; s<Slaveno; s++)
        {
            // copy descriptor if running on the current host
            if(Slaves[s].Hosttid==Pvmhosts[h].hi_tid)
            {
                // enlarge Slave array if neccessary
                if(Newslaveno>=Newslavesize)
                {
                    Newslavesize *= 2;
                    Slave_ *tmp = Newslaves;
                    Newslaves = new Slave_[Newslavesize];
                    for(int t=0; t<Newslavesize; t++)
                        Newslaves[t] = tmp[t];
                    delete [] tmp;
                } 
                Newslaves[Newslaveno]=Slaves[s];	// copy slave descriptor
                Newslaveno++;   // count slaves
                Soh++; // count slaves running on that host
            }
        }
        
        // check if number of processors of host is already known      
        int hh; // index in Hosts array
        for(hh=0; hh<Nhosts && Hosts[hh].Tid!=0; hh++)
        {
            if(Hosts[hh].Tid==Pvmhosts[h].hi_tid)
                break; // host known
        }
        if(hh==Nhosts)
        {
            // reached end of Hosts array and host not found
            // enlarge array
            Host_ *tmp = Hosts;
            Nhosts++;
            Hosts = new Host_[ Nhosts ];
            // copy old info
            for(int t=0; t<Nhosts-1; t++)
                Hosts[t] = tmp[t];
            delete [] tmp;
            // store host
            Hosts[hh].Tid = Pvmhosts[h].hi_tid;
            Hosts[hh].Ncpu = 1;
        }

        // spawn slaves on host if not enough (==number of processors) Dragon jobs there
        while(Soh < Hosts[hh].Ncpu)
        {
            // spawn a new slave on h-th host
            Result=pvm_spawn((Slavenm==NULL)? Slavexec: (char *)Slavenm, 
		    Slaveargs,  PvmTaskHost, 
		    Pvmhosts[h].hi_name, 1, &Newtid);
            if (Result<0)
            {
                prt_error("spawn_slaves(): [pvm_spawn]", Result);   // sys error
                continue;
            }
            if (!Result)    // sys OK but task didn't start, explanation in Newtid
            {
                prt_error("spawn_slaves(): [pvm_spawn]", Newtid);
                continue;
            }

            // launched successfully
            // enlarge Newslaves array if neccessary
            if(Newslaveno>=Newslavesize)
            {
                Newslavesize *= 2;
                Slave_ *tmp = Newslaves;
                Newslaves = new Slave_[Newslavesize];
                for(int t=0; t<Newslavesize; t++)
                  Newslaves[t] = tmp[t];
                delete [] tmp;
            } 
            Newslaves[Newslaveno]=Slave_(Newtid, Pvmhosts[h].hi_tid);    // Jobstatus=-1
            Newslaveno++;   
            Fresh++;

            // wait for the new slave to determine number of cpus on host
            // and send this information back
            // store information for next iteration
            Hosts[hh].Ncpu = recv_ncpus(Newtid);
            cout << "Slave ("<<Newslaveno<<"/"<<Hosts[hh].Ncpu<<") [" 
		<< hex << Newslaves[Newslaveno-1].Tid << dec 
		<< "] started on host "
                << Pvmhosts[h].hi_name << endl;
            Soh++;
        }
    }
    
    // replace old list with new
    delete [] Slaves; Slaves=Newslaves;
    Slaveno=Newslaveno;
    return(Fresh);  // no. of freshly created slaves
}
// END of spawn_slaves()

// ---- Signal handling ----

/* signal_pvm(): catches signals in the master task in PVM
 * and sends them on to the slaves. The slaves use the
 * "standard" signal_handler() function in "Sigproc".
 * The logic is that most signals terminate the slaves, 
 * whereupon the master will automatically exit; when
 * SIGINT <Ctrl-C> is caught, the slaves will be terminated.
 * NOTE: there must be a global Pvmtask_ called Pvmtask for this to work.
 */
void signal_pvm(int Sigtype) throw(Sigexcept_)
{
    extern Pvmtask_ Pvmtask;
    
    if (!Pvmtask.is_master() || !Pvmtask.check_slaves()) return;
    
    // print an explanatory message
    signal_message(Sigtype);
    
    // just pass the caught signal to the slaves
    for (unsigned int Sl=0; Sl<Pvmtask.Slaveno; Sl++)
	pvm_sendsig(Pvmtask.Slaves[Sl].Tid, Sigtype);
    
    // reinstall itself as handler for the signal just caught
    signal(Sigtype, (SIG_PF)signal_pvm);
    throw(Sigexcept_(Sigtype));	// get out
}
// END of signal_pvm()

// ---- Communication ----

/* NOTE: Parameter consistency is achieved as follows. The 
 * master process is responsible for communicating with the
 * user, either via the command prompt or through a command
 * file (see "Clip" for details). Just before running a
 * simulation, a parameter file is constructed from the
 * values of the changed parameters and this is multicast
 * to the slaves as a string. The slaves' Params_ object can
 * read this via istrstream. Additionally, the contents of the
 * changed data files are converted to strings and multicast
 * as separate messages (tagged with the corresponding public Filetags_)
 * to the slaves. Finally, RUN messages are sent to each of
 * the available slaves by send_jobs(). The slaves send back
 * DONE messages to the master upon completion. When slaves
 * wait for instructions, they send SLAVE_READY messages to
 * the master.
 */

/* send_params(P): first sends a string containing a modified
 * parameter file (the output of P.list_changed())
 * to all active slaves, telling them about all parameters
 * which were changed since the last run. This message is
 * tagged as PARAMS. If no modifications were done, then
 * just an int=0 is sent. If any of the
 * data file names were modified, then the contents of the
 * new files will be converted to strings and sent to the
 * slaves with the appropriate Filetags_.
 * Each string is preceded by its length to allow dynamic
 * allocation on the receiving end.
 * No action is taken if invoked in a non-master process.
 * NOTE: all parameters will be marked as not-changed after the call.
 * Return value: <0 on error, total no. of bytes transmitted otherwise.
 * 
 * send_params(P, Newslaves, Newno): sends the complete parameter set
 * (not only the changed ones) to all TIDs in Newslaves
 * (all in all Newno).
 * This version of the method is used when slave(s) had to be
 * re-started (due to config changes during a run).
 */
int Pvmtask_::send_params(Params_& P)
{
    static bool Chainchg=true;

    if (!is_master()) return(0);
    
    int s, Result=0, Msglen, Bytes=0;
    
    // send list of changed parameters
    ostrstream Chg;
    unsigned int Chgno=P.list_changed(Chg);
    char *Msgstr=Chg.str(); // "frozen"
    
    pvm_initsend(PvmDataDefault);
    Msglen=Chgno? strlen(Msgstr): 0;
    pvm_pkint(&Msglen, 1, 1);
    if (Chgno) pvm_pkstr(Msgstr);  // send only if non-empty
    delete [] Msgstr;
    
    check_slaves(); // update living slaves' list
    if (!Slaveno)
    {
	prt_error("send_params(P)", "All slaves are dead");
	return(0);
    }
    int *Slavetids=new int [Slaveno];	// simple array for PVM's benefit
    for (s=0; s<Slaveno; s++)
    {
	Slavetids[s]=Slaves[s].Tid;
	if (Slaves[s].Jobs<0) Slaves[s].Jobs=0;	// lose "virginity"
    }
    
    /* NOTE: there are lots of pvm_mcast() operations below,
     * all of which can go wrong. We just jump to the
     * end of the routine using  "goto PVMSEND_ERROR"
     * when an error is detected.
     */
    
    Result=pvm_mcast(Slavetids, Slaveno, PARAMS);   // send to all active slaves
    if (Result<0) goto PVMSEND_ERROR; else Bytes+=Msglen;
    
    /* Send contents of changed data files. Follow the
     * dependency logic in "init_dragon()" in Dragon.c++.
     * Note that the P.s_value() methods invoked in 
     * send_files() will reset the Changed status var
     */
    Chainchg=true; // def moved upwards: goto crossed it 12-Jul-1997
    if (P.changed("Alnfnm") || P.changed("Masterno"))
    {
	Result=send_files(Slavetids, Slaveno, P, "Alnfnm", ALN);
	if (Result<0) goto PVMSEND_ERROR; else Bytes+=Result;
	P.i_value("Masterno");	// touch Masterno as well
	Chainchg=true;
    }
    
    if (P.changed("Phobfnm"))
    {
	Result=send_files(Slavetids, Slaveno, P, "Phobfnm", PHO);
	if (Result<0) goto PVMSEND_ERROR; else Bytes+=Result;
    }
    if (P.changed("Volfnm"))
    {
	Result=send_files(Slavetids, Slaveno, P, "Volfnm", VOL);
	if (Result<0) goto PVMSEND_ERROR; else Bytes+=Result;
    }
    if (P.changed("Adistfnm"))
    {
	Result=send_files(Slavetids, Slaveno, P, "Adistfnm", ACD);
	if (Result<0) goto PVMSEND_ERROR; else Bytes+=Result;
    }
    if (P.changed("Simfnm"))
    {
	Result=send_files(Slavetids, Slaveno, P, "Simfnm", SIM);
	if (Result<0) goto PVMSEND_ERROR; else Bytes+=Result;
    }
    if (Chainchg || P.changed("Accfnm"))
    {
	Result=send_files(Slavetids, Slaveno, P, "Accfnm", ACC);
	if (Result<0) goto PVMSEND_ERROR; else Bytes+=Result;
    }
    
    if (Chainchg || P.changed("Restrfnm") || P.changed("Homfnm")
	    || P.changed("Maxdist") || P.changed("Minsepar"))
    {
	Result=send_files(Slavetids, Slaveno, P, "Restrfnm", RESTR);
	if (Result<0) goto PVMSEND_ERROR; else Bytes+=Result;
	Result=send_files(Slavetids, Slaveno, P, "Homfnm", HOM);
	if (Result<0) goto PVMSEND_ERROR; else Bytes+=Result;
    }
    
    if (Chainchg)
    {
	Result=send_files(Slavetids, Slaveno, P, "Sstrfnm", SSTR);
	if (Result<0) goto PVMSEND_ERROR; else Bytes+=Result;
	Chainchg=false;
    }
    else if (P.changed("Sstrfnm"))
    {
	Result=send_files(Slavetids, Slaveno, P, "Sstrfnm", SSTR);
	if (Result<0) goto PVMSEND_ERROR; else Bytes+=Result;
    }
    
    // mark the rest of the parameters unchanged
    P.reset_changed();
    
    PVMSEND_ERROR:
    delete [] Slavetids;
    if (Result<0)   // one of the send-s were bad
    {
	prt_error("send_params(P)", Result);
	return(Result);
    }
    else return(Bytes);  // apparently everything went fine
}

int Pvmtask_::send_params(Params_& P, const int *Newslaves, int Newno)
{
    if (!is_master() || Newslaves==NULL || !Newno) return(0);
    
    int Result, Msglen, Bytes=0;
    
    // send parameter list
    ostrstream Parstr;
    Parstr<<P; // write the complete parameter list 
    char *Msgstr=Parstr.str(); // "frozen"
    
    pvm_initsend(PvmDataDefault);
    Msglen=strlen(Msgstr);
    pvm_pkint(&Msglen, 1, 1);
    pvm_pkstr(Msgstr);
    delete [] Msgstr;
    
    Result=pvm_mcast((int *)Newslaves, Newno, PARAMS);   // send to slave
    if (Result<0)
    {
	prt_error("send_params(P, Newslaves, Newno)", Result);
	return(Result);
    }
    else Bytes+=Msglen;
    
    // Send contents of data files.
    Result=send_files(Newslaves, Newno, P, "Alnfnm", ALN);
    if (Result<0) return(Result); else Bytes+=Result;
    Result=send_files(Newslaves, Newno, P, "Phobfnm", PHO);
    if (Result<0) return(Result); else Bytes+=Result;
    Result=send_files(Newslaves, Newno, P, "Volfnm", VOL);
    if (Result<0) return(Result); else Bytes+=Result;
    Result=send_files(Newslaves, Newno, P, "Adistfnm", ACD);
    if (Result<0) return(Result); else Bytes+=Result;
    Result=send_files(Newslaves, Newno, P, "Simfnm", SIM);
    if (Result<0) return(Result); else Bytes+=Result;
    Result=send_files(Newslaves, Newno, P, "Accfnm", ACC);
    if (Result<0) return(Result); else Bytes+=Result;
    Result=send_files(Newslaves, Newno, P, "Restrfnm", RESTR);
    if (Result<0) return(Result); else Bytes+=Result;
    Result=send_files(Newslaves, Newno, P, "Homfnm", HOM);
    if (Result<0) return(Result); else Bytes+=Result;
    Result=send_files(Newslaves, Newno, P, "Sstrfnm", SSTR);
    if (Result<0) return(Result); else Bytes+=Result;
    
    return(Bytes);  // apparently everything went fine
}

// END of send_params()

/* send_files(): converts the contents of a file whose name is
 * the parameter called Pname in the parameter object P to a string
 * and sends the result to all slave processes 
 * whose TID-s are listed in the array Slavetids (Slno, all in all).
 * The message consists of an integer holding the length of the string and the
 * string itself (which is not sent if the length was 0). The
 * message is tagged with Tag.
 * Returns <=-2 on file errors, -1 if all slaves are dead, 
 * or the length of the file sent (>=0) if OK.
 */
int Pvmtask_::send_files(const int* Slavetids, int Slno, Params_& P,
	const char *Pname, Filetags_ Tag)
{
    if (!is_master() || Slavetids==NULL || !Slno) return(-1);
    
    // get the filename (may be "")
    const char *Fname=P.s_value(Pname);	// Changed reset (1st access)
    int Fsize=0;	// length of file
    char *Fstr=NULL;	// contents of the file
    
    if (Fname==NULL)
	prt_error("send_files(): Cannot find parameter ", Pname);
    else if (strlen(Fname))
    {
	// Fname!="", get the size of the given file 
	struct stat Statbuf;
	if (stat(Fname, &Statbuf))
	    prt_error("send_files(): Cannot stat ", Fname);
	else Fsize=Statbuf.st_size;
    }
    
    // open the file if available and read contents into Fstr
    if (Fsize)
    {
	ifstream Inf(P.s_value(Pname));
	if (!Inf)
	    prt_error("send_files(): Cannot open ", Fname);
	else
	{
	    Fstr=new char [Fsize+1];  // big input buffer
	    Inf.read(Fstr, Fsize);   // unformatted input
	    Fsize=Inf.gcount();	    // actual no. of bytes read
	    Fstr[Fsize]='\0';	    // terminate 
	    Inf.close();
	}
    }
    
    // pack and send to slaves in Slavetids
    pvm_initsend(PvmDataDefault);
    pvm_pkint(&Fsize, 1, 1);
    if (Fsize) pvm_pkstr(Fstr);	// send only if contained something
    delete [] Fstr;
    
    int Result=pvm_mcast((int *)Slavetids, Slno, Tag);
    if (Result<0) prt_error("send_files()", Result);
    return(Fsize);    // OK (>=0)
}
// END of send_files()

/* recv_params(): in a slave process, this method should be
 * called when wait_master(PARAMS) returns successfully
 * from the listening. The PARAMS message holds a parameter
 * list string which will be read into P directly.
 * Note that PARAMS is in the active buffer and must be read
 * immediately after receiving with wait_master().
 * Return value: the length of the message.
 */
int Pvmtask_::recv_params(Params_& P) const
{
    if (!is_slave()) return(0);
    
    // message is here, unpack and read
    int Result=pvm_nrecv(Mastertid, PARAMS);
    if (Result<0)
    {
	prt_error("recv_params(): [pvm_nrecv]", Result);
	return(0);
    }
    if (!Result)
    {
	prt_error("recv_params()", "Message did not arrive");
	return(0);
    }
    
    int Len=0;
    char *Pstr=NULL;
    pvm_upkint(&Len, 1, 1);
    if (Len)	// Len==0 means that no further data are sent
    {
	Pstr=new char [Len+1];
	pvm_upkstr(Pstr);
	istrstream Ifs(Pstr, Len+1); Ifs>>P; // read changed params from Pstr
	delete [] Pstr;
    }
    return(Len);
}
// END of recv_params()

/* recv_filestr(): receives one of the data files as a string,
 * indicated by the Tag. Waits until the message arrives
 * or the master dies. Returns a char ptr pointing
 * to a string containing the corresponding ASCII data, or ""
 * on error. Empty files also return a "". NOTE: the 
 * returned pointer points to a static area which will
 * be overwritten upon each invocation of the method.
 * Does nothing if not invoked in a slave.
 */
char* Pvmtask_::recv_filestr(int Tag) const
{
    static char *Pstr=NULL;
    static int Plen=0;
    
    if (Pstr==NULL)	// static initialisation (done once)
    {
	Pstr=new char [1]; 
	Plen=0;
    }
    *Pstr='\0';  // make the "" string (always done)
    
    if (!is_slave()) return(Pstr);
    
    int Bufid=wait_master(Tag);	// return if message received or master died
    if (Bufid<=0) return(Pstr);
    
    // message is here, unpack and read
    int Result=pvm_nrecv(Mastertid, Tag);
    if (Result<0)
    {
	prt_error("recv_filestr(): [pvm_nrecv]", Result);
	return(Pstr);
    }
    if (!Result)
    {
	prt_error("recv_filestr()", "Message did not arrive");
	return(Pstr);
    }
    
    int Len=0;
    pvm_upkint(&Len, 1, 1);
    if (Len)	// message was non-empty (otherwise "" will be returned)
    {
	// do we need to grow Pstr?
	if (Len>Plen)
	{
	    delete [] Pstr; Pstr=new char [Len+1];
	    Plen=Len;
	}
	pvm_upkstr(Pstr);
    }
    return(Pstr);	    // normal return
}
// END of recv_filestr()

// ---- Job distribution ----

/* send_jobs(): when invoked on a master, distributes the Runno
 * jobs among the available slaves so that once a slave finishes, 
 * it will receive the next job on the queue ("Pool-of-Tasks"
 * paradigm). Parameters will be supplied in Params.
 * Slaves will be re-spawned in every cycle if necessary
 * but it is assumed that they have a consistent parameter set
 * upon entering this method. Signal catching is done here.
 * Return value: the signal caught during execution, 0 if OK.
 * Also returns the actual number of jobs done in Runno.
 */
int Pvmtask_::send_jobs(Params_& Params, unsigned int& Runno)
{
    if (!is_master() || !Runno) return(0);
    
    /* Every job must know whether it is done or not, and
     * if being processed, then on which slave.
     */
    struct Jobstatus_
    {
        // job status symbols
        enum Jstat_ {TOBE_DONE, BEING_DONE, DONE};

	int Tid;    // the slave process TID
	Jstat_ Stat;	// the job status
	
	Jobstatus_(): Tid(0), Stat(TOBE_DONE) {}    // ctor
	// set status to TOBE_DONE if processing slave is  N/A any more
	void check_dead()
	{
	    if (Stat==BEING_DONE && 
		(PvmSysErr==pvm_pstat(Tid) || PvmNoTask==pvm_pstat(Tid)))
		    Stat=TOBE_DONE;
	}
    };
    
    // signal processing facilities
    extern Sigproc_ Sigproc;
    
    int Newno, Jobno, Done=0, Sl, Stid, Prevjobno, Cycno, Signal=0;
    Jobstatus_* Jobs=new Jobstatus_ [Runno];
    int *Newsl=NULL;
    
    Sigproc.set_signal((SIG_PF)signal_pvm);    // set the PVM master signal trap
    try		// do the job, watch out for signal exceptions
    {
	while (Done<Runno)
	{
	    /* set status to TOBE_DONE for all jobs which
	     * are not completed and the processing slave died.
	     * These will be re-sent to other, healthy slaves.
	     */
	    for (Jobno=0; Jobno<Runno; Jobno++)
		Jobs[Jobno].check_dead();
	    check_slaves();
	    
	    Newno=spawn_slaves();	// re-spawn slaves if necessary
	    if (Newno<0) break;	// problem
    
	    // update parameters for all new slaves
	    if (Newno>0)
	    {
		Newsl=new int [Newno];
		int i;
		for (i=Sl=0; Sl<Slaveno; Sl++)
		    if (Slaves[Sl].Jobs==-1)	// needs full params
		    {
			Newsl[i++]=Slaves[Sl].Tid;
			Slaves[Sl].Jobs=0; // lost virginity...
			cout<<"New slave: ["<<hex<<Slaves[Sl].Tid<<dec<<"]\n";
		    }
		int Result=send_params(Params, Newsl, Newno);	// update new slaves
		delete [] Newsl;
		if (Result<0) break;    // something went wrong
	    }
	    
	    // pick an idle slave and give it a new job
	    for (Sl=0; Sl<Slaveno; Sl++)
	    {
		Stid=Slaves[Sl].Tid;
		// consume the SLAVE_READY message sent by this slave
		if (0<pvm_nrecv(Stid, SLAVE_READY))
		{
		    // check if this slave finished a job
		    if (pvm_nrecv(Stid, SLAVE_DONE))
		    {
			Prevjobno=0;
			pvm_upkint(&Prevjobno, 1, 1);
			Jobs[Prevjobno-1].Stat=Jobstatus_::DONE;
			Cycno=0;
			while(0<pvm_nrecv(Stid, SLAVE_RUNNING)) pvm_upkint(&Cycno, 1, 1);
			cout<<"Job "<<Prevjobno<<" completed: ["<<hex<<Jobs[Prevjobno-1].Tid<<dec
			    <<"]";
			if (Cycno) cout<<": "<<Cycno;
			cout<<endl;
			Slaves[Sl].Jobs++; Done++;
		    }
		    
		    // find first not-done job and send it
		    for (Jobno=0; Jobno<Runno && Jobs[Jobno].Stat!=Jobstatus_::TOBE_DONE; Jobno++);
		    if (Jobno<Runno)    // found a not-done job
		    {
			Jobs[Jobno].Stat=Jobstatus_::BEING_DONE;    // Jobno will be processed
			Jobs[Jobno].Tid=Stid;  // by this slave
			
			// send the job number
			Jobno++;	// jobs are [1..Runno]
			pvm_initsend(PvmDataDefault);
			pvm_pkint(&Jobno, 1, 1);	// send the job number only
			pvm_send(Stid, RUN);
			cout<<"Job "<<Jobno<<" sent: ["<<hex<<Stid<<dec<<"]"<<endl;
		    }
		    continue;
		}	    // if (SLAVE_READY)
		
		// if this slave is running, get info about cycle number
		Cycno=0;
		while(0<pvm_nrecv(Stid, SLAVE_RUNNING)) pvm_upkint(&Cycno, 1, 1);
		if (Cycno) cout<<"["<<hex<<Stid<<dec<<"]: "<<Cycno<<endl;
    
	    }	// for (Sl)
	    if (Done<Runno) sleep(1);   // wait 1 sec until next check
	}	    // while
    }	    // end of try block
    catch(Sigexcept_ Sigexc)	    // signals land us here
    {
	Signal=Sigexc.sigval();
    }
    Sigproc.set_signal(SIG_DFL);    // reset default traps
    
    delete [] Jobs; Runno=Done;
    return(Signal);
}
// END of send_jobs()

/* recv_job(): when invoked by a slave immediately after wait_master()
 * returned with the news that a RUN message has been sent, this
 * method reads the message and returns the job number (>0)
 * or a negative int on error.
 */
int Pvmtask_::recv_job() const
{
    if (!is_slave()) return(-1);
    
    int Result=pvm_nrecv(Mastertid, RUN);
    if (Result<0)
    {
	prt_error("recv_job() [pvm_nrecv]", Result);
	return(Result);
    }
    if (!Result)
    {
	prt_error("recv_job()", "No RUN message received");
	return(-1);
    }
    
    int Jobno=0;
    pvm_upkint(&Jobno, 1, 1);
    return(Jobno);
}
// END of recv_job()

/* job_status(): when invoked by a slave, this method sends a
 * message with Tag to the master. The message consists of an
 * integer Num. If Tag==SLAVE_RUNNING, then it is the number of cycles
 * already done. If Tag==SLAVE_DONE, then it is the job number.
 * Return value: <0 indicates an error.
 */
int Pvmtask_::job_status(int Tag, int Num) const
{
    if (!is_slave()) return(-1);
    
    pvm_initsend(PvmDataDefault);
    pvm_pkint(&Num, 1, 1);	// send the number
    return(pvm_send(Mastertid, Tag));
}
// END of job_status()

// ---- Auxiliaries ----

/* check_slaves(): when invoked by the master, this method
 * checks the process status of the slaves.
 * The dead tid-s are removed from the Slavetids
 * array and Slaveno (which is returned) is adjusted.
 * Slavetids[0..Slaveno-1] always contains the compact list
 * of all active slave tids after the call.
 * Private
 */
int Pvmtask_::check_slaves()
{
    if (!is_master()) return(Slaveno);
    if (Slaves==NULL) return(0); // virgin master
    
    register unsigned int i, k, Deadno=0;
    int Stid, Result;
    
    for (i=k=0; i<Slaveno; i++)
    {
	    Stid=Slaves[i].Tid;  
	    Result=pvm_pstat(Stid);
	    if (Result==PvmSysErr || Result==PvmNoTask)
	    {
	        cout<<"Slave ["<<hex<<Stid<<dec<<"] died: Jobno="
		        <<Slaves[i].Jobs<<endl;
	        Deadno++;
	    }
	    else if (Result<0)
	        prt_error("check_slaves()", Result);
	    else
	    {
	        Slaves[k]=Slaves[i];   // compress array
	        k++;
	    }
    }
    Slaveno-=Deadno;
    return(Slaveno);
}
// END of check_slaves()

/* slave_ready(): this method tells the master that the slave
 * has received all the parameters and is ready to run the
 * simulations.
 * Return value: <0 on error.
 */
int Pvmtask_::slave_ready() const
{
    if (!is_slave()) return(-1);
    
    // tell the master that this slave is ready
    pvm_initsend(PvmDataDefault);
    pvm_pkint((int*)&Tid, 1, 1);
    int Result=pvm_send(Mastertid, SLAVE_READY);
    if (Result<0)
    {
	prt_error("slave_ready()", Result);
	return(Result);
    }
    return(Tid);
}
// END of slave_ready()

/* recv_ncpus(): for the master only
 * waits for the freshly spawned slave to send the 
 * number of cpus found on its node
 */
int Pvmtask_::recv_ncpus(int Tid)
{
    int Ncpus=1;
    int Result=pvm_recv(Tid, CPUCNT);
    if (Result<0)
    {
        prt_error("recv_ncpus()", Result);
        return(Ncpus);
    }
    Result=pvm_upkint(&Ncpus, 1, 1);
    if (Result<0)
    {
        prt_error("recv_ncpus()", Result);
        return(Ncpus);
    }
    return (Ncpus);
}
// END of recv_ncpus()

/* send_ncpus(): for the slaves only (machine dependant)
 * determines the number of cpus on the machine
 */
int Pvmtask_::send_ncpus()
{
    if (!is_slave()) return(-1);
    int Ncpus=1, avgload=0;
    
// Here begins an architecture-dependent section:
// figuring out number of CPU-s and average system load
#if defined(__sun)
    Ncpus = sysconf(_SC_NPROCESSORS_ONLN);
#ifdef DETERMINE_AVGLOAD
    kstat_ctl_t *kc;
    kstat_t *ks;
    kstat_named_t *kn;

    kc = kstat_open();
    if(kc!=NULL)
    {
        ks = kstat_lookup(kc, "unix", -1, "system_misc");
        if(ks!=NULL)
        {
            if(kstat_read(kc, ks, 0) != -1)
            {
                // get load averages of last 1, 5 or 15 minutes
                /*
                kn = (kstat_named_t*)kstat_data_lookup(ks, "avenrun_1min");
                kn = (kstat_named_t*)kstat_data_lookup(ks, "avenrun_15min");
                */
                kn = (kstat_named_t*)kstat_data_lookup(ks, "avenrun_5min");
                if (kn)
                    avgload = (int)(loaddouble(kn->value.ul));
                else
            	    prt_error("send_ncpus", "kstat_data_lookup failed");
            }
            else
                prt_error("send_ncpus", "kstat_read failed");
        }
        else
            prt_error("send_ncpus", "kstat_lookup failed");
    }
    else
        prt_error("send_ncpus", "kstat_open failed");
    kstat_close(kc);
#endif // DETERMINE_AVGLOAD
#elif defined(__linux)
#define LINELEN 128
    char Line[LINELEN+1];
    ifstream Info;
    
    // read "/proc/cpuinfo" and count the "processor" lines
    Info.open("/proc/cpuinfo");
    if (!Info)	// could not open, but
    { 
	Ncpus=1;    // there must be a CPU if this is running :-)
	prt_error("send_ncpus", "cannot open\"/proc/cpuinfo\"");
    }
    else
    {
	Ncpus=0;
	while (Info.getline(Line, LINELEN, '\n').good())
	    if (!strncmp(Line, "processor", 9)) Ncpus++;
	Info.close();
    }
#undef LINELEN    
#ifdef DETERMINE_AVGLOAD
    float Load5;    // 5-minute load average
    
    Info.open("/proc/loadavg");
    if (!Info)
	prt_error("send_ncpus", "cannot open\"/proc/loadavg\"");
    else
    {
	Info>>Load5>>Load5; // skip 1-minute average
	Info.close();
	avgload=(int)Load5;
    }
#endif // DETERMINE_AVGLOAD
#elif defined(__sgi)
    Ncpus = sysconf(_SC_NPROC_ONLN);
    #ifdef DETERMINE_AVGLOAD
    FILE *Uptimepipe;
    
    // Jan Hungershoefer's idea to parse the output of the 'uptime'
    // command on SGI-s to get the average load. The user can
    // specify the location of the command in the UPTIME env var.
    char *Uptime = getenv("UPTIME");
    // open pipe to uptime cmd (let's hope it is available)
    if((Uptimepipe = popen(Uptime==NULL? "uptime": Uptime, "r")) != NULL)
    {
        int i;
        float ld;
#ifdef BSD41
        i = fscanf(Uptimepipe,"%*[^l] load average: %*f %f", &ld);
#else
        i = fscanf(Uptimepipe,"%*[^l] load average: %*f, %f", &ld);
#endif /* BSD41 */
        avgload = (int)ld;
        i = pclose(Uptimepipe); // ignore errors when closing
    }
    #endif  // DETERMINE_AVGLOAD
#endif	// architecture-dependent section
    
    cout << "#CPU="<<Ncpus<<", average load here " << avgload << endl;
    // subtract average load from number of CPUs
    Ncpus -= avgload;
    // make sure Ncpus stays >= 1 because one slave already there
    // user should remove host from virtual machine if overcrowded...
    Ncpus = (Ncpus>1)?Ncpus:1;
    int Result=pvm_initsend(PvmDataDefault);
    if (Result<0)
    {
        prt_error("send_ncpus()", Result);
        return(Result);
    }
    Result=pvm_pkint(&Ncpus, 1, 1);
    if (Result<0)
    {
        prt_error("send_ncpus()", Result);
        return(Result);
    }
    Result=pvm_send(Mastertid, CPUCNT);
    if (Result<0)
    {
        prt_error("send_ncpus()", Result);
        return(Result);
    }
    return(0);
}
// END of send_ncpus()

/* wait_master(): when invoked by a slave, this method
 * waits for messages from the master. If Tag==ANY (the default), 
 * then the method returns when a message is received.
 * Otherwise, it returns only if a message with the 
 * specified Tag arrived or the master died (MASTER_EXIT).
 * The message is not read, just peeked into.
 * Return value: the buffer ID or <0 on error (-1 on master's death).
 * The value in Tag is replaced by the actual tag.
 */
int Pvmtask_::wait_master(int& Tag) const
{
    if (!is_slave()) return(-1);
    
    // block and wait for any message
    int Result=0, Bufid=0;
    while(PvmOk==(Result=pvm_pstat(Mastertid)))
    {
	// check if anything relevant has arrived
	Bufid=pvm_probe(Mastertid, (Tag==ANY? -1: Tag));
	if (!Bufid)
	{
	    sleep(1); continue;	// wait a bit then check again
	}
	if (Bufid<0)	// error
	{
	    prt_error("wait_master() [pvm_probe]", Bufid);
	    return(Bufid);
	}
	
	// peek into the message
	int Bytes=0, Msgtag=0, Mtid=0;
	Result=pvm_bufinfo(Bufid, &Bytes, &Msgtag, &Mtid);
	if (Result<0)	// error
	{
	    prt_error("wait_master() [pvm_bufinfo]", Result);
	    return(Result);
	}
	
	// accept messages having the right tag
	if (Tag==ANY || Tag==Msgtag)
	{
	    Tag=Msgtag; return(Bufid);
	}
	
	// ignore all other messages
    }	    // while
    
    // master not available
    prt_error("wait_master() [pvm_probe]", Result);
    return(Result);
}
// END of wait_master()

/* prt_error(Methodname, Pvminfo): prints PVM error messages to cerr.
 * Not all PVM error conditions are implemented: RTFM for the missing ones.
 * The format is: "\n? Pvmtask_::Methodname (M|S)[tid]@host: PVM error (message)\n"
 * prt_error(Methodname, Infostr): prints general error messages to cerr.
 * The format: "\n? Pvmtask_::Methodname (M|S)[tid]@host: Infostr\n"
 * (M|S) is for master/slave distinction.
 * Private
 */
void Pvmtask_::prt_error(const char *Methodname, int Pvminfo) const
{
    if (Pvminfo>=PvmOk) return;	// no error
    
    cerr<<"\n? Pvmtask_::"<<Methodname<<" TID="
	<<(is_master()? 'M': 'S')<<Idstr<<": PVM error(";
    switch(Pvminfo)
    {
	case PvmBadParam: cerr<<"bad parameter"; break;
	case PvmNoData: cerr<<"read past end of buffer"; break;
	case PvmNoHost: cerr<<"unknown host"; break;
	case PvmNoFile: cerr<<"cannot find executable"; break;
	case PvmNoMem: cerr<<"out of memory"; break;
	case PvmBadMsg: cerr<<"cannot decode message"; break;
	case PvmSysErr: cerr<<"daemon not responding"; break;
	case PvmNoBuf: cerr<<"no current buffer"; break;
	case PvmNoSuchBuf: cerr<<"bad message ID"; break;
	case PvmHostFail: cerr<<"host failed"; break;
	case PvmNoParent: cerr<<"no parent task"; break;
	case PvmDSysErr: cerr<<"daemon system error"; break;
	case PvmOutOfRes: cerr<<"out of resources"; break;
	case PvmNoTask: cerr<<"nonexistant task"; break;
	default: cerr<<"code "<<Pvminfo<<", see manual"; break;
    }
    cerr<<")\n"<<flush;
}

void Pvmtask_::prt_error(const char *Methodname, const char *Infostr) const
{
    cerr<<"\n? Pvmtask_::"<<Methodname<<" TID="
	<<(is_master()? 'M': 'S')<<Idstr<<": "<<Infostr<<endl;
}
// END of prt_error()

// ==== END OF METHODS Pvmtask.c++ ====
#endif  /* USE_PVM */
