#ifdef USE_PVM
#ifndef PVMTASK_CLASS
#define PVMTASK_CLASS

// ==== PROJECT DRAGON: HEADER Pvmtask.h ====

/* Task management on the Parallel Virtual Machine (PVM).
 * PVM is free software originally developed at the
 * Oak Ridge National Laboratory, Oak Ridge, Tennessee, USA.
 */

// 5-Aug-2000. Andras Aszodi

// ---- PVM HEADER ----

#include <iostream.h>
// this is because of a bug in PVM 3.4:
#include <stdio.h> 

#include "pvm3.h"

// ---- MODULE HEADERS ----

#include "Params.h"
#include "Sigproc.h"

// ==== CLASSES ====

/* Class Pvmtask_: this class manages a set of DRAGON tasks on
 * the Parallel Virtual Machine. A master task spawns a number
 * of slave tasks which perform the simulations while the master
 * provides a front-end to the outside world and sends updated
 * parameter values and run commands to the slaves. The paradigm
 * is similar to the multiprocess option implemented in
 * "Sigproc" but here the slaves live as long as the master.
 * There is one object per process.
 */
class Pvmtask_
{
    // enums
    public:
    
    // Filetags_: data file contents (cf. "Params.h")
    enum Filetags_ {ALN, PHO, VOL, ACD, SIM, RESTR, SSTR, ACC, HOM};
    
    // Msgtags_: PVM message tags for process status check.
    enum Msgtags_ {SLAVE_READY=6500, SLAVE_DONE, SLAVE_RUNNING, PARAMS, RUN, ANY, CPUCNT};

    // Pvmstat_: decides whether PVM is running and who's the boss
    enum Pvmstat_ {NO_PVM, MASTER, SLAVE};
    
    // data
    private:
    
    /* Slave_: this classlet stores the status of a slave. */
    struct Slave_
    {
	int Tid;    // the task ID
	int Hosttid;	// the TID of the slave's PVM daemon
	int Jobs;  // the jobs (>0) done on the slave: 0 if idle, -1 if param sync is needed
	
	Slave_(int T=0, int H=0, int J=-1): 
	    Tid(T), Hosttid(H), Jobs(J) {}
    };
    
    // needed to remember information about every PVM host
    // identification by host's TID 
    // currently only information stored is the number of CPUs of that host
    struct Host_
    {
        int Tid;        // the TID of the host's PVM demon
        int Ncpu;       // the number of processors of the host
        Host_(int T=0, int N=1):
          Tid(T), Ncpu(N) {}
    };

    static unsigned int Objno;	// only 1 object per pgm is allowed
    
    Host_ *Hosts;   // array of Host_ (NULL in slaves), dynamically allocated
    int Nhosts;     // size of array Hosts
    Slave_ *Slaves; // slave data in master (NULL in slaves)
    int Slaveno;    // no. of slaves in master (must not be equal to no. of hosts, 0 in slaves)
    char *Idstr; // "[Tid]@hostname"
    char *Slavexec; // slave executable name in master, NULL in slaves
    int Tid, Mastertid;	// task ID and master's ID (-1 for master)
    Pvmstat_ Pvmstat;	// master/slave status (or no PVM)
    
    // methods
    public:
    
	// initialisation
	    
    // Inits the object to empty.
    Pvmtask_();
	
    // enrol_pvm(Slexenm):
    // enrols the process in PVM and finds out if it runs as a master or slave.
    // The optional Slexenm argument specifies the slave executable name.
    void enrol_pvm(const char *Slexenm=NULL);
    
	// destructor
    /* When the master exits, it attempts to kill all remaining slaves
     * before dying. This assumes that the calling object was global
     * and dies only with the program.
     */
    ~Pvmtask_();
    
	// access
    int is_master() const { return(Pvmstat==MASTER); }
    int is_slave() const { return(Pvmstat==SLAVE); }
    int no_pvm() const { return(Pvmstat==NO_PVM); }
    int tid() const { return(Tid); }
    const char* id_str() const { return(Idstr); }
    
	// task management
    /* spawn_slaves(): attempts to launch one slave on each node
     * of the virtual machine. Checks old slaves and virtual machine
     * status, and spawns slaves on the nodes which have none.
     * The slave executable's name is Slavenm, no command-line
     * arguments are passed: if left unspecified, then the
     * internal Slavexec string is used. Can be invoked by the master only.
     * Return value: if PVM is not available, then <0 is returned, 
     * otherwise the no. of newly created slaves is returned.
     */
    int spawn_slaves(const char *Slavenm=NULL);

	// signals
    /* signal_pvm(): catches signals in the master task in PVM
     * and sends them on to the slaves. The slaves use the
     * "standard" signal_handler() function in "Sigproc".
     * The logic is that most signals terminate the slaves, 
     * whereupon the master will automatically exit; when
     * SIGINT <Ctrl-C> is caught, the slaves will be terminated.
     * NOTE: there must be a global Pvmtask_ called Pvmtask for this to work.
     */
    friend void signal_pvm(int Sigtype) throw(Sigexcept_);
	
    /* recv_ncpus(): for the master only
     * waits for the freshly spawned slave to send the 
     * number of cpus found on its node
     */
    int recv_ncpus(int Tid);
    /* send_ncpus(): for the slaves only (machine dependant)
     * determines the number of cpus on the machine
     */
    int send_ncpus();

	// Communication
    /* send_params(): first sends a string containing a modified
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
    int send_params(Params_& P);
    int send_params(Params_& P, const int *Newslaves, int Newno);

    /* recv_params(): in a slave process, this method should be
     * called when wait_master(PARAMS) returns successfully
     * from the listening. The PARAMS message holds a parameter
     * list string which will be read into P directly.
     * Note that PARAMS is in the active buffer and must be read
     * immediately after receiving with wait_master().
     * Return value: the length of the message.
     */
    int recv_params(Params_& P) const;
	
    /* recv_filestr(): receives one of the data files as a string,
     * indicated by the Tag. Waits until the message arrives
     * or the master dies. Returns a char ptr pointing
     * to a string containing the corresponding ASCII data, or ""
     * on error. Empty files also return a "". NOTE: the 
     * returned pointer points to a static area which will
     * be overwritten upon each invocation of the method.
     * Does nothing if not invoked in a slave.
     */
    char* recv_filestr(int Tag) const;
	
	// Job distribution
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
    int send_jobs(Params_& Params, unsigned int& Runno);
    
    /* recv_job(): when invoked by a slave immediately after wait_master()
     * returned with the news that a RUN message has been sent, this
     * method reads the message and returns the job number (>0)
     * or a negative int on error.
     */
    int recv_job() const;
    
    /* job_status(): when invoked by a slave, this method sends a
     * message with Tag to the master. The message consists of an
     * integer Num. If Tag==SLAVE_RUNNING, then it is the number of cycles
     * already done. If Tag==SLAVE_DONE, then it is the job number.
     * Return value: <0 indicates an error.
     */
    int job_status(int Tag, int Num) const;
    
    /* slave_ready(): this method tells the master that the slave
     * has received all the parameters and is ready to run the
     * simulations.
     * Return value: <0 on error.
     */
    int slave_ready() const;
    
    /* wait_master(): when invoked by a slave, this method
     * sends a SLAVE_READY message to the master and then
     * waits for messages from the master. If Tag==ANY (the default), 
     * then the method returns when a message is received.
     * Otherwise, it returns only if a message with the 
     * specified Tag arrived or the master died (MASTER_EXIT).
     * The message is not read, just peeked into.
     * Return value: the buffer ID or <0 on error (-1 on master's death).
     * The value in Tag is replaced by the actual tag.
     */
    int wait_master(int& Tag) const;
    
    // hidden methods
    private:
    
    int check_slaves();
    int send_files(const int* Slavetids, int Slno, Params_& P,
	const char *Pname, Filetags_ Tag);
    
    void prt_error(const char *Methodname, int Pvminfo) const;
    void prt_error(const char *Methodname, const char *Infostr) const;
    
    // forbidden methods
    Pvmtask_(const Pvmtask_&);
    Pvmtask_& operator=(const Pvmtask_&);
};
// END OF CLASS Pvmtask_

// ==== END OF HEADER Pvmtask.h ====
#endif	/* PVMTASK_CLASS */
#endif  /* USE_PVM */
