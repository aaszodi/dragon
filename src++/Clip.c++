// ==== PROJECT DRAGON: METHODS Clip.c++ ====

/* The command-line interpreter. */

// SGI C++ 6.2 (o32) , IRIX 6.2, 21. Aug. 1996. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <string.h>
#include <ctype.h>
#include <unistd.h>	// for isatty() and std fildes nos
#include <signal.h>
#include <iostream.h>
#include <iomanip.h>
#include <strstream.h>

// ---- MODULE HEADER ----

#include "Clip.h"

// ---- DEFINITIONS ----

/* The standard I/O file descriptors should be defined in
 * "unistd.h". If not, then we give them here
 */
#ifndef STDIN_FILENO
#define	STDIN_FILENO	0
#endif
#ifndef STDOUT_FILENO
#define STDOUT_FILENO	1
#endif

// ==== Clip_ METHODS ====

/* put_prompt(): writes the prompt string to cout.
 * Appends Cmdlevel '>' characters after the prompt string.
 */
void Clip_::put_prompt() const
{
    cout<<'\n'<<Prompt;
    for (int p=0; p<Cmdlevel; p++) cout.put('>');
    cout<<' '<<flush;
}

/* get_command(): tries to obtain the next command from the file
 * named Cmd (cin if Cmd==NULL). Most commands
 * can be dealt with internally. If a "r[un] x"
 * command is seen which means "run run_func() x times
 * with the current parameters" then x>0 is returned. 
 * run_func() takes an unsigned int specifying the number of cycles to
 * be executed. It is supposed to return 0 at the end of a
 * successful run or the value of a signal caught inside.
 * Must catch SIGINT (Ctrl-C) for run interrupts.
 * If "q[uit]" has been seen then -1 is returned, meaning "finish up everything".
 * If a script was processed and the end has been reached then
 * get_command() returns with 0. 
 */
int Clip_::get_command(const char *Cmdfnm, Runfnc_ run_func)
{
    static const unsigned int LINELEN=132, CMDLEN=20, MAX_CMDLEVEL=16;
    static const int RET_OK=0, RET_QUIT=-1;	// return values: OK and quit
    
    static char Line[LINELEN], Fname[LINELEN], Cmd[CMDLEN];
    static istrstream Instr(Line, LINELEN);
    
    ifstream In;
    char *Shell;
    int Interact=0, Cycno=1, Retval=RET_OK;
    
    Cmdlevel++;	    // new level of recursion

    // open the command file Cmdfnm or read from cin
    if (Cmdfnm!=NULL && strlen(Cmdfnm))
    {
	In.open(Cmdfnm);	// open for input
	if (!In)
	{
	    cerr<<"\n? Clip_::get_command(): Cannot open command file \""<<Cmdfnm<<"\"\n";
	    Cmdlevel--;
	    return(RET_OK);
	}
	else cout<<"Command script \""<<Cmdfnm<<"\", level "<<Cmdlevel<<endl;
    }
    else In.attach(STDIN_FILENO);	// use std input
    
    // interactive mode if In and cout are both associated with a terminal
    Interact=isatty(In.rdbuf()->fd()) && isatty(STDOUT_FILENO);
    if (Interact) cout<<"Interactive mode (press \'h\' for help)\n";
    In.tie(&cout);	// output flushed before input
    
    while (Retval==RET_OK && (put_prompt(), memset(Line, '\0', LINELEN), 
	    In.getline(Line, LINELEN, '\n').good()))
    {
	if (!strlen(Line) || Line[0]=='#')  // skip empty and comment
	    continue;
	
	Instr.seekg(0); Instr.clear(); // attempt interpretation
	Instr>>setw(CMDLEN)>>Cmd;
	if (Instr.bad() || Instr.fail())
	{
	    cerr<<"\n? Clip_::get_command(): Cannot read\n";
	    continue;
	}
	
	// parameter names start with uppercase
	if (isupper(Cmd[0]))
	{
	    Instr.seekg(0); Instr.clear();  // set good bit again,13-Apr-1998
	    Instr>>Params;  // send for parameter input
	    continue;
	}
	
	// only the first char of the command is significant
	switch(Cmd[0])
	{
	    case 'c':	// c[ommand] <file>
	    
	    if (Cmdlevel>=MAX_CMDLEVEL)	    // nested calls too deep
	    {
		cerr<<"\n? Clip_::get_command(): Only "<<MAX_CMDLEVEL<<" nested calls are allowed\n";
		break;
	    }
	    
	    // get new command file name: don't nest interactive modes
	    Instr>>setw(LINELEN)>>Fname;
	    if (Interact && (Fname==NULL || !strlen(Fname)))
		break;
	    
	    // from now on, process Fname (recursive)
	    Retval=get_command(Fname, run_func);
	    if (Retval<0) Retval=RET_OK; // low-level 'q' is swallowed
	    break;
	    
	    case 'd':	// default values
	    Params.reset_default();
	    cout<<"Parameters reset to default\n";
	    break;
	    
	    case 'h':	// h[elp]
	    cout<<"c[ommand] <file>: execute commands in <file>\n";
	    cout<<"d[efault]: reset all parameters to default\n";
	    cout<<"h[elp]: print this help\n";
	    cout<<"l[ist]: list all parameters to stdout\n";
	    cout<<"l[ist] <Param>: list parameter <Param> to stdout\n";
	    cout<<"o[s]: OS shell\n";
	    cout<<"p[aram] <file>: read parameters from <file>\n";
	    cout<<"q[uit]: quit\n";
	    cout<<"r[un] <int>: run "<<Prompt<<" int times (default 1)\n";
	    cout<<"s[ave] <file>: save parameters to <file>\n";
	    cout<<"<Param> <value>: set parameter <Param> to <value>\n";
	    break;
	    
	    case 'l':	// l[ist] a parameter 
	    Instr>>setw(LINELEN)>>Fname;
	    if (Fname[0]=='\0')	    // no param name
		cout<<Params;	// list all parameters
	    else
		Params.list_param(Fname);   // list param in Fname to cout
	    break;
	    
	    case 'o':	// invoke operating system shell
	    Shell=getenv("SHELL");
	    if (Shell==NULL)
		cerr<<"\n? Sorry, OS shell is unavailable\n";
	    else
	    {
		cout<<"\nType \'exit\' to return to "<<Prompt<<endl;
		system(Shell);
	    }
	    break;
	    
	    case 'p':	// p[aram] <file>
	    Instr>>setw(LINELEN)>>Fname;
	    if (Fname[0]=='\0')	    // no file
		cerr<<"\n? Please specify parameter file\n";
	    else
		Params.read_file(Fname);    // try to read from Fname
	    break;
	    
	    case 'q':	// quit
	    if (Interact)
	    {
		cout<<"Do you really wish to exit "<<Prompt<<" (y/n)? ";
		In.getline(Fname, LINELEN, '\n');
		Retval=(Fname[0]=='y' || Fname[0]=='Y')? RET_QUIT: RET_OK;
	    }
	    else Retval=RET_QUIT;
	    break;
	    
	    case 'r':	// r[un] <int>
	    Instr>>Cycno;   // try to get no. of cycles
	    if (!Cycno || !Instr) Cycno=1;
	    Cycno=abs(Cycno); 
	    cout<<"\nRun "<<Prompt;
	    if (Cycno==1) cout<<" once\n";
	    else cout<<" "<<Cycno<<" times\n";
	    Retval=run_func(Cycno); // launch simulation: Retval>=0
	    
	    if (Retval==SIGINT)
		Retval=RET_OK;   // Ctrl-C caught, pretend OK status
	    break;
	    
	    case 's':	// s[ave] parameters to file
	    Instr>>setw(LINELEN)>>Fname;
	    if (Fname[0]=='\0')	    // no output file name
		cout<<Params;	// list all parameters to stdout like l[ist]
	    else
	    {
		ofstream Poutf(Fname);
		if (!Poutf)
		{
		    cerr<<"\n? Cannot save parameters to file \""<<Fname<<"\"\n";
		    break;
		}
		Poutf<<Params;
		Poutf.close();
	    }
	    break;
	    
	    default:
	    cerr<<"\n? Unrecognised command\n";
	    break;
	}	// switch
	
    }	    // while
    
    if (In.rdbuf()->fd()!=STDIN_FILENO) In.close();   // close cmd file if not cin
    In.tie(NULL);	// undo tie
    Cmdlevel--;	    // return to previous command level
    return(Retval);	// end of commands
}
// END of get_command()

// ==== END OF METHODS Clip.c++ ====
