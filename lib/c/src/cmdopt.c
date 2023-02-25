/* ==== FUNCTIONS cmdopt.c ==== */

/* One-letter command line option processing. */

/* ANSI C, 27. Apr. 1998. Andris */

/* ---- HEADERS ---- */

/* NOTE: check where getopt() is declared on the SGI. */
#ifdef __linux__
#include <unistd.h>
#endif

#include "cmdopt.h"

/* ---- DEFINITIONS ---- */

#ifndef GETOPTHUH
#define GETOPTHUH '?'
#endif
#define OPTMAXLEN 63	/* 2*26+10 alphanum + '#' */
#define DESCRMAXLEN 32  /* arbitrary limit to option descriptor strings */

/* ---- TYPEDEFS ---- */

/* Argtype_: the command option argument type. CMDOPT_BOOL is a
 * Boolean value (stored as an int just like CMDOPT_INT) for switches
 * such as "-x" with no arguments. The rest are for int,  double and 
 * string arguments.
 */
typedef enum {CMDOPT_BOOL, CMDOPT_INT, CMDOPT_DBL, CMDOPT_STR} Argtype_;

/* Argval_: a union that stores the value of the option arguments. */
typedef union
{
    int i;	/* integer, double and string */
    double d;
    char *s;
}
Argval_;

/* Option_: a struct that describes a command-line option. */
typedef struct
{
    char Ch;	/* the option character */
    int Idx;	/* argv[] index, 0 if absent */
    Argtype_ Type;  /* type of argument */
    Argval_ Val;    /* the value */
    char Descr[DESCRMAXLEN+1];    /* descriptor string */
}
Option_ ;

/* ---- PROTOTYPES ---- */

static int token_type(char *Tok);
static int good_optchar(char Och);

/* ---- FILE-SCOPE VARIABLES ---- */

/* Cmdopts is a static array of Option_ struct-s and holds the necessary
 * information for each command-line option. Cmdoptno is its length.
 * Only one option array is needed in a program so Cmdopts and Cmdoptno have
 * local file scope and are accessed from the routines as globals.
 */
static Option_ Cmdopts[OPTMAXLEN];
static int Cmdoptno=0;

/* ==== FUNCTIONS ==== */

/* ---- Parse option string ---- */

/* parse_optstr(): constructs and returns the hidden Option_ array from the
 * string Cmdoptstr. The length of the array is also stored in a hidden location.
 * The option string is composed of the following tokens separated
 * by whitespaces:-
 * 
 * "x": stands for the Boolean command line option -x
 * "xYz" : the command line options -x -Y -z (may be grouped together)
 * "x%d<name>" : option -x expecting a mandatory integer argument
 *               which is described as "name" in the help string.
 * "x%f<name>" : as above, expecting a floating-point argument
 * "x%s<name>" : as above, expecting a string argument
 * 
 * 'x' can be the characters 'a'..'z', 'A'..'Z', '0'..'9', '#'.
 * Invalid options are ignored. If nothing was found then Cmdopts==NULL and 
 * Cmdoptno==0. To be called only once.
 */
void parse_optstr(char *Cmdoptstr)
{
    /* NOTE: there is a bizarre bug in the GNU C library (version 1) */
    /* as supplied with Linux 2.0.x: strtok() dumps core */
    /* if its input string is a parameter of a function. */
    /* Workaround: make a local copy of Cmdoptstr. */
    /* To activate the workaround, define GLIBC_STRTOK_BUG */
    /* on the cc command line. */
    /* 27-Mar-1998. */

    const char *Wsp=" \t\n";	/* whitespaces */
    char *Tok=NULL;

    #ifdef GLIBC_STRTOK_BUG
    char *Cmdoptstr_copy=NULL;
    #endif

    int i, Atyp;
    
    /* check virginity */
    if (Cmdoptno)
    {
	fputs("? parse_optstr(): not first call, ignored\n", stderr);
	return;
    }
    
    /* get first token */
    #ifdef GLIBC_STRTOK_BUG
    Cmdoptstr_copy=strdup(Cmdoptstr);  /* workaround: local copy */
    if (NULL==(Tok=strtok(Cmdoptstr_copy, Wsp)))
    {
      free(Cmdoptstr_copy);
      return;
    }
    #else
    if (NULL==(Tok=strtok(Cmdoptstr, Wsp))) return; /* no tokens at all */
    #endif

    do	/* process tokens */
    {
	Atyp=token_type(Tok);	/* type check: bad syntax is <0 */
	if (Atyp<0)
	{
	    fprintf(stderr, "? parse_optstr(): Bad token \"%s\"\n", Tok);
	    continue;
	}
	
	if (Atyp==CMDOPT_BOOL)	/* Boolean option(s) */
	    for (i=0; i<strlen(Tok); i++)
	    {
		if (!good_optchar(Tok[i]))
		{
		    fprintf(stderr, 
			"? parse_optstr(): Boolean option \'%c\' is duplicate or invalid\n",
			Tok[i]);
		    continue;
		}
		
		Cmdopts[Cmdoptno].Ch=Tok[i]; Cmdopts[Cmdoptno].Type=CMDOPT_BOOL;
		Cmdopts[Cmdoptno].Idx=0; Cmdopts[Cmdoptno].Descr[0]='\0';
		Cmdoptno++;
	    }
	else	    /* mandatory argument option: one per token */
	{
	    if (!good_optchar(Tok[0]))
	    {
		fprintf(stderr, 
		    "? parse_optstr(): Arg option \'%c\' is duplicate or invalid\n",
		    Tok[0]);
		continue;
	    }
		
	    Cmdopts[Cmdoptno].Ch=Tok[0]; Cmdopts[Cmdoptno].Type=(Argtype_)Atyp;
	    Cmdopts[Cmdoptno].Idx=0;
	    if (strlen(Tok+4)>DESCRMAXLEN)
	    {
	      fprintf(stderr,
		      "\n? parse_optstr(): Option descriptor \"%s\" longer than %d chars, truncated\n",
		      Tok+4,DESCRMAXLEN);
	      strncpy(Cmdopts[Cmdoptno].Descr,Tok+4,DESCRMAXLEN);
	      Cmdopts[Cmdoptno].Descr[DESCRMAXLEN]='\0';
	    }
	    else
	      strcpy(Cmdopts[Cmdoptno].Descr,Tok+4);
	    Cmdoptno++;
	}
    }
    while(NULL!=(Tok=strtok(NULL, Wsp)));
    
    #ifdef GLIBC_STRTOK_BUG
    free(Cmdoptstr_copy);
    #endif
}
/* END of parse_optstr() */

/* token_type(): parses the token Tok and returns -1 if it is not
 * of the form defined for the option string (see above) or returns
 * the argument type (see Argtype_ ). For mandatory argument tokens, 
 * the first '>' is converted to '\0'.
 */
static int token_type(char *Tok)
{
    /* arg token? */
    if (NULL!=strchr(Tok, '%'))
    {
	if (Tok[1]!='%' || Tok[3]!='<' || 
		strlen(Tok)<5 || Tok[strlen(Tok)-1]!='>')
	    return(-1);	    /* malformed arg token */
	Tok[strlen(Tok)-1]='\0';    /* chop off '>' */
	switch (Tok[2])
	{
	    case 'd': return(CMDOPT_INT);	/* get type */
	    case 'f': return(CMDOPT_DBL);
	    case 's': return(CMDOPT_STR);
	    default: return(-1);	/* malformed */
	}
    }
    else return(CMDOPT_BOOL);  /* non-arg Boolean token */
}
/* END of token_type() */

/* good_optchar(): checks whether the option char Och is acceptable.
 * It should be alphanumeric or '#' and it should not be duplicated.
 * Duplication is checked by keeping a static record of characters tested
 * so far (that's why calling parse_optstr() twice in a program would
 * screw up everything). 
 * Return value: 0 if Och is "bad" or duplicated, non-0 otherwise.
 */
static int good_optchar(char Och)
{
    static char Found[OPTMAXLEN+1];
    static int Len=0;
    
    /* unacceptable char */
    if (!isalnum(Och) && Och!='#') return(0);
    
    if (!Len) Found[0]='\0';	/* "init"  array for the first time */
    
    if (NULL!=strchr(Found, Och)) return(0);	/* has been seen already */
    Found[Len++]=Och; Found[Len]='\0';
    return(Len);
}
/* END of good_optchar() */

/* ---- Parse command line ---- */

/* get_options(): processes the command line (argc,argv) to find options
 * encoded in the hidden array Cmdopts (length Cmdoptno). If an option is found then
 * the corresponding Cmdopts[] member is updated to show which argv[] member
 * contained it and if it had an argument, then it is attempted to be
 * stored as well. Options not found or having malformed arguments will
 * be set as "absent".
 * Return value: the argv[] index of the first non-option member.
 * The index is multiplied by -1 if an error occurred.
 */
int get_options(int argc, char *const *argv)
{
    extern int optind;	    /* next argv[] index (set by getopt()) */
    extern int opterr;    /* suppress error messages by setting to 0 */
    extern int optopt;	    /* bad arg character stored in this */
    extern char *optarg;    /* option argument ptr set by getopt() */
    signed char Opt;
    char Err, *Endconv, *Cmdoptstr;
    int i, Errsign=1;
    long Ltmp;
    double Dtmp;
    
    /* check option array, build getopt() control string */
    Cmdoptstr=(char *) calloc(2*Cmdoptno+1, sizeof(char));
    for (i=0; i<Cmdoptno; i++)
    {
	/* a colon is appended after option chars w/ arguments */
	sprintf(Cmdoptstr+strlen(Cmdoptstr), 
	    ((Cmdopts[i].Type!=CMDOPT_BOOL)? "%c:":"%c"), Cmdopts[i].Ch);
	Cmdopts[i].Idx=0;	/* set all to 'absent' */
	
	switch (Cmdopts[i].Type)	/* zero the values */
	{
	    case CMDOPT_BOOL:
	    case CMDOPT_INT:
	    Cmdopts[i].Val.i=0; break;
	    
	    case CMDOPT_DBL:
	    Cmdopts[i].Val.d=0.0; break;
	    
	    case CMDOPT_STR:
	    Cmdopts[i].Val.s=NULL; break;
	}
    }

    /* process the command line */
    opterr=0;
    while (EOF!=(Opt=getopt(argc, argv, Cmdoptstr)))
    {
	if (Opt==GETOPTHUH)	/* unrecognised option */
	{
	    fprintf(stderr, "\n? %s: %s option %c\n", argv[0], 
		(strchr(Cmdoptstr, optopt))? 
		    "Bad argument for": "Unknown", optopt);
	    Errsign=-1;
	    continue;	/* skip */
	}
	
	/* get option index (Opt must be there) */
	for (i=0; i<Cmdoptno; i++)
	    if (Cmdopts[i].Ch==Opt) break;
	    
	Cmdopts[i].Idx=optind;   /* save argv[] index: corrected 19-Apr-1998 */
	Err=0;
	switch (Cmdopts[i].Type)
	{
	    case CMDOPT_INT:   /* integer (orig. long) conversion */
	    Ltmp=strtol(optarg, &Endconv, 10);
	    if (*Endconv!='\0') Err=Opt;
	    else Cmdopts[i].Val.i=(int)Ltmp;
	    break;
	    
	    case CMDOPT_DBL:   /* double conversion */
	    Dtmp=strtod(optarg, &Endconv);
	    if (*Endconv!='\0') Err=Opt;
	    else Cmdopts[i].Val.d=Dtmp;
	    break;
	    
	    case CMDOPT_STR:   /* simply store the ptr */
	    Cmdopts[i].Val.s=optarg;
	    break;
	    
	    case CMDOPT_BOOL:
	    default: ; break;	/* do nothing */
	}
	
	/* error during conversion? set offending option to absent */
	if (Err)
	{
	    fprintf(stderr, "\n? %s: Malformed argument for option %c, default used\n", 
		argv[0], Err);
	    Cmdopts[i].Idx=0;  /* reset */
	    Errsign=-1;
	    continue;
	}
    }
    
    free(Cmdoptstr);
    return(Errsign*optind);
}
/* END of get_options() */

/* ---- Return values ---- */

/* optval_[bool,int,dbl,str](): these functions query the Cmdopts[] array
 * after the command line has been processed. Och is the option char
 * and the value is returned in *Val. If the option char is invalid
 * then a warning is printed. If the option was not present on the
 * command line then *Val is not updated and the functions return 0, 
 * otherwise the Idx field from the corresponding Cmdopts record is returned.
 * Val may be set to NULL.
 */
int optval_bool(char Och)
{
    int i;
    
    for (i=0; i<Cmdoptno && Cmdopts[i].Ch!=Och; i++);
    if (i>=Cmdoptno)
    {
	fprintf(stderr, "? optval_bool(): invalid option \'%c\'\n", Och);
	return(0);
    }
    if (Cmdopts[i].Type!=CMDOPT_BOOL)
    {
	fprintf(stderr, "? optval_bool(): option \'%c\' not Boolean\n", Och);
	return(0);
    }
    return(Cmdopts[i].Idx);
}

int optval_int(char Och, int *Val)
{
    int i;
    
    for (i=0; i<Cmdoptno && Cmdopts[i].Ch!=Och; i++);
    if (i>=Cmdoptno)
    {
	fprintf(stderr, "? optval_int(): invalid option \'%c\'\n", Och);
	return(0);
    }
    if (Cmdopts[i].Type!=CMDOPT_INT)
    {
	fprintf(stderr, "? optval_int(): option \'%c\' not integer\n", Och);
	return(0);
    }
    if (Val!=NULL && Cmdopts[i].Idx) *Val=Cmdopts[i].Val.i;
    return(Cmdopts[i].Idx);
}

int optval_dbl(char Och, double *Val)
{
    int i;
    
    for (i=0; i<Cmdoptno && Cmdopts[i].Ch!=Och; i++);
    if (i>=Cmdoptno)
    {
	fprintf(stderr, "? optval_dbl(): invalid option \'%c\'\n", Och);
	return(0);
    }
    if (Cmdopts[i].Type!=CMDOPT_DBL)
    {
	fprintf(stderr, "? optval_dbl(): option \'%c\' not double\n", Och);
	return(0);
    }
    if (Val!=NULL && Cmdopts[i].Idx) *Val=Cmdopts[i].Val.d;
    return(Cmdopts[i].Idx);
}

int optval_str(char Och, char **Val)
{
    int i;
    
    for (i=0; i<Cmdoptno && Cmdopts[i].Ch!=Och; i++);
    if (i>=Cmdoptno)
    {
	fprintf(stderr, "? optval_str(): invalid option \'%c\'\n", Och);
	return(0);
    }
    if (Cmdopts[i].Type!=CMDOPT_STR)
    {
	fprintf(stderr, "? optval_str(): option \'%c\' not string\n", Och);
	return(0);
    }
    if (Val!=NULL && Cmdopts[i].Idx) *Val=Cmdopts[i].Val.s;
    return(Cmdopts[i].Idx);
}

/* END of optval_[...] */

/* opt_defval(): this routine provides an alternative to the optval_[...]
 * functions. Och is the option character, Val is a void pointer to
 * the variable which should hold the option value, Defval is a
 * void pointer to the default value. If the option is found, 
 * then *Val is set to its value if Val!=NULL. If not found or on error, 
 * then *Val will be set to *Defval (to the default) if both
 * Val and Defval !=NULL. 
 * Be careful, there is NO POINTER TYPE CHECKING!
 * Return value: 0 on error or if the option was not found, 
 * the option array index otherwise.
 */
int opt_defval(char Och, void *Val, const void *Defval)
{
    int i;
    int *Vip;	/* temp ptrs with types */
    double *Vdp;
    char **Vsp;
    
    for (i=0; i<Cmdoptno && Cmdopts[i].Ch!=Och; i++);
    if (i>=Cmdoptno)
    {
	fprintf(stderr, "\n? opt_defval(): invalid option \'%c\'\n", Och);
	return(0);
    }
    if (Val!=NULL)
    {
	switch(Cmdopts[i].Type)
	{
	    case CMDOPT_BOOL:
		Vip=(int *)Val;
		*Vip=(Cmdopts[i].Idx)? 1: 0;
	    break;
	    case CMDOPT_INT:
		Vip=(int *)Val;
		if (Cmdopts[i].Idx) *Vip=Cmdopts[i].Val.i;
		else if (Defval!=NULL) *Vip=*((int*)Defval);
	    break;
	    case CMDOPT_DBL:
		Vdp=(double *)Val;
		if (Cmdopts[i].Idx) *Vdp=Cmdopts[i].Val.d;
		else if (Defval!=NULL) *Vdp=*((double*)Defval);
	    break;
	    case CMDOPT_STR:
		Vsp=(char **)Val;
		if (Cmdopts[i].Idx) *Vsp=Cmdopts[i].Val.s;
		else if (Defval!=NULL) *Vsp=*((char**)Defval);
	    break;
	}
    }
    return(Cmdopts[i].Idx);
}
/* END of opt_defval() */

/* ---- Help string ---- */

/* opt_helpstr(): generates a "help string" from the option list in
 * Cmdopts[]. Boolean options are collected together, the "argumented"
 * options are listed separately. If -x and -y are switches, -i expects
 * an integer option and -D a double, then the following string will
 * be returned: "[-xy] [-i<name>] [-D<name>]" where "name" stands
 * for the string stored in the Descr field. Space for the help string is
 * allocated within.
 */
char *opt_helpstr(void)
{
    char *Bs, *As;	/* strings for Boolean and Argument options */
    char *Hs;	/* the help string */
    int Alen, Blen, i;
    
    /* allocate big chunks */
    Bs=(char *) calloc(Cmdoptno+5, sizeof(char));
    As=(char *) calloc(10*Cmdoptno+3, sizeof(char));
    Alen=Blen=0;
    
    /* scan all options */
    for (i=0; i<Cmdoptno; i++)
	if (Cmdopts[i].Type==CMDOPT_BOOL)
	{
	    sprintf(Bs+Blen, (Blen)? "%c": "[-%c", Cmdopts[i].Ch);
	    Blen=strlen(Bs);
	}
	else
	{
	    sprintf(As+Alen, "[-%c %s] ", Cmdopts[i].Ch, Cmdopts[i].Descr);
	    Alen+=strlen(Cmdopts[i].Descr)+6;
	}
    
    if (Blen) { strcat(Bs, "] "); Blen+=2; }	/* close bracket */
    if (Alen) { As[Alen-1]='\0'; Alen--; }  /* chop off last space */
    
    /* join the two together */
    Hs=(char *) calloc(Alen+Blen+2, sizeof(char));
    strcpy(Hs, Bs); strcat(Hs, As);
    free(Bs); free(As);
    return(Hs);
}
/* END of opt_helpstr() */

#undef OPTMAXLEN 

/* ==== END OF FUNCTIONS cmdopt.c ==== */
