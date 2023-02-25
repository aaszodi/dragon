/* ==== FUNCTIONS tstamp.c ==== */

/* Useful formatted time routines. */

/* ANSI C, IRIX 6.2, 2. Sept. 1996. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <sys/times.h>

/* ---- MODULE HEADER ---- */

#include "tstamp.h"

/* ---- DEFINITIONS ---- */

/* In System V, the struct tms holds the time measurements in
 * clock tick units: CLK_TCK units make up a second. If this is
 * not defined, we assume the time measurements are in 1/100 seconds.
 * Apparently, this is POSIX but not pure ANSI
 */
#define _POSIX_SOURCE
#include <limits.h>
#ifndef CLK_TCK
#define CLK_TCK 100
#endif
#undef _POSIX_SOURCE

#define TSTAMP_LEN 30	/* length of time stamp string */
#define TIMESTR_LEN 80	/* length of "day-hr-min-.." string */

/* ---- LOCAL VARIABLES ---- */

static struct tms Start, Stop;	    /* process time structures, see <sys/times.h> */

/* ==== FUNCTIONS ==== */

/* time_stamp(): returns a pointer to an internal static string
 * that contains a time string in the following format:-
 * "Thu 02-Jun-1994 18:24:23". The string is evaluated at the
 * time of the call and is overwritten between calls.
 */
char *time_stamp(void)
{
    static char Tstr[TSTAMP_LEN];
    
    time_t Now;
    struct tm *Tm;
    
    Now=time(NULL);
    Tm=localtime(&Now);
    strftime(Tstr, TSTAMP_LEN, "%a %d-%b-%Y %X", Tm);
    return(Tstr);
}
/* END of time_stamp */

/* greeting(): returns a pointer to an internal static string
 * that contains a greeting appropriate for the time of the day:
 * "Good morning" from 6am till noon, "Good afternoon" from 12:01pm
 * until 6pm, "Good evening" from 6:01pm until 10pm and "Good night"
 * for want of anything more appropriate otherwise. The string
 * is evaluated at the time of the call and is overwritten between calls.
 */
char *greeting(void)
{
    static char Tstr[TSTAMP_LEN];
    
    time_t Now;
    struct tm *Tm;
    
    Now=time(NULL);
    Tm=localtime(&Now);
    if (6<=Tm->tm_hour && (Tm->tm_hour<=11 || 12==Tm->tm_hour && 0==Tm->tm_min))
	strcpy(Tstr, "Good morning");
    else if (12<=Tm->tm_hour && (Tm->tm_hour<=17 || 18==Tm->tm_hour && 0==Tm->tm_min))
	strcpy(Tstr, "Good afternoon");
    else if (18<=Tm->tm_hour && (Tm->tm_hour<=21 || 22==Tm->tm_hour && 0==Tm->tm_min))
	strcpy(Tstr, "Good evening");
    else strcpy(Tstr, "Good night");
    return(Tstr);
}
/* END of greeting() */

/* start_timer(): starts a timer by saving the current time in an
 * internal variable.
 */
void start_timer(void) { times(&Start); }

/* stop_timer(): stops the timer by saving the current time in
 * another internal variable. The start time is not affected.
 */
void stop_timer(void) { times(&Stop); }

/* timer_results(): returns the time (in integer seconds)
 * elapsed between the last calls to start_timer() and stop_timer().
 * User and system times for the current process and its children
 * are tracked separately. E.g. timer_results(TS_UTIME) returns the
 * user time of the process, timer_results(TS_CUTIME|TS_CSTIME)
 * returns the total time of the children. The four times can be
 * summed in any combination by OR-ing the corresponding constants
 * (see TsSelector_) together.
 */
long timer_results(int Sel)
{
    long T=0;
    
    if (Sel & TS_UTIME)	    /* user time */
    T+=Stop.tms_utime-Start.tms_utime;
    if (Sel & TS_STIME)	    /* system time */
    T+=Stop.tms_stime-Start.tms_stime;
    if (Sel & TS_CUTIME)    /* children user time */
    T+=Stop.tms_cutime-Start.tms_cutime;
    if (Sel & TS_CSTIME)    /* children system time */
    T+=Stop.tms_cstime-Start.tms_cstime;
    return(T/CLK_TCK);
}
/* END of timer_results() */

/* time_string(): when given a time (interval) in seconds in T,
 * this routine returns a ptr to an internal character buffer holding a
 * formatted string like "26 days 1 hour 3 mins 55 secs".
 * This string is overwritten by each call.
 * Fractions of seconds are not supported and T<0 is interpreted as T==0.
 */
const char *time_string(long T)
{
    static char Tstr[TIMESTR_LEN];
    char *Tcur;
    long D, H, M, S;
    ldiv_t Div;
    *Tstr='\0';
    
    if (T<=0)
    {
	sprintf(Tstr, "0 seconds");
	return(Tstr);
    }
    
    /* convert to days, hours, minutes */
    Div=ldiv(T, 86400L);
    D=Div.quot; T=Div.rem;
    Div=ldiv(T, 3600L);
    H=Div.quot; T=Div.rem;
    Div=ldiv(T, 60L);
    M=Div.quot; S=Div.rem;
    
    /* format to string */
    Tcur=Tstr;
    if (D)
	sprintf(Tcur, "%ld day%s ", D, ((D==1)? "": "s"));
    Tcur=Tstr+strlen(Tstr);
    if (H)
	sprintf(Tcur, "%ld hour%s ", H, ((H==1)? "": "s"));
    Tcur=Tstr+strlen(Tstr);
    if (M)
	sprintf(Tcur, "%ld min%s ", M, ((M==1)? "": "s"));
    Tcur=Tstr+strlen(Tstr);
    if (S)
	sprintf(Tcur, "%ld sec%s", S, ((S==1)? "": "s"));
    return(Tstr);    
}
/* END of time_string() */

/* ==== END OF FUNCTIONS tstamp.c ==== */
