#ifndef TSTAMP_H
#define TSTAMP_H

/* ==== HEADER tstamp.h ==== */

/* Useful formatted time routines. */

/* ANSI C, IRIX 5.3, 30. July 1996. Andris Aszodi */

#ifdef __cplusplus
extern "C" {
#endif

/* ---- TYPEDEFS ---- */

/* TsSelector_: user and sys times of the current process and
 * its children may be queried after each start_timer() / stop_timer()
 * call pair. These enums decide which is returned by timer_results().
 * By OR-ing them together, sums of these times may be obtained.
 */
typedef enum { TS_UTIME=1, TS_STIME=2, TS_CUTIME=4, TS_CSTIME=8 } TsSelector_;

/* ---- PROTOTYPES ---- */

/* time_stamp(): returns a pointer to an internal static string
 * that contains a time string in the following format:-
 * "Thu 02-Jun-1994 18:24:23". The string is evaluated at the
 * time of the call and is overwritten between calls.
 */
char *time_stamp(void);

/* greeting(): returns a pointer to an internal static string
 * that contains a greeting appropriate for the time of the day:
 * "Good morning" from 6am till noon, "Good afternoon" from 12:01pm
 * until 6pm, "Good evening" from 6:01pm until 10pm and "Good night"
 * for want of anything more appropriate otherwise. The string
 * is evaluated at the time of the call and is overwritten between calls.
 */
char *greeting(void);

/* start_timer(): starts a timer by saving the current time in an
 * internal variable.
 */
void start_timer(void);

/* stop_timer(): stops the timer by saving the current time in
 * another internal variable. The start time is not affected.
 */
void stop_timer(void);

/* timer_results(): returns the time (in integer seconds)
 * elapsed between the last calls to start_timer() and stop_timer().
 * User and system times for the current process and its children
 * are tracked separately. E.g. timer_results(TS_UTIME) returns the
 * user time of the process, timer_results(TS_CUTIME|TS_CSTIME)
 * returns the total time of the children. The four times can be
 * summed in any combination by OR-ing the corresponding constants
 * (see TsSelector_) together.
 */
long timer_results(int Sel);

/* time_string(): when given a time (interval) in seconds in T,
 * this routine returns a ptr to an internal character buffer holding a
 * formatted string like "26 days 1 hour 3 mins 55 secs".
 * This string is overwritten by each call.
 * Fractions of seconds are not supported and T<0 is interpreted as T==0.
 */
const char *time_string(long T);

#ifdef __cplusplus
}
#endif

/* ==== END OF HEADER tstamp.h ==== */
#endif /* TSTAMP_H */
