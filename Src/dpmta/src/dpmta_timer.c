/*
*  dpmta_timer.c - procedures to collect and report timing information
*
*  w. t. rankin
*
*  Copyright (c) 1997 Duke University
*  All rights reserved
*
*/

static char rcsid[] = "$Id: dpmta_timer.c,v 2.7 1998/01/05 13:42:01 wrankin Exp $";

/*
 * Revision history:
 *
 * $Log: dpmta_timer.c,v $
 * Revision 2.7  1998/01/05 13:42:01  wrankin
 * added more timing code.
 *
 * Revision 2.6  1997/11/07 16:54:25  wrankin
 * massive cleanup of code.
 *  - ansi-fication and inclusion of prototypes
 *  - removed unused variables
 *  - all (except the test) code compiles with minimal warnings under gcc.
 *
 * Revision 2.5  1997/04/10 17:26:15  wrankin
 * improved particle sorting algorithm from serial runs.
 * removed extraneous timing code frm within sort proc.
 *
 * Revision 2.4  1997/03/26  20:29:56  wrankin
 * cleaned up timing code outputs
 *
 * Revision 2.3  1997/03/12  20:15:18  wrankin
 * updates to timer codes
 *
 * Revision 2.2  1997/03/04  19:18:18  wrankin
 * updates to timing codes
 *
 * Revision 2.1  1997/02/26  20:43:34  wrankin
 * updated timing measurements and placed routines in single module
 *
 *
*/

#ifdef TIME

/* include files */

#include <stdio.h>
#include "dpmta_timer.h"
#ifndef SERIAL
#include "pvm3.h"
#include "dpmta_pvm.h"
#endif


/*
 * internal prototypes
 */

void Print_Common( double * );

/*
*  some timing and profiling storage, if needed
*  external declarations in dpmta_timer.h
*/

struct tms startbuf, endbuf;
struct timeval runstruct;

long   time_tmpint = 0;          /* temp timing value */
double times_arr[TIMINGSIZE];    /* array of timing values for master */


#ifdef SERIAL

/****************************************************************
 *
 * Print_Times() -  printout slave timing information
 *   for serial DPMTA.
 *
 */

void Print_Times( double *timings )
{

   fprintf(stderr, "Timing Information:\n");
   Print_Common(timings);
   fflush(stderr);
}

#else

/****************************************************************
*
*  Recv_Slave_Times() - receive and printout slave timing
*     information.
*
*/

void Recv_Slave_Times( int mypid, int nprocs )

{
   int    i;                    /* loop counter */
   int    pid;                  /* processor id of slave */
   double timings[MAXPROC][16]; /* slave timing data */

   if ( mypid != 0 ) {
       return;
   }

   /* receive slave timing data */
   for (i=0; i<nprocs; i++) {
      pvm_recv(-1, MSG_TIME2);

      pvm_upkint(&pid,1,1);
      pvm_upkdouble(&(timings[pid][0]),TIMINGSIZE,1);
   }

   for (i=0; i<nprocs; i++) {
      fprintf(stderr, "PID %d\n", i);
      Print_Common(&(timings[i][0]));
   } /* for i */

   fflush(stderr);
}

/****************************************************************
 *
 *  Send_Slave_Times() - send timing results to master
 *  these are later picked up by the Recv_Slave_Times() 
 *
 */


void Send_Slave_Times( int pid, double *timings, int tid )
{
   pvm_initsend(DATA_NORMAL_PVM);
   pvm_pkint(&(pid),1,1);
   pvm_pkdouble(timings,TIMINGSIZE,1);
   pvm_send(tid, MSG_TIME2);
}


#endif /* not SERIAL */


/****************************************************************
 *
 *  Print_Common() - common output format for timing results
 *
 *  Here are the array fields and what they contain
 *
 *  0 - wall clock start time right after Recv_Particles()
 *  1 - elapse time to after Send_Results()
 *  2 - elapse time to after Slave_Send_Multipole()
 *  3 - elapse time to after Slave_Recv_Multipole()
 *  4 - elapse time to after Slave_MPE_Force()
 *
 *  7 - cpu time for particle sorting.
 *  8 - cpu time for entire resizing code.
 *  9 - cpu time for Slave_Mpole_Exp() - upward pass
 * 10 - cpu time for Slave_Direct_Calc() - direct comp
 * 11 - cpu time for Slave_MPE_Calc() - downward pass
 * 12 - cpu time for Slace_MPE_Force() - local expansion forces
 * 13 - cpu time for MacroPreComp() - macroscopic expansion
 * 14 - cpu time for actual macroscopic compuation.
 *
 */

void Print_Common(double *timings)
{

   fprintf(stderr, "CPU Times:\n");
   fprintf(stderr, "  ");
   fprintf(stderr, "upward: %g  ", timings[9]);
   fprintf(stderr, "direct: %g  ", timings[10]);
   fprintf(stderr, "downward: %g  ", timings[11]);
   fprintf(stderr, "mpe force: %g  \n", timings[12]);

   fprintf(stderr, "  ");
   fprintf(stderr, "sorting: %g  ", timings[7]);
   fprintf(stderr, "resize: %g  ", timings[8]);
   fprintf(stderr, "macro pre: %g  ", timings[13]);
   fprintf(stderr, "macro: %g \n", timings[14]);

   fprintf(stderr, "Elapse Times:\n");
   fprintf(stderr, "  ");
   fprintf(stderr, "Send_Init_Particles: %g  ", timings[5]);
   fprintf(stderr, "Recv_Init_Particles: %g\n", timings[6]);
   fprintf(stderr, "  ");
   fprintf(stderr, "Send_Multipoles: %g  ", timings[2]);
   fprintf(stderr, "Recv_Multipoles: %g \n", timings[3]);
   fprintf(stderr, "  ");
   fprintf(stderr, "MPE_Force: %g  ", timings[4]);
   fprintf(stderr, "Send_Results: %g \n", timings[1]);
   fprintf(stderr, "\n");
}

#endif /* TIME */
