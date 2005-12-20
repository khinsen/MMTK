/****************************************************************
*
*  dpmta_timer.h - global variables for dpmta timer 
*    routines
*
*  w. t. rankin
* 
*  these are the external references to the  global variables
*  used by the slave process for timing
*
*  Copyright (c) 1994 Duke University
*  All rights reserved
*
*  RCS info: $Id: dpmta_timer.h,v 2.8 1997/11/07 16:54:26 wrankin Exp $
*/

/*
 *  revision history:
 *
 *  $Log: dpmta_timer.h,v $
 *  Revision 2.8  1997/11/07 16:54:26  wrankin
 *  massive cleanup of code.
 *   - ansi-fication and inclusion of prototypes
 *   - removed unused variables
 *   - all (except the test) code compiles with minimal warnings under gcc.
 *
 *  Revision 2.7  1997/03/12 20:12:56  wrankin
 *  added support for SGIMP64 architecture (IRIX 6.2)
 *
 * Revision 2.6  1997/03/04  19:18:21  wrankin
 * updates to timing codes
 *
 * Revision 2.5  1997/02/26  20:43:35  wrankin
 * updated timing measurements and placed routines in single module
 *
 * Revision 2.4  1997/02/24  19:12:31  wrankin
 * fixes to timing code
 *
 * Revision 2.3  1997/01/28  17:10:47  wrankin
 * added new timing codes for macroscopic expansion calculations
 *
 * Revision 2.2  1996/10/18  17:05:13  wrankin
 * improved PVM message passing structures.
 * made current and cleaned up T3D library code
 * added additional test codes
 * implement new heirarchical make structure
 *
 * Revision 2.1  1995/06/13  04:26:26  wrankin
 * No updates.
 * Sources placed under CVS.
 *
 * Revision 1.1  1995/04/24  04:19:52  wrankin
 * Initial revision
 *
 *
*/

#include <sys/times.h>

#if defined HPPA
#   include<time.h>
#elif defined LINUX
#   include <sys/time.h>
#elif defined CRAY
#   include <sys/time.h>
#elif defined SGI64
#   include <sys/time.h>
#elif defined SGIMP64
#   include <sys/time.h>
#elif defined RS6K
#   include <sys/time.h>
#   define CLK_TCK 100
#endif

#ifndef CLK_TCK
#   define CLK_TCK 60
#endif

/* size of timing array */

#define TIMINGSIZE 16

extern struct tms startbuf, endbuf;
extern struct timeval runstruct;

extern long time_tmpint;     /* temp timing value */

extern double times_arr[];    /* array of timing values for master */


/*
 * dpmta_timer.c
 */

void Recv_Slave_Times( int, int );
void Send_Slave_Times( int, double *, int );
void Print_Times( double * );

