/****************************************************************
*
*  dpmta_pvm.h - include file for parallel FMA PVM declarations
*
*  Copyright (c) 1994 Duke University
*  All rights reserved
*
*/

/* $Id: dpmta_pvm.h,v 2.11 1997/05/08 20:12:23 chumphre Exp $
 *
 * revision history:
 *
 * $Log: dpmta_pvm.h,v $
 * Revision 2.11  1997/05/08 20:12:23  chumphre
 * Fixed DataInPlace for SGIMP
 *
 * Revision 2.10  1997/02/24  19:12:23  wrankin
 * fixes to timing code
 *
 * Revision 2.9  1997/01/13  22:05:15  wrankin
 * added code to perform parallel macroscopic expansion computation
 *
 * Revision 2.8  1996/11/14  17:50:18  wrankin
 * performance enhancements to multipole M2L routine.
 * additions to make slaves exit gracefully.
 *
 * Revision 2.7  1996/10/18  17:04:47  wrankin
 * improved PVM message passing structures.
 * made current and cleaned up T3D library code
 * added additional test codes
 * implement new heirarchical make structure
 *
 * Revision 2.6  1995/01/01  21:40:43  wrankin
 * removed T3DBUG processing since it never really worked
 *
 * Revision 2.5  1994/11/30  18:08:43  wrankin
 * changed PVM message IDs to the range 0x1000-0x10ff
 *
 * Revision 2.4  1994/10/14  04:46:26  wrankin
 * added duke univ copyright notice
 *
 * Revision 2.3  1994/08/25  06:12:25  wrankin
 * added message type definitions for T3D message forwarding patch
 *
 * Revision 2.2  1994/07/25  23:10:42  wrankin
 * increased maximum number of processors to 256
 *
 * Revision 2.1  1994/06/02  12:52:27  wrankin
 * No change.  Update version number for Release 2
 *
 * Revision 1.6  1994/05/27  07:14:29  wrankin
 * added macro def for pvm_pkint() and pvm_upkint() for cray-t3d
 *
 * Revision 1.5  1994/04/09  03:46:08  wrankin
 * added dan gray's performance timing code
 *
 * Revision 1.4  1994/03/26  16:20:52  wrankin
 * added multipole expansion msg types for upward and downward passes
 *
 * Revision 1.3  1994/02/01  19:54:38  wrankin
 * increased MAXPROC from 16 to 64
 *
 * Revision 1.2  1994/01/31  10:19:24  wrankin
 * filled in message passing sections with dummy calls
 * code is still untested
 *
 * Revision 1.1  1994/01/30  17:31:01  wrankin
 * Initial revision
 *
*/


/* defines *****************/

/* maximum number of processors */
#define MAXPROC 256

/* define maximum message buffer size */
#define MAXBUFSIZE 32768

/*
*  define message IDs 
*
*  all messages will be in the range 0x1000-0x10ff and
*  the calling application should not use these message id's
*
*/

#define MSG_INIT1 0x1001
#define MSG_INIT2 0x1002
#define MSG_INTR1 0x1011
#define MSG_INTR2 0x1012
#define MSG_INTR3 0x1012
#define MSG_INTR4 0x1012
#define MSG_PART1 0x1021
#define MSG_PDIST 0x1031
#define MSG_MDIST 0x1032
#define MSG_MACRO 0x1033
#define MSG_M2P   0x1041
#define MSG_M2C   0x1042
#define MSG_L2P   0x1043
#define MSG_L2C   0x1044
#define MSG_RESLT 0x1051
#define MSG_TIME1 0x1061
#define MSG_TIME2 0x1062
#define MSG_EXIT  0x1071

/*
 *  patch to get around broken CRAY pvm_pkint routine
 *  note that stride must equal one
 */

#ifdef CRAY

#define pvm_pkint(iblk,size,stride) \
        pvm_pkbyte((char *)(iblk), (size) * sizeof(int), (stride))

#define pvm_upkint(iblk,size,stride) \
        pvm_upkbyte((char *)(iblk), (size) * sizeof(int), (stride))

#endif

/*
 * we need to define our data encoding strategies for message
 * passing.  this is because of the limitations in the use of
 * DataInPlace encoding for shared memory pvm versions, as well
 * as differences if we know we are going to run on a homogeneous
 * platform.
 *
 * the issue of using DataInPlace on a T3E, (where all the buffer
 * data must remain set until the message is received, as opposed
 * to just sent) will have to be decided locally.
 *
 */

#ifdef PVMHOMO 
#define DATA_NORMAL_PVM  PvmDataRaw

#if defined PVMSHMEM || defined SGIMP64
#define DATA_INPLACE_PVM PvmDataRaw
#else
#define DATA_INPLACE_PVM PvmDataInPlace
#endif

#else
#define DATA_NORMAL_PVM  PvmDataDefault
#define DATA_INPLACE_PVM PvmDataDefault
#endif

