/* 
*  dpmta_distmisc.h - prototypes for DPMTA internal functions
*
*  w. t. rankin
*
*  Copyright (c) 1997 Duke University
*  All rights reserved
*
*  this files contains the prototype definitions for the external
*  functions provided by the corresponding '.c' file.
*
*/

/*
 * $Id: dpmta_distmisc.h,v 2.2 1998/04/01 20:08:06 wrankin Exp $
 *
 * RCS History:
 *
 * $Log: dpmta_distmisc.h,v $
 * Revision 2.2  1998/04/01 20:08:06  wrankin
 * added support for HILBERT ordering
 * general cleanup of code structure (more OOP if you will)
 *
 * Revision 2.1  1997/11/07 16:49:01  wrankin
 * massive cleanup of code.
 *  - ansi-fication and inclusion of prototypes
 *  - removed unused variables
 *  - all (except the test) code compiles with minimal warnings under gcc.
 *
 *
 */


/*
 * dpmta_distmisc.c
 */

void Dist_Init( int );
void Dist_Delete( int );

int geteindex( int, int, int );
int getfirstchild( int );
int getparent( int );
int getsindex( int, int, int );
int getslvpid( int, int, int );
int getslvpid_indx( int, int, int );
int index2cell( int, int );
int cell2index( int, int );


