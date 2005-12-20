/* 
*  dpmta_slvmkhl.h - prototypes for DPMTA internal functions
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
 * $Id: dpmta_slvmkhl.h,v 2.1 1997/11/07 16:49:41 wrankin Exp $
 *
 * RCS History:
 *
 * $Log: dpmta_slvmkhl.h,v $
 * Revision 2.1  1997/11/07 16:49:41  wrankin
 * massive cleanup of code.
 *  - ansi-fication and inclusion of prototypes
 *  - removed unused variables
 *  - all (except the test) code compiles with minimal warnings under gcc.
 *
 *
 */

void Compute_Hlist( int );
void Init_Hlist();
void Make_Hlist();
void Make_RelVec( int );
void Delete_Hlist();
