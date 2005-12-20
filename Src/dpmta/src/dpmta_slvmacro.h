/* 
*  dpmta_macro.h - prototypes for DPMTA internal functions
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
 * $Id: dpmta_slvmacro.h,v 2.1 1997/11/07 16:49:31 wrankin Exp $
 *
 * RCS History:
 *
 * $Log: dpmta_slvmacro.h,v $
 * Revision 2.1  1997/11/07 16:49:31  wrankin
 * massive cleanup of code.
 *  - ansi-fication and inclusion of prototypes
 *  - removed unused variables
 *  - all (except the test) code compiles with minimal warnings under gcc.
 *
 *
 */

void MacroCleanup();
void MacroCompute( Mtype, Mtype );
void MacroInit( int, int, int, Real, int, int, int *, int );
void MacroPreComp( Vector, Vector, Vector, Real );

