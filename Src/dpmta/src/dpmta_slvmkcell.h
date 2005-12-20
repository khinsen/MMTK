/* 
*  dpmta_slvmkcell.h - prototypes for DPMTA internal functions
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
 * $Id: dpmta_slvmkcell.h,v 2.2 1997/11/12 13:49:56 wrankin Exp $
 *
 * RCS History:
 *
 * $Log: dpmta_slvmkcell.h,v $
 * Revision 2.2  1997/11/12 13:49:56  wrankin
 * updates to communications routines and general cleanup of code in prep
 *   for introducing load balancing functionality in hte near future
 *
 * Revision 2.1  1997/11/07 16:49:38  wrankin
 * massive cleanup of code.
 *  - ansi-fication and inclusion of prototypes
 *  - removed unused variables
 *  - all (except the test) code compiles with minimal warnings under gcc.
 *
 *
 */


void Alloc_Cell_Table();
void Alloc_Ilist_Cells();
void Make_Cell_Centers();
void Delete_Cell_Table();
void Realloc_Cell_Table();
