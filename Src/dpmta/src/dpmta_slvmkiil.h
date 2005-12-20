/* 
*  dpmta_slvmkiil.h - prototypes for DPMTA internal functions
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
 * $Id: dpmta_slvmkiil.h,v 2.1 1997/11/07 16:49:44 wrankin Exp $
 *
 * RCS History:
 *
 * $Log: dpmta_slvmkiil.h,v $
 * Revision 2.1  1997/11/07 16:49:44  wrankin
 * massive cleanup of code.
 *  - ansi-fication and inclusion of prototypes
 *  - removed unused variables
 *  - all (except the test) code compiles with minimal warnings under gcc.
 *
 *
 */


/*
 * dpmta_slvmkiil.c
 */

void Init_Inv_Ilist();
void Make_Inv_Ilist();
void Delete_Inv_Ilist();


