/* 
*  dpmta_compute.h - prototypes for DPMTA internal functions
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
 * $Id: dpmta_slvcompute.h,v 2.2 1998/03/10 22:22:04 wrankin Exp $
 *
 * RCS History:
 *
 * $Log: dpmta_slvcompute.h,v $
 * Revision 2.2  1998/03/10 22:22:04  wrankin
 * folded start/cleanup functionality into dpmta_slvcompute
 *
 * Revision 2.1  1997/11/07 16:49:25  wrankin
 * massive cleanup of code.
 *  - ansi-fication and inclusion of prototypes
 *  - removed unused variables
 *  - all (except the test) code compiles with minimal warnings under gcc.
 *
 *
 */

void Slave_Init();
void Slave_Start();
void Slave_Compute();
void Slave_Cleanup();
void Slave_Delete();
