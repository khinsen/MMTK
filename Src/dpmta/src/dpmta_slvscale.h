/* 
*  dpmta_slvscale.h - prototypes for DPMTA internal functions
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
 * $Id: dpmta_slvscale.h,v 2.1 1997/11/07 16:54:14 wrankin Exp $
 *
 * RCS History:
 *
 * $Log: dpmta_slvscale.h,v $
 * Revision 2.1  1997/11/07 16:54:14  wrankin
 * massive cleanup of code.
 *  - ansi-fication and inclusion of prototypes
 *  - removed unused variables
 *  - all (except the test) code compiles with minimal warnings under gcc.
 *
 *
 */

double Vec_Mag(Vector *);
double Max_CellLength();
void Rescale_Particles();
void Rescale_Results();
void MultipoleResize();

