/* 
*  dpmta_mastiter.h - prototypes for DPMTA internal functions
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
 * $Id: dpmta_mastiter.h,v 2.1 1997/11/07 16:49:10 wrankin Exp $
 *
 * RCS History:
 *
 * $Log: dpmta_mastiter.h,v $
 * Revision 2.1  1997/11/07 16:49:10  wrankin
 * massive cleanup of code.
 *  - ansi-fication and inclusion of prototypes
 *  - removed unused variables
 *  - all (except the test) code compiles with minimal warnings under gcc.
 *
 *
 */


/*
 * dpmta_mastiter.c
 */

void Send_Slave_Info( int, int *, int, int, int, int, int, int,
   double, int, int, int *);
void Delete_Local_Buffers();
void Init_Local_Buffers();
void Recv_Slave_Results();
void Return_Virial( double *, PmtaVector *, double *, PmtaVector * );
void Send_Slave_Particles( int, PmtaParticlePtr, int, int *, int, int,
   PmtaVector *, PmtaVector *, PmtaVector *, PmtaVector *, int );


