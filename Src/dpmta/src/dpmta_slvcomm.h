/* 
*  dpmta_slvcomm.h - prototypes for DPMTA internal functions
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
 * $Id: dpmta_slvcomm.h,v 2.3 1998/04/01 20:08:15 wrankin Exp $
 *
 * RCS History:
 *
 * $Log: dpmta_slvcomm.h,v $
 * Revision 2.3  1998/04/01 20:08:15  wrankin
 * added support for HILBERT ordering
 * general cleanup of code structure (more OOP if you will)
 *
 * Revision 2.2  1997/11/12 13:49:51  wrankin
 * updates to communications routines and general cleanup of code in prep
 *   for introducing load balancing functionality in hte near future
 *
 * Revision 2.1  1997/11/07 16:49:22  wrankin
 * massive cleanup of code.
 *  - ansi-fication and inclusion of prototypes
 *  - removed unused variables
 *  - all (except the test) code compiles with minimal warnings under gcc.
 *
 *
 */

void Recv_Master_Info();

void Comm_Init();
void Comm_Delete();

int Recv_Particles();
void Send_Results();

void Recv_Mpe_from_Child( int );
void Recv_Lcl_from_Parent( int );
void Send_Mpe_to_Parent( int );
void Send_Lcl_to_Child( int );

void Slave_Send_SDirect();
void Slave_Recv_SDirect();
void Slave_Send_Multipole();
void Slave_Recv_Multipole();

