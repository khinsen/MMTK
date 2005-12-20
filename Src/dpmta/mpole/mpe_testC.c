/*
 *  test program for Macroscopic potential
 *
 *  w. t. rankin
 *
 *  9/1/95
 *
 */

static char RCSid[] = "$Id: mpe_testC.c,v 1.1 1996/09/24 19:05:52 wrankin Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpe.h"


main(int argc, char **argv)
{

  int     i,j,k,l,m,n;
  int     x,y,z;
  int     p;
  Real    q;
  Real    pot;
  Vector  v, vp[20];
  Vector  f;
  Mtype   M1, M2, M3, M4;
  Mtype   T1, T2;
  Mtype   L1, L2, L3, L4;

  p = 8;
  q = 1.0;


  /* initialize arrays and alloc multipoles */
  Cinit(p);
  Calloc(&M1,p);
  Calloc(&M2,p);
  Calloc(&M3,p);
  Calloc(&M4,p);
  Calloc(&T1,p);
  Calloc(&T2,p);
  Calloc(&L1,p);
  Calloc(&L2,p);
  Calloc(&L3,p);
  Calloc(&L4,p);

  /* create a set of random particles */

  for ( i=0; i<10; i++) {
    vp[i].x = drand48() - 0.5;
    vp[i].y = drand48() - 0.5;
    vp[i].z = drand48() - 0.5;
  }
  for ( i=0; i<10; i++) {
    printf("part = (%lg, %lg, %lg)\n",vp[i].x,vp[i].y,vp[i].z);
  }


  /* add eight particles to a multipole */
  CMclear(M1,p);
  for ( i=0; i<10; i++ ) {
    AddMultipoleC(M1,p,q,vp[i]);
  }

  /* construct the M2M Xfer matrix */
  CMclear(T1,p);
  for ( i=-1; i<=1; i+=2 ) {
    for ( j=-1; j<=1; j+=2 ) {
      for ( k=-1; k<=1; k+=2 ) {
         v.x = (double)i;
         v.y = (double)j;
         v.z = (double)k;
	 addF(T1,p,v);
      }
    }
  }

  /* create the central multipole from 8 copies of
   * the original
   */

  CMclear(M2,p);
  M2M_Cshort(M1,M2,T1,p);

  /*
   * print out this multipole
   */
  printf("Macro-Multipole:\n==============\n");
  for ( n=0; n<p; n++ ) {
    for ( m=0; m<=n; m++ ) {
      printf("(%lg,%lg) ",M2[n][m].x,M2[n][m].y);
    }
    printf("\n");
  }
  printf("\n");


  /*
   * now construct the multipole
   * directly from eh particle data
   */

  CMclear(M3,p);
  for ( i=0; i<10; i++ ) {
    for ( x=-1; x<=1; x+=2 ) {
      for ( y=-1; y<=1; y+=2 ) {
	for ( z=-1; z<=1; z+=2 ) {
	  v.x = vp[i].x + (double)x;
	  v.y = vp[i].y + (double)y;
	  v.z = vp[i].z + (double)z;

	  AddMultipoleC(M3,p,q,v);
	}
      }
    }
  }

  /*
   * print out this multipole
   */
  printf("Directly Constructed-Multipole:\n==============\n");
  for ( n=0; n<p; n++ ) {
    for ( m=0; m<=n; m++ ) {
      printf("(%lg,%lg) ",M3[n][m].x,M3[n][m].y);
    }
    printf("\n");
  }
  printf("\n");


  /* construct the multipole directly from M2M calls */
  CMclear(M4,p);
  for ( i=-1; i<=1; i+=2 ) {
    for ( j=-1; j<=1; j+=2 ) {
      for ( k=-1; k<=1; k+=2 ) {
         v.x = (double)i;
         v.y = (double)j;
         v.z = (double)k;

	 M2M_C(M1,M4,p,v);
      }
    }
  }

  /*
   * print out this multipole
   */
  printf("M2M multipole-Multipole:\n==============\n");
  for ( n=0; n<p; n++ ) {
    for ( m=0; m<=n; m++ ) {
      printf("(%lg,%lg) ",M4[n][m].x,M4[n][m].y);
    }
    printf("\n");
  }
  printf("\n");


  /*
   * now check the potentials and forces at distant regions
   */

  v.x = 2.0;
  v.y = 3.7;
  v.z = 0.2;

  ForceM_C(M2,p,q,v,&pot,&f);
  printf("M2 Pot/Force = %lg (%lg,%lg,%lg)\n",pot,f.x,f.y,f.z);

  ForceM_C(M3,p,q,v,&pot,&f);
  printf("M3 Pot/Force = %lg (%lg,%lg,%lg)\n",pot,f.x,f.y,f.z);

  ForceM_C(M4,p,q,v,&pot,&f);
  printf("M4 Pot/Force = %lg (%lg,%lg,%lg)\n",pot,f.x,f.y,f.z);


  /*
   * now we will construct several local expansions.
   */

  /* first use some vanilla M2L calls */

  CMclear(L1,p);
  for ( i=0; i<3; i++ ) {
    v.x = (double)i;
    v.y = 3.0;
    v.z = 0.2;
    M2L_C(M2,L1,p,v);
  }

  printf("Normal Local:\n==============\n");
  for ( n=0; n<p; n++ ) {
    for ( m=0; m<=n; m++ ) {
      printf("(%lg,%lg) ",L1[n][m].x,L1[n][m].y);
    }
    printf("\n");
  }
  printf("\n");


  /* now use M2L Xfer matrices instead */
  CMclear(T2,p);
  CMclear(L2,p);
  for ( i=0; i<3; i++ ) {
    v.x = (double)i;
    v.y = 3.0;
    v.z = 0.2;
    addG(T2,p,v);
  }
  M2L_Cshort(M2,L2,T2,p);

  printf("Xfer'd Local:\n==============\n");
  for ( n=0; n<p; n++ ) {
    for ( m=0; m<=n; m++ ) {
      printf("(%lg,%lg) ",L2[n][m].x,L2[n][m].y);
    }
    printf("\n");
  }
  printf("\n");

}

