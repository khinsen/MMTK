/*
 *  test program for lennaard jones potential
 *
 *  w. t. rankin
 *
 *  9/1/95
 *
 */

static char RCSid[] = "$Id: mpe_testLJ.c,v 1.2 1996/09/24 19:05:57 wrankin Exp $";

#include <stdio.h>
#include <math.h>
#include "mpe.h"


main(int argc, char **argv)
{

  int     p;
  Vector  f;
  Vector  v;
  Real    pot;
  Real    b_lj;
  MtypeLJ M1, M2, L1, L2, T1;


  p = 8;
  b_lj = 1.0;


  /* initialize arrays and alloc multipoles */
  LJinit(p);
  LJalloc(&M1,p);
  LJalloc(&M2,p);
  LJalloc(&L1,p);
  LJalloc(&L2,p);
  LJalloc(&T1,p);


  /* position of part1 wrt the multipole at <3,2,1> */
  /* add particle to multipole */

  v.x = -0.5; v.y = 0.0; v.z = 0.0;
  AddMultipoleLJ(M1,p,b_lj,v);

  /* shift the MPE to a new location */

  v.x = -0.5; v.y = 0.5; v.z = 0.0;
  M2M_LJ(M1,M2,p,v);

  /* shift multipole to a local expansion at <0,0.5,0> */

  v.x = -2.0; v.y = -2.0; v.z = -1.0;
  copyYI(T1,p,v);
  M2L_LJshort(M2, L2, T1, p);
  /* M2L_LJ(M2,L2,p,v); */

  /* shift local expansion  to <-0.5,0,0> */

  v.x = -0.5; v.y = -0.5; v.z = 0.0;
  L2L_LJ(L2,L1,p,v);

  /* evaluate the multipole at the position <0,0,0> */

  v.x = 0.5; v.y = 0.0; v.z = 0.0;
  Force_LJ(L1,p,b_lj,v,&pot,&f);

  printf("LJ Potential = %lg\n",pot);
  printf("LJ Force = [%lg,%lg,%lg]\n",f.x,f.y,f.z);

}

