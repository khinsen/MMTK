/*
 * PMTAlegendre.h
 *
 * Copyright (c) 1995 Duke University
 * All Rights Reserved.
 *
 * based on: Version 3.0, February 20, 1994 by Jim Leathrum
 *
 */

/*
 * RCS Id:
 *
 * $Id: mpe_legendre.h,v 1.1.1.1 1995/07/10 13:11:46 wrankin Exp $
 *
 * RSC History:
 *
 * $Log: mpe_legendre.h,v $
 * Revision 1.1.1.1  1995/07/10 13:11:46  wrankin
 * Initial release of the Multipole Library
 * Based upon W. Elliott's MDMA codes
 * Implements Colomb Potentials only (no LJ-potentials yet)
 *
 *
 */

/*
 * Legendre() sets up the Legendre Polynomial
 * adapted from Numerical Recipes in C
 *
 * P is type Mtype, xval is type double, and p is type unsigned int
 * p is the number of terms in the multipole expansion
 */

#define Legendre(P,p,xval)                                              \
{                                                                       \
    int Li, Lj;                                                         \
    double negterm, oddfact, nextodd, sqroot, sqrtterm;                 \
                                                                        \
    negterm = 1.0;                                                      \
    oddfact = 1.0;                                                      \
    nextodd = 1.0;                                                      \
    sqroot = sqrt(1.0 - xval*xval);                                     \
    sqrtterm = 1.0;                                                     \
    for(Li=0;Li < p;Li++){                                              \
        P[Li][Li] = negterm*oddfact*sqrtterm;                           \
        negterm *= -1.0;                                                \
        oddfact *= nextodd;                                             \
        nextodd += 2.0;                                                 \
        sqrtterm *= sqroot;                                             \
        if(Li < p-1){                                                     \
        P[Li+1][Li] = xval * (double)(2*Li+1) * P[Li][Li];              \
            for(Lj=Li+2;Lj < p;Lj++){                                  \
                P[Lj][Li] = (xval*(double)(2*Lj-1)*P[Lj-1][Li] -        \
                    (double)(Lj+Li-1)*P[Lj-2][Li])/(double)(Lj-Li);     \
            }                                                           \
        }                                                               \
    }                                                                   \
}

