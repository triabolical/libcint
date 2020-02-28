//
//  g1e_lim.h
//  
//
//  Created by Tri Le on 2/14/20.
//

#ifdef I8
#include <stdint.h>
#define FINT int64_t
#else
#define FINT int
#endif

#ifndef g1e_lim_h
#define g1e_lim_h

#include <stdio.h>
#include "g1e.h"
#include "cint_const.h"
#include <math.h>
#include "cint_bas.h"
#include "misc.h"

void SCsum3D(double * abc, double * Pij, double * cr, double aijk, double* v, double fac, int* nt);
void SCsum2D(double * abc, double * Pij, double * cr, double aijk, double* v, double fac, int* nt);

void CINTg_nuc_all(double *g, double aij, double *rij, double *cr, double ak, double fac, CINTEnvVars *envs);
void CINTg_nuc_lim(double *g, double aij, double *rij, double *cr, double ak, double fac, CINTEnvVars *envs);

void CINTg_nuc_far(double *g, double aij, double *rij, double *cr, double ck, double fac, CINTEnvVars *envs);


#endif /* g1e_lim_h */
