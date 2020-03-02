//
//  g1e_exp.h
//  libcint
//
//  Created by Tri Le on 2/24/20.
//  Copyright © 2020 Tri Le. All rights reserved.
//

#ifdef I8
#include <stdint.h>
#define FINT int64_t
#else
#define FINT int
#endif

#ifndef g1e_exp_h
#define g1e_exp_h

#include <stdio.h>
#include "g1e.h"
#include "cint_const.h"
#include <math.h>
#include <complex.h>

void CINTg_nuc_exp(double *g, double* rij, double* cr, double aij, double ak, double fac, CINTEnvVars *envs);
FINT num_kterms(double ak, double wk, double epsilion);
#endif /* g1e_exp_h */
