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
void CINTg_nuc_exp_ksum(double *g, double* rij, double* cr, double aij, double ak, double fac, CINTEnvVars *envs);
FINT num_kterms(double ak, double wk, double epsilion = 1E-10);
double JTsum(double* gk_real, double* gk_img, double* Ri, double* Rj, double* Rc, double* alpha, FINT Nk, FINT D);
double JTsum_x0(double* gk_real, double* gk_img, double* R, double* alpha, FINT Nk);
double JTsum_x1(double* gk_real, double* gk_img, double* R, double* alpha, FINT Nk);
double JTsum_x2(double* gk_real, double* gk_img, double* R, double* alpha, FINT Nk);
#endif /* g1e_exp_h */
