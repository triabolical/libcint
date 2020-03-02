//
//  g1e_exp.c
//  libcint
//
//  Created by Tri Le on 2/24/20.
//  Copyright © 2020 Tri Le. All rights reserved.
//

#include "g1e_exp.h"

void CINTg_nuc_exp(double *g, double* rij, double* cr, double aij, double ak, double fac, CINTEnvVars *envs){
    
    /*
     rij = ( ai*ri + aj*rj ) / aij
     cr = Rc nuc position
     fac = exp(-ai*aj/aij*Rij^2) * pk
     */
    
    FINT gsize = envs->g_size;
    double gnear[3*gsize];
    double gper[3*gsize];
    double g0[3*gsize];
    double* gN[3], *gP[3],*gK0[3];
    
    for (FINT k = 0; k < 3; ++k){
        gN[k] = gnear + k*gsize;
        gP[k] = gper + k*gsize;
        gK0[k] = g0 + k*gsize;
    }
    
    const FINT nmax = envs->li_ceil + envs->lj_ceil;
    const FINT lj = envs->lj_ceil;
    const FINT dj = envs->g_stride_j;
    const double *ri = envs->ri;
    const double *rj = envs->rj;
    FINT i, j, ptr, k;
    double rirk[3], rirj[3], ririj[3], rijk[3], ririjk[3];
    double *gx = g;
    double *gy = g + envs->g_size;
    double *gz = g + envs->g_size * 2;
    double pijk[3];
    double aijk = aij + ak;
    double rijrk = 0, dr;
    
    for (FINT k = 0; k < 3; ++k){
        rijk[k] = ( rij[k]*aij + cr[k]*ak ) / aijk;
        pijk[k] = cr[k] - rij[k];
        rirk[k] = pijk[k] - ri[k];
        rirj[k] = ri[k] - rj[k];
        ririj[k] = ri[k] - rij[k];
        dr = rij[k] - cr[k];
        rijrk += dr*dr;
    }
    
    // gnear and gk=0
    
    double Cnear = exp(-aij*ak*rijrk/aijk) / (aijk * sqrt(aijk) ) ;
    double Ck = PI * PI * PI / (ak * sqrt(ak) * sqrt(aij) * aij );
    
    for (k = 0; k < 3; ++k){
        ririjk[k] = ri[k] - rijk[k];
        ririj[k] = ri[k] - rij[k];
    }
    
    gN[0][0] = 1;
    gN[1][0] = 1;
    gN[2][0] = SQRTPI * PI * fac * Cnear;
    
    gK0[0][0] = 1;
    gK0[1][0] = 1;
    gK0[2][0] = fac * Ck;
    
    if (nmax > 0) {
        for (k = 0; k < 3; ++k){
            gN[k][1] = -ririjk[k] * gN[k][0];
            gK0[k][1] = -ririj[k] * gK0[k][0];
        }
    }
    
    for (i = 1; i < nmax; i++) {
        for (k = 0; k < 3; ++k){
            gN[k][i+1] = 0.5 * i / aij * gN[k][i-1] - ririjk[0] * gN[k][i];
            gK0[k][i+1] = 0.5 * i / aij * gK0[k][i-1] - ririj[0] * gK0[k][i];
        }
    }
    
    for (j = 1; j <= lj; j++) {
        ptr = dj * j;
        for (i = ptr; i <= ptr + nmax - j; i++) {
            for (k = 0; k < 3; ++k){
                gN[k][i] = gN[k][i+1-dj] + rirj[k] * gN[k][i-dj];
                gK0[k][i] = gK0[k][i+1-dj] + rirj[k] * gK0[k][i-dj];
            }
        }
    }
    
    // k sum
    FINT nk = num_kterms(ak, fac, 1E-9);
    double complex gk[3*envs->g_size];
    double complex* gkx = gk;
    double complex* gky = gk + envs->g_size;
    double complex* gkz = gk + envs->g_size * 2;
    double expfac = PI*PI*aijk/(ak*aij);
    double complex ikpi[3];
    for (i = 0; i < 3*envs->g_size; ++i) g[i] = 0;
    
    for (FINT k = 0; k < nk; ++k){
        
        gkx[0] = exp(-k*k*expfac)*exp(2*I*PI*k*pijk[0]);
        gky[0] = exp(-k*k*expfac)*exp(2*I*PI*k*pijk[1]);
        gkz[0] = exp(-k*k*expfac)*exp(2*I*PI*k*pijk[2]);
        
        if (nmax > 0) {
            
            ikpi[0] = -I*PI*k/aij + rirk[0];
            ikpi[1] = -I*PI*k/aij + rirk[1];
            ikpi[2] = -I*PI*k/aij + rirk[2];
            
            gkx[1] = ikpi[0] * gkx[0];
            gky[1] = ikpi[1] * gky[0];
            gkz[1] = ikpi[2] * gkz[0];
            
        }
        
        for (i = 1; i < nmax; i++) {
            gkx[i+1] = 0.5 * i / aij * gkx[i-1] + ikpi[0] * gkx[i];
            gky[i+1] = 0.5 * i / aij * gky[i-1] + ikpi[1] * gky[i];
            gkz[i+1] = 0.5 * i / aij * gkz[i-1] + ikpi[2] * gkz[i];
        }
        for (j = 1; j <= lj; j++) {
            ptr = dj * j;
            for (i = ptr; i <= ptr + nmax - j; i++) {
                
                gkx[i] = gkx[i+1-dj] + rirj[0] * gkx[i-dj];
                gky[i] = gky[i+1-dj] + rirj[1] * gky[i-dj];
                gkz[i] = gkz[i+1-dj] + rirj[2] * gkz[i-dj];
                
                gx[i] += cabs(gkx[i]);
                gy[i] += cabs(gky[i]);
                gz[i] += cabs(gkz[i])*fac * Ck;
            }
        }
        for (i = 0; i <= nmax; i++){
            gx[i] += cabs(gkx[i]);
            gy[i] += cabs(gky[i]);
            gz[i] += cabs(gkz[i])*fac * Ck;
        }
    }
    
    double *gxyz[3] = {gx, gy, gz};
    for (i = 0; i < gsize; ++i){
        for (k = 0; k < 3; ++k)
        gxyz[k][i] = gN[k][i] - gK0[k][i] + gP[k][i];
    }
}
FINT num_kterms(double ak, double wk, double epsilion){
    FINT N = 1;
    double exp_fac = -PI*PI / ak;
    double fac = wk;
    while (fac > epsilion){
        ++N;
        fac = wk*exp(N*N*exp_fac);
    }
    return N;
}

