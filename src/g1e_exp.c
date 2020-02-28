//
//  g1e_exp.c
//  libcint
//
//  Created by Tri Le on 2/24/20.
//  Copyright © 2020 Tri Le. All rights reserved.
//

#include "g1e_exp.h"

void CINTg_nuc_exp(double *g, double* rij, double* cr, double aij, double ak, double fac, CINTEnvVars *envs){
        const FINT nmax = envs->li_ceil + envs->lj_ceil;
        const FINT lj = envs->lj_ceil;
        const FINT dj = envs->g_stride_j;
        const double *ri = envs->ri;
        const double *rj = envs->rj;
        FINT i, j, ptr;
        double ririj[3], rirj[3];
        double *gx = g;
        double *gy = g + envs->g_size;
        double *gz = g + envs->g_size * 2;
        double pij[3];
        double aijk = aij + ak;
        double rijrk = 0, dr;
        for (FINT k = 0; k < 3; ++k){
            rirj[k] = ri[k] - rj[k];
            pij[k] = ( aij*rij[k] + ak*cr[k] ) / aijk;
            ririj[k] = ri[0] - pij[k];
            dr = rij[k] - cr[k];
            rijrk += dr*dr;
        }
        double fac2 = exp(-(aij*ak/ aijk)*rijrk)*fac;
        gx[0] = 1;
        gy[0] = 1;
        gz[0] = SQRTPI * PI * fac2;
        if (nmax > 0) {
                gx[1] = -ririj[0] * gx[0];
                gy[1] = -ririj[1] * gy[0];
                gz[1] = -ririj[2] * gz[0];
        }

        for (i = 1; i < nmax; i++) {
                gx[i+1] = 0.5 * i / aij * gx[i-1] - ririj[0] * gx[i];
                gy[i+1] = 0.5 * i / aij * gy[i-1] - ririj[1] * gy[i];
                gz[i+1] = 0.5 * i / aij * gz[i-1] - ririj[2] * gz[i];
        }
 
        for (j = 1; j <= lj; j++) {
                ptr = dj * j;
                for (i = ptr; i <= ptr + nmax - j; i++) {
                        gx[i] = gx[i+1-dj] + rirj[0] * gx[i-dj];
                        gy[i] = gy[i+1-dj] + rirj[1] * gy[i-dj];
                        gz[i] = gz[i+1-dj] + rirj[2] * gz[i-dj];
                }
        }
}
void CINTg_nuc_exp_ksum(double *g, double* rij, double* cr, double aij, double ak, double fac, CINTEnvVars *envs){
    const FINT nmax = envs->li_ceil + envs->lj_ceil;
    const FINT lj = envs->lj_ceil;
    const FINT dj = envs->g_stride_j;
    const double *ri = envs->ri;
    const double *rj = envs->rj;
    FINT i, j, ptr;
    double rirk[3], rirj[3];
    double *gx = g;
    double *gy = g + envs->g_size;
    double *gz = g + envs->g_size * 2;
    double pijk[3];
    double aijk = aij + ak;
    double rijrk = 0, dr;
    
    FINT nk = num_kterms(ak, fac, 1E-9);
    
    for (FINT k = 0; k < 3; ++k){
        pijk[k] = cr[k] - rij[k];
        rirk[k] = pijk[k] - ri[k];
        rirj[k] = ri[k] - rj[k];
        dr = rij[k] - cr[k];
        rijrk += dr*dr;
    }
    double complex gk[3*envs->g_size];
    double complex* gk_x = gk;
    double complex* gk_y = gk + envs->g_size;
    double complex* gk_z = gk + envs->g_size * 2;
    double expfac = PI*PI*ak*aij / (ak + aij);
    double complex ikpi[3];
    for (i = 0; i < 3*envs->g_size; ++i) g[i] = 0;
    
    double fac2 = exp(-(aij*ak/aijk)*rijrk)*fac;
    
    for (FINT k = 0; k < nk; ++k){
        
        gk_x[0] = exp(-k*k*expfac)*exp(2*I*PI*k*pijk[0]);
        gk_y[0] = exp(-k*k*expfac)*exp(2*I*PI*k*pijk[1]);
        gk_z[0] = exp(-k*k*expfac)*exp(2*I*PI*k*pijk[2]);
        
        if (nmax > 0) {
            
            ikpi[0] = -I*PI*k/aij + rirk[0];
            ikpi[1] = -I*PI*k/aij + rirk[1];
            ikpi[2] = -I*PI*k/aij + rirk[2];
            
            gk_x[1] = ikpi[0] * gk_x[0];
            gk_y[1] = ikpi[1] * gk_y[0];
            gk_z[1] = ikpi[2] * gk_z[0];
            
        }
        for (i = 1; i < nmax; i++) {
            gk_x[i+1] = 0.5 * i / aij * gk_x[i-1] + ikpi[0] * gk_x[i];
            gk_y[i+1] = 0.5 * i / aij * gk_y[i-1] + ikpi[1] * gk_y[i];
            gk_z[i+1] = 0.5 * i / aij * gk_z[i-1] + ikpi[2] * gk_z[i];
        }
        for (j = 1; j <= lj; j++) {
            ptr = dj * j;
            for (i = ptr; i <= ptr + nmax - j; i++) {
                
                gk_x[i] = gk_x[i+1-dj] + rirj[0] * gk_x[i-dj];
                gk_y[i] = gk_y[i+1-dj] + rirj[1] * gk_y[i-dj];
                gk_z[i] = gk_z[i+1-dj] + rirj[2] * gk_z[i-dj];
                
                gx[i] += cabs(gk_x[i]);
                gy[i] += cabs(gk_y[i]);
                gz[i] += cabs(gk_z[i])*fac2;
            }
        }
        for (i = 0; i <= nmax; i++){
            gx[i] += cabs(gk_x[i]);
            gy[i] += cabs(gk_y[i]);
            gz[i] += cabs(gk_z[i])*fac2;
        }
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
double JTsum_x0(double* gk_real, double* gk_img, double* R, double* alpha, FINT Nk){
    double aijk_inv = ( 1.0/alpha[0] + 1.0/alpha[1] )*M_PI*M_PI;
    double P = R[1] - R[0];
    double cosT, sinT;
    FINT k = 0;
    double out = 1;
    double fac;
    for (k = 0; k < Nk; ++k){
        fac = 2*exp(-k*k*aijk_inv);
        cosT = cos(2*M_PI*k*P)*fac;
        sinT = sin(2*M_PI*k+P);
        gk_real[k] = cosT;
        gk_img[k] = sinT;
        out += sqrt(cosT*cosT + sinT*sinT)*fac;
    }
    return out;
}
double JTsum_x1(double* gk_real, double* gk_img, double* R, double* alpha, FINT Nk){
    
    
}
double JTsum_x2(double* gk_real, double* gk_img, double* R, double* alpha, FINT Nk){
    
}
double JT_1D_sum(double aij, double ak, double wk, double pij_k, double cr_k){
    double fac = 1;
    /*
    exp [ - k^2 pi^2 / ak ] * integral exp[-aij*(x-Pij)^2 * exp[2*Pi*I*(x-R)]
    x' = x - R -> x' + R = x'
     
     exp [ - k^2 pi^2 / ak ] * integral exp[-aij*(x' + R -Pij )^2 * exp[2*Pi*I*x*k]
     
     = exp[- k^2 pi^2 / ak] * ( exp [ - k^2 * pi^2 / aij ] * exp [2*i*p]
     
     Ck = exp[ - k^2 pi^2 ( 1.0/ak + 1.0/ aij ] * exp[2*i*p] pi^2 / [sqrt(ak) ]
     
     x^j, j = 0 : I00k = Ck / sqrt(aij)
     
     */
    
    double aijk_inv = ( 1.0/ak + 1.0/aij )*M_PI*M_PI;
    double P = cr_k - pij_k;
    double cosT, sinT;
    FINT k = 0;
    double out = 1;
    while (fac > 1E-9){
        ++k;
        fac = 2*exp(-k*k*aijk_inv);
        cosT = cos(2*M_PI*k*P);
        sinT = sin(2*M_PI*k+P);
        
        out += sqrt(cosT*cosT + sinT*sinT)*fac;
    }
    return out;
}
/*
    Exp [ - pi^2 * k^2 / aij ] * exp[2*i*p] * (a*(p-r1) - i*pi*w)
    
    I_10k = I_00k*(p-r1) - i*pi*k / aij
 
    I_20k = I
 
 
 
 
 
 
 
 
 
 
 
 
 */
