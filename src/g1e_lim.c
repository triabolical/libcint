//
//  g1e_lim.c
//  
//
//  Created by Tri Le on 2/14/20.
//

#include "g1e_lim.h"
/*
void SCsum3D(double * abc, double * Pij, double * cr, double aijk, double* v, double fac, int* nt){
    int n[3] = {0, 0, 0}, k, q;
    double Rn[3];
    double eij, r2;
    nt = 0;
    int N = 0;
    for ( n[2] = -3 ; n[2] <= 3; ++n[2]){
        for (n[1] = -3; n[1] <= 3; ++n[1]){
            for (n[0] = -3; n[0] <= 3; ++n[0]){
            //    if (n[0] == 0 && n[1] == 0 && n[2] == 0) continue;
                for (k = 0; k < 3; ++k){
                    Rn[k] = 0;
                    for (q = 0; q < 3; ++q) Rn[k] += abc[3*k+q]*n[q];
                }
                eij = exp(-aijk* CINTsquare_dist(Pij, Rn) ) * fac;
                if (eij > EXPCUTOFF){
                    for (k = 0; k < 3; ++k) v[3*N+k] = Rn[k];
		    v[3*N+3] = eij;
                    ++N;
                }
            }
        }
    }
    *nt = N;
}
void SCsum2D(double * abc, double * Pij, double * cr, double aijk, double* v, double fac, int* nt){
    int n[3] = {0, 0, 0}, k, q;
    double Rn[3];
    double eij, r2;
    int N = 0;
    for ( n[2] = 0; n[2] <= 0; ++n[2]){
        for (n[1] = -3; n[1] <= 3; ++n[1]){
            for (n[0] = -3; n[0] <= 3; ++n[0]){
               // if (n[0] == 0 && n[1] == 0 ) continue;
                for (k = 0; k < 3; ++k){
                    Rn[k] = 0;
                    for (q = 0; q < 3; ++q) Rn[k] += abc[3*k+q]*n[q];
                }
                eij = exp(-aijk* CINTsquare_dist(Pij, Rn) ) * fac;
                if (eij > EXPCUTOFF){
                    for (k = 0; k < 3; ++k) v[3*N+k] = Rn[k];
                    v[3*N+3] = eij;
                    ++N;
                }
            }
        }
    }
    *nt = N;
}
void CINTg_nuc_all(double *g, double aij, double *Pij, double *cr, double ak, double fac, CINTEnvVars *envs){
    // fac = exp[-ai*aj/aij*Rij^2]*Ck
    double v[4*7*7*7];
    int gsize = envs->g_size * 3 * ((1<<envs->gbits)+1);
    double* gtemp = malloc(sizeof(double)*gsize);

    int nterms = 0, k;
    double z = 0;
    double Rcn[3];
    double fac2;
    double aijk = aij*ak / (aij + ak);
    for (k = 0; k < 3; ++k) z += envs->abc[2+k*3];

    if (z == 0 )
        SCsum2D(envs->abc, Pij, cr, aijk, v, fac, &nterms);
    else
        SCsum3D(envs->abc, Pij, cr, aijk, v, fac, &nterms);
    
    for (k = 0; k < gsize; ++k) gtemp[k] =0;
    
    for (int n = 0; n < nterms; ++n){
        
        for (k = 0; k < 3; ++k)
            Rcn[k] = cr[k] + v[3*n+k];
        
        CINTg_nuc_lim(gtemp, aij, Pij, Rcn, ak, v[3*n+3], envs);
        
        for (k = 0; k < gsize; ++k)
            g[k] += gtemp[k];
        
    }
    free(gtemp);
}
void CINTg_nuc_lim(double *g, double aij, double *Pij, double *cr, double ak, double fac, CINTEnvVars *envs){
    // fac = exp [-ai*aj/aij * r2ij ] * exp [ -(aij)*(ak)/(aijk)* |P-Rc|^2 ] * Ck
    double aijk = aij + ak;
    const FINT nmax = envs->li_ceil + envs->lj_ceil;
    const FINT lj = envs->lj_ceil;
    const FINT dj = envs->g_stride_j;
    const double *ri = envs->ri;
    const double *rj = envs->rj;
    
    FINT i, j, ptr, k;
    
    double *gx = g;
    double *gy = g + envs->g_size;
    double *gz = g + envs->g_size * 2;
    double p[3], rij[3], x1[3], x2[3];
    double bounds[6], dR;
    double gupper[3], glower[3];
    double saijk = sqrt(aijk);
    double aijk_inv = 0.5/aijk;
    
    for (k = 0; k < 3; ++k) {
        x1[k] = envs->int_lower[k] - ri[k];
        x2[k] = envs->int_upper[k] - ri[k];
        dR = ( aij * Pij[k] + cr[k]*ak ) / aijk;
        p[k] = dR - ri[k];
        rij[k] = ri[k] - rj[k];
        bounds[k] = x1[k] - p[k];
        bounds[k+3] = x2[k] - p[k];
        gupper[k] = exp(-aijk*x2[k]*x2[k]);
        glower[k] = exp(-aijk*x1[k]*x1[k]);
    }
    double fac2 = fac * SQRTPI * PI / (saijk* aijk);
    gx[0] = erf(bounds[0]) - erf(bounds[3]);
    gy[0] = erf(bounds[1]) - erf(bounds[4]);
    gz[0] = (erf(bounds[2]) - erf(bounds[5]))*fac2;
    
    if (nmax > 0){
        gx[1] = aijk_inv*(glower[0] - gupper[0]) + gx[0]*p[0];
        gy[1] = aijk_inv*(glower[1] - gupper[1]) + gy[0]*p[1];
        gz[1] = aijk_inv*(glower[2] - gupper[2])*fac2 + gz[0]*p[2];
    }
    
    for (i = 1; i < nmax; i++){
        for (k = 0; k < 3; ++k){
            glower[k] *= x1[k];
            gupper[k] *= x2[k];
        }
        gx[i+1] = aijk_inv*(i*gx[i-1] + (glower[0] - gupper[0])) + gx[i-1]*p[0];
        gy[i+1] = aijk_inv*(i*gy[i-1] + (glower[1] - gupper[1])) + gy[i-1]*p[1];
        gz[i+1] = aijk_inv*(i*gz[i-1] + (glower[2] - gupper[2])*fac2)+ gz[i-1]*p[2];
    }
    for (j = 1; j <= lj; j++) {
        ptr = dj * j;
        for (i = ptr; i <= ptr + nmax - j; i++) {
            gx[i] = gx[i+1-dj] + rij[0]*gx[i-dj];
            gy[i] = gy[i+1-dj] + rij[1]*gx[i-dj];
            gz[i] = gz[i+1-dj] + rij[2]*gx[i-dj];
        }
    }
}
void CINTg_nuc_far(double *g, double aij, double *rij, double *cr, double ck, double fac, CINTEnvVars *envs){
    
    
    
    
}
*/
