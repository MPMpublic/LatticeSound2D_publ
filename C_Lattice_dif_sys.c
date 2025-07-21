/*
 * LatticeSound 2D open-source 
 *
 * Copyright (C) 2025  Giorgio Lo Presti (MPMpublic)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>  // Include this header to define GSL_SUCCESS
#include <gsl/gsl_matrix.h>
#include "C_Lattice_dif_sys.h"
#include "h_shape.h"
#include <math.h>
#define pi 3.14159265358979323846
// Function for the Jacobian matrix (not necessary for the linear case, but we leave the structure)
int jac(double t, const double y[], double *dfdy, double dfdt[], void *params)
{
    double gamma = *(double *)params;
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 2, 2);
    gsl_matrix *m = &dfdy_mat.matrix;
    
    // Jacobian of the system
    gsl_matrix_set(m, 0, 0, 0.0);  // d(x')/dx = 0
    gsl_matrix_set(m, 0, 1, 1.0);  // d(x')/dv = 1
    gsl_matrix_set(m, 1, 0, -1.0); // d(v')/dx = -omega^2
    gsl_matrix_set(m, 1, 1, -gamma); // d(v')/dv = -gamma

    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    return GSL_SUCCESS;
}

//The real solver
int Osc_smorzato(double t, const double y[], double f[], void *params)
{
    InternalData *internaldata = (InternalData *)params;
    double* gamma = internaldata->gamma_zero;  // Damping coefficient
    double* masses = internaldata->mu;  // Mass
    int index_i=(internaldata->indexes)[0];
    int index_j=(internaldata->indexes)[1];
    double* under=(internaldata->sotto);
    double* upper=(internaldata->sopra);
    double* lf=(internaldata->sn);
    double* rg=(internaldata->dx);
    double* kappas=(internaldata->km_pt);//order: right left down up
    double* v_lf=internaldata->v_sn;
    double* v_rg=internaldata->v_dx;
    double* v_upper=internaldata->v_sopra;
    double* v_under=internaldata->v_sotto;

    int trasmission=(internaldata->transverse); //If it is solid or fluid
    double* rand_interaction=(internaldata->u_rand_vect);
    // printf("Inside diffe system: %i %i, k r l d u %lf %lf %lf %lf\n ",index_i,index_j,kappas[0], kappas[1], kappas[2], kappas[3]);

    //The code evaluates the actual length of the springs, and therefore their direction. Therefore, it must solve the system's 
    // trigonometry to evaluate the correct elongation.
    double x_u, y_u, x_d, y_d, x_l, y_l, x_r, y_r;
    double ipotenusa_u, ipotenusa_d, ipotenusa_l, ipotenusa_r;
    double CosTheta_u, CosTheta_d, CosPhi_l, CosPhi_r, SenTheta_u, SenTheta_d, SenPhi_l, SenPhi_r;
    double deltax_u, deltay_u, deltax_d, deltay_d, deltax_l, deltay_l, deltax_r, deltay_r;
    double deltaSx_u, deltaSy_u, deltaSx_d, deltaSy_d, deltaSx_l, deltaSy_l, deltaSx_r, deltaSy_r;

    x_u=y[0]-upper[0];
    y_u=upper[1]-y[2];
    x_d=y[0]-under[0];
    y_d=y[2]-under[1];
    x_l=y[0]-lf[0];
    y_l=y[2]-lf[1];
    x_r=rg[0]-y[0];
    y_r=y[2]-rg[1];

    ipotenusa_u=sqrt(pow(y_u,2)+pow(x_u,2));
    ipotenusa_d=sqrt(pow(y_d,2)+pow(x_d,2));
    ipotenusa_l=sqrt(pow(y_l,2)+pow(x_l,2));
    ipotenusa_r=sqrt(pow(y_r,2)+pow(x_r,2));

    CosTheta_u=y_u/ipotenusa_u;
    CosTheta_d=y_d/ipotenusa_d;
    CosPhi_l=x_l/ipotenusa_l;
    CosPhi_r=x_r/ipotenusa_r;
    SenTheta_u=x_u/ipotenusa_u;
    SenTheta_d=x_d/ipotenusa_d;
    SenPhi_l=y_l/ipotenusa_l;
    SenPhi_r=y_r/ipotenusa_r;

    ipotenusa_u-=1;
    ipotenusa_d-=1;
    ipotenusa_l-=1;
    ipotenusa_r-=1;

    deltax_u=SenTheta_u*ipotenusa_u;
    deltay_u=CosTheta_u*ipotenusa_u;
    deltax_d=SenTheta_d*ipotenusa_d;
    deltay_d=CosTheta_d*ipotenusa_d;
    deltax_l=CosPhi_l*ipotenusa_l;
    deltay_l=SenPhi_l*ipotenusa_l;
    deltax_r=CosPhi_r*ipotenusa_r;
    deltay_r=SenPhi_r*ipotenusa_r;

    deltaSx_u=abs(y[1]-v_upper[0]);
    deltaSy_u=abs(y[3]-v_upper[1]);
    deltaSx_d=abs(y[1]-v_under[0]);
    deltaSy_d=abs(y[3]-v_under[1]);
    deltaSx_l=abs(y[1]-v_lf[0]);
    deltaSy_l=abs(y[3]-v_lf[1]);
    deltaSx_r=abs(y[1]-v_rg[0]);
    deltaSy_r=abs(y[3]-v_rg[1]);
    //order: right left down up
    f[0] = y[1];  // x' = v
    f[1] = (gamma[0]*abs(rand_interaction[0]*y[1])-(gamma[0] * deltaSx_r + gamma[1] * deltaSx_l + gamma[2] * deltaSx_d + gamma[3] * deltaSx_u) - (trasmission*kappas[3]*deltax_u+trasmission*kappas[2]*deltax_d+kappas[1]*deltax_l-kappas[0]*deltax_r))/masses[0];  // v' = -gamma * v - omega^2 * x 
    f[2] = y[3];  // y' = v
    f[3] = (gamma[0]*abs(rand_interaction[1]*y[3])-(gamma[0] * deltaSy_r + gamma[1] * deltaSy_l + gamma[2] * deltaSy_d + gamma[3] * deltaSy_u) - (-kappas[3]*deltay_u+kappas[2]*deltay_d+trasmission*kappas[1]*deltay_l+trasmission*kappas[0]*deltay_r))/masses[0];  // v' = -gamma * v - omega^2 * x

    return GSL_SUCCESS;
}

//For the forcing, the system does not actually have to solve the differential equation, but read and report the dynamic positions of 
//the boundary oscillators.
int forzante(double t, const double y[], double f[], void *forzanti)
{
    // printf("Inside diffe system: %i %i, k r l d u %lf %lf %lf %lf\n ",index_i,index_j,kappas[0], kappas[1], kappas[2], kappas[3]);

    f[0] = y[1];  // x' = v
    f[1] = ((double *)forzanti)[0];  // v' = -gamma * v - omega^2 * x
    f[2] = y[3];  // y' = v
    f[3] = ((double *)forzanti)[1];  // v' = -gamma * v - omega^2 * x
    return GSL_SUCCESS;
}