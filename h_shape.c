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
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

double h_shape_g(double mass, double kappa_0, double eff_k_const)//double csi csi is the damping ratio,int n_basis n_basis: number of dinstinct damping ratios
{
    double omega,omega_p;
    omega_p=pow(eff_k_const/mass,0.5);
    omega=pow(kappa_0/mass,0.5);
    return 2*omega_p*omega/(pow(omega_p,2)+pow(omega,2));//Qui ho corretto un * con un + al denominatore
}

int h_shape_f(int force_points,double* force_pt, double omega_inside[], int freq_numbers,double f0[], int f0_len, double time)
{
    int i,j;//sizeof(omega_inside);
    // printf("freq_numbers è %i\n", freq_numbers);
    
    for(i=0;i<force_points;i++)
    {
        force_pt[i]=0;
        for(j=0;j<freq_numbers;j++)
        {
            force_pt[i]+=f0[j]*cos(omega_inside[j]*time);
        }
    }
    return 0;
}


//Old pieces of code
//#include <malloc.h>
// double d;
// pt = (double *) malloc(sizeof(double) * 3*a*b*c);
// if(pt == NULL)
// {
// printf("WARNING: Memoria esaurita in chiamata malloc k\n");
// exit(1);
// }
// *pt=d;
// printf("turiiii %lf \n",*pt);

// d=57.98;
// pt[0]=d;
// pt[1]=2.6;
// printf("turi %lf %lf\n",pt[0],pt[1]);


// int h_shape_k(int a, int b, int c, double* pt,double km)
// {
//     int g;
//     for(g=0;g<3*a*b*c;g++){
//         pt[g]=km;
//     }
//     return 0;
// }