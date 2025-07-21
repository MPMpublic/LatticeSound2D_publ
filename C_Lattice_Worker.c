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

//For each thread there is a worker who is responsible for running the simulations and passing the information from the shared memory 
//to the temporary one which will have to be copied with the barriers
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>  // Include this header to define GSL_SUCCESS
#include <gsl/gsl_matrix.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "h_shape.h"
#include "C_Lattice_dif_sys.h"
#define pi 3.14159265358979323846
#define debug 0

//###########################################  Function representing the integration method  ###########################################
void *solver_method(void *arg)
{   
    //----------------------------- Receives data and decodes the structure and work plan ---------------------------
    ThreadData *data = (ThreadData *)arg;
    int id = data->id;
    double t1 = data->tau;    
    double h = data->passo;
    int num_steps_internal = data->num_steps_internal;
    unsigned char* address=data->address;
    unsigned char* info_pt=data->dual_lattice_pt;
    int cham_dim_x=data->cham_dim_x;
    int cham_dim_y=data->cham_dim_y;
    int n_thread=data->n_thread_tot;
    int dimension=data->dimension;
    int riga_riferimento=data->start;
    int fine_riga_riferimento=data->stop;
    double* kxp_pt=data->kxp_position;
    double thermal_coupling=data->thermal_coupling;

    //----------------------------- It must take care of the memory management of its simulator block ---------------------------
    int chamb_x=(cham_dim_x-1), chamb_y=cham_dim_x-1;//On y there are definitely 2 boundaries in 2D

    double up[dimension], down[dimension], left[dimension], right[dimension], center[dimension], mom[dimension];
    double spring[4], mu[4], gamms[4], transverse;
    double S_up[dimension], S_down[dimension], S_left[dimension], S_right[dimension], S_center[dimension],myfloat_plot[dimension];
    double first_h=h, t = 0.0, t_internal, first_fixed_h=h;  // Time variables
    unsigned char* y_lat_interno[1];
    double u_rand;
    double u_rand_vect[2];
    double y[2*dimension];
    double* force_pt_x=data->forcing_x;
    double* force_pt_y=data->forcing_y;
    double* force_pt_z=data->forcing_z;

    double* omega_insidex=data->omega_insidex;
    double* omega_insidey=data->omega_insidey;
    double* omega_insidez=data->omega_insidez;

    double* f0x=data->f0x;
    double* f0y=data->f0y;
    double* f0z=data->f0z;

    int len_omega_insidex=data->len_omega_insidex;
    int len_omega_insidey=data->len_omega_insidey;
    int len_omega_insidez=data->len_omega_insidez;

    int len_f0x=data->len_f0x;
    int len_f0y=data->len_f0y;
    int len_f0z=data->len_f0z;

    double forcing_pointer[2];

    srand48((long)time(NULL) ^ (id << 16));//To initialize the randomizers in a non-fixed way
    
    InternalData argomenti;
    if (debug){
        printf("TEST ON LATTICE WORKER INFO\n");
        printf("id %i righe rif  %i %i\n",id,riga_riferimento,fine_riga_riferimento);
        printf("id %i chamb_x e chamb_y %i %i\n",id,chamb_x,chamb_y);
    }

    //........................... Selecting particles and running RK

    //boundary points
    //For each group of particles, identified as bulk and boundary, there is a different condition on the differential equation. 
    // The lower boundary is the one with the forcing. The other boundaries do not move.
    if(id==0){//lower bound: column changes (index i), j is zero
        for(int i=0;i<cham_dim_x;i++){
            // y_lat_interno[0]=address+(i)*cham_dim_y*2*dimension*sizeof(double);
            center[0]=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)))[0];
            center[1]=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+sizeof(double)))[0];
            S_center[0]=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+2*sizeof(double)))[0];
            S_center[1]=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+3*sizeof(double)))[0];
            if(debug) printf("\n\nINSIDE LATTICE WORKER, center x %lf center y %lf Scenter x %lf Scenter y %lf\n\n",center[0],center[1],S_center[0],S_center[1]);
            //LAUNCHING THE OSCILLATOR WITH FORCING HERE!!!
            forcing_pointer[0]=force_pt_x[i];
            forcing_pointer[1]=force_pt_y[i];
            y[0]=center[0];
            y[1]=S_center[0];
            y[2]=center[1];
            y[3]=S_center[1];

            const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd;  // Runge-Kutta Method
            gsl_odeiv_step *s = gsl_odeiv_step_alloc(T, 4);  // Allocation of the structure for the step
            gsl_odeiv_control *c = gsl_odeiv_control_y_new(1e-7, 1.e-5);  // Accuracy controll: try between 1.e-3 and 1.e-8
            gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc(4);  // Allocation for evolution
            gsl_odeiv_system sys = {forzante, jac, 4, forcing_pointer};  // Definition of ODE system

            // Numerical integration loop: the solve_ivp system is defined and called every time step to ensure a constant number of time steps 
            //so that all oscillators update together
            for ( int b=1; b < num_steps_internal+1; b++)
            {
                // printf("h vale %lf, first h vale %lf\n",h,first_h);
                t_internal=(b)*first_h;
                while (t < t_internal)
                {                    
                    int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t_internal, &h, y);

                    if (status != GSL_SUCCESS)
                        break;
                }
                if (debug) printf("Dati FORZANTE: id:%i, deltatime:%.10e, x:%.10e, speed:%.10e \n", id, t_internal, y[0],y[1]);  // Stampa i risultati
            }
            h=first_fixed_h;
            t=0.0;
            center[0]=y[0];
            S_center[0]=y[1];
            center[1]=y[2];
            S_center[1]=y[3];

            if (debug) printf("IN QUESTA POSIZIONE ij SONO estremo INF Center index: %i %i Center Position:%lf %lf\n",i,0,center[0],center[1]);

            *(kxp_pt+(i)*cham_dim_y*2*dimension) = center[0];
            *(kxp_pt+(i)*cham_dim_y*2*dimension+1) = center[1];
            *(kxp_pt+(i)*cham_dim_y*2*dimension+2) = S_center[0];
            *(kxp_pt+(i)*cham_dim_y*2*dimension+3) = S_center[1];
            gsl_odeiv_evolve_free(e);
            gsl_odeiv_control_free(c);
            gsl_odeiv_step_free(s);
        }
    }
    if(id==n_thread-1){//upper bound: column changes (index i), j is the upper one
        for(int i=0;i<cham_dim_x;i++){

            center[0]=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(cham_dim_y-1)*2*dimension*sizeof(double)))[0];
            center[1]=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(cham_dim_y-1)*2*dimension*sizeof(double)+sizeof(double)))[0];
            S_center[0]=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(cham_dim_y-1)*2*dimension*sizeof(double)+2*sizeof(double)))[0];
            S_center[1]=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(cham_dim_y-1)*2*dimension*sizeof(double)+3*sizeof(double)))[0];

            if (debug) printf("IN QUESTA POSIZIONE ij SONO estremo SUP Center index: %i %i Center Position:%lf %lf\n",i,cham_dim_y-1,center[0],center[1]);
            *(kxp_pt+(i)*cham_dim_y*2*dimension+(cham_dim_y-1)*2*dimension)=center[0];
            *(kxp_pt+(i)*cham_dim_y*2*dimension+(cham_dim_y-1)*2*dimension+1)=center[1];
            *(kxp_pt+(i)*cham_dim_y*2*dimension+(cham_dim_y-1)*2*dimension+2)=S_center[0];
            *(kxp_pt+(i)*cham_dim_y*2*dimension+(cham_dim_y-1)*2*dimension+3)=S_center[1];

        }
    }

    //internal points and lateral bounds
    for(int j=riga_riferimento;j<fine_riga_riferimento;j++){
        for(int i=1;i<chamb_x;i++){//i is the index of x: outer index in loop, innermost in matrix
            
            center[0]=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)))[0];
            center[1]=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+sizeof(double)))[0];
            S_center[0]=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+2*sizeof(double)))[0];
            S_center[1]=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+3*sizeof(double)))[0];

            up[0]=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j+1)*2*dimension*sizeof(double)))[0];
            up[1]=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j+1)*2*dimension*sizeof(double)+sizeof(double)))[0];
            S_up[0]=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j+1)*2*dimension*sizeof(double)+2*sizeof(double)))[0];
            S_up[1]=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j+1)*2*dimension*sizeof(double)+3*sizeof(double)))[0];

            down[0]=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j-1)*2*dimension*sizeof(double)))[0];
            down[1]=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j-1)*2*dimension*sizeof(double)+sizeof(double)))[0];
            S_down[0]=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j-1)*2*dimension*sizeof(double)+2*sizeof(double)))[0];
            S_down[1]=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j-1)*2*dimension*sizeof(double)+3*sizeof(double)))[0];

            left[0]=((double*)(address+(i-1)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)))[0];
            left[1]=((double*)(address+(i-1)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+sizeof(double)))[0];
            S_left[0]=((double*)(address+(i-1)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+2*sizeof(double)))[0];
            S_left[1]=((double*)(address+(i-1)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+3*sizeof(double)))[0];

            right[0]=((double*)(address+(i+1)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)))[0];
            right[1]=((double*)(address+(i+1)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+sizeof(double)))[0];
            S_right[0]=((double*)(address+(i+1)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+2*sizeof(double)))[0];
            S_right[1]=((double*)(address+(i+1)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+3*sizeof(double)))[0];

            mu[0]=((double*)(info_pt+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)))[0];//RIGHT... Che dovrebbe essere quella dell'elemento
            mu[1]=((double*)(info_pt+(i-1)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)))[0];//LEFT
            mu[2]=((double*)(info_pt+(i)*cham_dim_y*2*dimension*sizeof(double)+(j-1)*2*dimension*sizeof(double)))[0];//DOWN
            mu[3]=((double*)(info_pt+(i)*cham_dim_y*2*dimension*sizeof(double)+(j+1)*2*dimension*sizeof(double)))[0];//UP
            if (debug) printf("Len indici: %i %i, k r %lf l %lf d %lf u %lf \n",i,j,mu[0],mu[1],mu[2],mu[3]);
            spring[0]=((double*)(info_pt+sizeof(double)+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)))[0];//RIGHT
            spring[1]=((double*)(info_pt+sizeof(double)+(i-1)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)))[0];//LEFT
            spring[2]=((double*)(info_pt+sizeof(double)+(i)*cham_dim_y*2*dimension*sizeof(double)+(j-1)*2*dimension*sizeof(double)))[0];//DOWN
            spring[3]=((double*)(info_pt+sizeof(double)+(i)*cham_dim_y*2*dimension*sizeof(double)+(j+1)*2*dimension*sizeof(double)))[0];//UP
            if (debug) printf("Molle indici: %i %i, k r %lf l %lf d %lf u %lf \n",i,j,spring[0],spring[1],spring[2],spring[3]);
            gamms[0]=((double*)(info_pt+2*sizeof(double)+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)))[0];//RIGHT
            gamms[1]=((double*)(info_pt+2*sizeof(double)+(i-1)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)))[0];//LEFT
            gamms[2]=((double*)(info_pt+2*sizeof(double)+(i)*cham_dim_y*2*dimension*sizeof(double)+(j-1)*2*dimension*sizeof(double)))[0];//DOWN
            gamms[3]=((double*)(info_pt+2*sizeof(double)+(i)*cham_dim_y*2*dimension*sizeof(double)+(j+1)*2*dimension*sizeof(double)))[0];//UP
            if (debug) printf("gamma indici: %i %i, k r %lf l %lf d %lf u %lf",i,j,gamms[0],gamms[1],gamms[2],gamms[3]);
            transverse=((double*)(info_pt+3*sizeof(double)+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)))[0];
            if (debug) printf("Transverse indici: %i %i, %lf\n",i,j,transverse);

            //ACQUISITION OF PHOTONS TO EXCITE THE SYSTEM
            if (thermal_coupling>0.001)//TO MAKE THE AUTOMATIC OPTION IN 3D
            {
                u_rand=rand() / (RAND_MAX + 1.0);
                u_rand_vect[0]=thermal_coupling*cos(u_rand);//speed normalized to 1
                u_rand_vect[1]=thermal_coupling*sin(u_rand);
            }
            else
            {
                u_rand_vect[0]=0;
                u_rand_vect[1]=0;
            }

            y[0]=center[0];
            y[1]=S_center[0];
            y[2]=center[1];
            y[3]=S_center[1];

            //Defining the variables to be passed to the differential solver
            int* indici=malloc(dimension*sizeof(int));
            indici[0]=i;
            indici[1]=j;
            
            argomenti.sn=left;
            argomenti.dx=right;
            argomenti.sopra=up;
            argomenti.sotto=down;
            argomenti.v_sn=S_left;
            argomenti.v_dx=S_right;
            argomenti.v_sopra=S_up;
            argomenti.v_sotto=S_down;
            argomenti.gamma_zero=gamms;
            argomenti.km_pt=spring;
            argomenti.indexes=indici;
            argomenti.mu=mu;
            argomenti.transverse=transverse;
            argomenti.u_rand_vect=u_rand_vect;

            const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd;  // Runge-Kutta Method
            gsl_odeiv_step *s = gsl_odeiv_step_alloc(T, 4);  // Allocation of the structure for the step
            gsl_odeiv_control *c = gsl_odeiv_control_y_new(1e-7, 1.e-5);  // Accuracy controll: try between 1.e-3 and 1.e-8
            gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc(4);  // Allocation for evolution
            gsl_odeiv_system sys = {Osc_smorzato, jac, 4, (void *)&argomenti};  // Definition of ODE system

            // Numerical integration loop
            for ( int b=1; b < num_steps_internal+1; b++)
            {
                t_internal=(b)*first_h;
                while (t < t_internal)
                {                    
                    int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t_internal, &h, y);

                    if (status != GSL_SUCCESS)
                        break;
                }

                if (debug) printf("id:%i, deltatime:%.10e, csi:%.10e, speed:%.10e \n", id, t_internal, y[0],y[1]);  // Print results
            }
            h=first_fixed_h;
            t=0.0;
            center[0]=y[0];
            S_center[0]=y[1];
            center[1]=y[2];
            S_center[1]=y[3];
            
            //Passing data back
            *(kxp_pt+(i)*cham_dim_y*2*dimension+(j)*2*dimension)=center[0];
            *(kxp_pt+(i)*cham_dim_y*2*dimension+(j)*2*dimension+1)=center[1];
            *(kxp_pt+(i)*cham_dim_y*2*dimension+(j)*2*dimension+2)=S_center[0];
            *(kxp_pt+(i)*cham_dim_y*2*dimension+(j)*2*dimension+3)=S_center[1];
            //printf("Nella locazione del kx del worker dopo rk c'è %lf\n",*(kxp_pt+(i+1+riga_riferimento)*cham_dim_y*2*dimension*sizeof(double)+(j+colonna_riferimento+1)*2*dimension*sizeof(double)));
            if (debug) printf("Nella locazione del kx del worker %i %i dopo rk c'è %.10lf %.10lf %.10lf %.10lf\n",i,j,center[0],center[1],S_center[0],S_center[1]);
            // Clean allocated structures
            gsl_odeiv_evolve_free(e);
            gsl_odeiv_control_free(c);
            gsl_odeiv_step_free(s);
            free(indici);
        }
    //LATERAL BOUNDARY
        //left bound

        y_lat_interno[0]=address+(j)*2*dimension*sizeof(double);
        memcpy(&center[0], y_lat_interno[0], sizeof(double));
        memcpy(&center[1], y_lat_interno[0]+sizeof(double), sizeof(double));
        if(debug) printf("IN QUESTA POSIZIONE ij SONO estremo SN Center index: %i %i Center Position:%lf %lf\n",0,j,center[0],center[1]);

        memcpy(&S_center[0], &center[0]+2*sizeof(double), sizeof(double));
        memcpy(&S_center[1], &center[1]+2*sizeof(double), sizeof(double));
        
        // printf("Nella locazione del kx del worker prim rk c'è %lf\n",*(kxp_pt+(i+1+riga_riferimento)*cham_dim_y*2*dimension*sizeof(double)+(j+colonna_riferimento+1)*2*dimension*sizeof(double)));
        // printf("%li",(i+1+riga_riferimento)*cham_dim_y*2*dimension*sizeof(double)+(j+colonna_riferimento+1)*2*dimension*sizeof(double));
        *(kxp_pt+(j)*2*dimension)=center[0];
        *(kxp_pt+(j)*2*dimension+1)=center[1];
        *(kxp_pt+(j)*2*dimension+2)=S_center[0];
        *(kxp_pt+(j)*2*dimension+3)=S_center[1];

        //right bound
        y_lat_interno[0]=address+(cham_dim_x-1)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double);
        memcpy(&center[0], y_lat_interno[0], sizeof(double));
        memcpy(&center[1], y_lat_interno[0]+sizeof(double), sizeof(double));
        if(debug) printf("IN QUESTA POSIZIONE ij SONO estremo DX Center index: %i %i Center Position:%lf %lf\n",cham_dim_x-1,j,center[0],center[1]);

        memcpy(&S_center[0], &center[0]+2*sizeof(double), sizeof(double));
        memcpy(&S_center[1], &center[1]+2*sizeof(double), sizeof(double));
         
        // printf("Nella locazione del kx del worker prim rk c'è %lf\n",*(kxp_pt+(i+1+riga_riferimento)*cham_dim_y*2*dimension*sizeof(double)+(j+colonna_riferimento+1)*2*dimension*sizeof(double)));
        // printf("%li",(i+1+riga_riferimento)*cham_dim_y*2*dimension*sizeof(double)+(j+colonna_riferimento+1)*2*dimension*sizeof(double));
        *(kxp_pt+(cham_dim_x-1)*cham_dim_y*2*dimension+(j)*2*dimension)=center[0];
        *(kxp_pt+(cham_dim_x-1)*cham_dim_y*2*dimension+(j)*2*dimension+1)=center[1];
        *(kxp_pt+(cham_dim_x-1)*cham_dim_y*2*dimension+(j)*2*dimension+2)=S_center[0];
        *(kxp_pt+(cham_dim_x-1)*cham_dim_y*2*dimension+(j)*2*dimension+3)=S_center[1];

    }
}
