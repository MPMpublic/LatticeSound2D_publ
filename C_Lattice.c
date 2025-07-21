/*This is the code that divides the work among the workers having received the info from LSM. Each worker will have his task called 
  assignment. The assignment starts in this same code because we need to synchronize the threads in reading and writing the shared memory 
  portions. The physical command that must execute the code is solver_method.
  So C_lattice must first interpret the forcings, then divide the work (in the main) and then launch the assignments.*/

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

//###################################################   LIBRARIES   ###################################################
#include <sys/ipc.h>		                            /* for system's IPC_xxx definitions */
#include <sys/shm.h>		                            /* for shmget, shmat, shmdt, shmctl */
#include <sys/sem.h>		                            /* for semget, semctl, semop */

#include <stdlib.h> 
#include <stdio.h> 
#include <errno.h> 
#include <unistd.h>                                     /* Contains Sleep */
#include <string.h>
#include <time.h>

#include <math.h>
#include <limits.h>
#include <malloc.h>

#include <pthread.h>                                     /*For multithread */

#include "h_shape.h"
#include "C_Lattice_dif_sys.h"
#include "C_Lattice_Worker.h"

#define pi 3.14159265358979323846
#define debug 0
#define MAX_NUM_FLOAT 100

long long counter = 0;                                    // Counter shared between workers
pthread_mutex_t mutex;                                    // Mutex to protect access to the counter in multithread
pthread_barrier_t barrier;                                // The barrier for synchronizing threads

//###################################################   FUNCTIONS   ###################################################
//------------------------------------------  FUNCTION TO READ THE FORCING  -------------------------------------------
void parse_vector_string(const char *input, double *output, int *count) {
    char buffer[1024];
    strncpy(buffer, input, sizeof(buffer));
    buffer[sizeof(buffer) - 1] = '\0';                                               // It must be null-terminated

    // Removes square brackets
    if (buffer[0] == '[') memmove(buffer, buffer + 1, strlen(buffer));                // remove '['
    if (buffer[strlen(buffer) - 1] == ']') buffer[strlen(buffer) - 1] = '\0';         // remove ']'

    // Tokenizzazione
    char *token = strtok(buffer, ",");
    *count = 0;

    while (token != NULL && *count < MAX_NUM_FLOAT) {
        output[*count] = strtof(token, NULL);
        (*count)++;
        token = strtok(NULL,",");
    }
}

//------------------------------------------  TASK THAT THREADS WILL EXECUTE  -------------------------------------------
/*Furthermore, each output_time must, for each worker, give us global and central oscillator information.*/
/*The workers must run all together for each step. When they work, they temporarly use the kxp_pt matrix which contains positions
and speed as the shared one (identified by "address" pointer). When they save, they must copy from  kxp_pt to address. Remember
that kxp_pt and address are two different kind of pointers (double and char), thus the algebra of pointers must contain this info
during operations*/
void *assignment(void *arg) {
    FILE *file;
    ThreadData *data = (ThreadData *)arg;                                              // Receives and transmits the data structure
    /*Time and id variables*/
    int id = data->id;
    double h = data->passo;
    double hinitial=h;

    /*Global analysis info*/
    double position_element_x, position_element_y, position_element_dx, 
           initial_len_x, initial_len_y, speed_element_x, speed_element_y;
    double deltaVx, deltaVy;
    double averagex, msqx, averagey, msqy, average_speed_x, v_msqx, average_speed_y, 
           v_msqy, deltaVxaverage, deltaVyaverage, msq_deltaVx, msq_deltaVy;
    int count_average;
    double* force_field;
    double** force_field_matrix;
    double *positionx, *positiony, *speed_x,*speed_y;

    /*Data struct variables*/
    int num_steps_ext = data->num_steps_ext;
    double* kxp_pt=data->kxp_position;//Momentary positions and velocities before re-updating the reticle
    unsigned char* address=data->address;
    unsigned char* dual_lattice_pt=data->dual_lattice_pt;
    int dimension=data->dimension; 
    int cham_dim_x=data->cham_dim_x;
    int cham_dim_y=data->cham_dim_y;
    int n_thread=data->n_thread_tot;
    int riga_start=data->start;
    int riga_stop=data->stop;
    ParticleData* data_array = (ParticleData*)address;

    
    int chamb_x=(cham_dim_x-1), chamb_y;                                            //Size information inside this code
    int x_sn_plot=data->start_plotX,                             //Choose the element to plot
        x_dx_plot=data->stop_plotX,
        y_giu_plot=data->start_plotY, 
        y_su_plot=data->stop_plotY;

    /*Data struct variables fon analysis*/
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

    int output_time=data->output_time;
    char* nome_file=data->nome_file;

    /*Control variable*/
    int check;

    /*Prepare the file name*/
    char risultato[256];
    #ifdef _WIN32 
        char end1[]="'\'Trd";
    #else  
        char end1[]="/Trd";
    #endif
    char num[20];
    char num2[100];
    char end2[]="_t";
    char txt[]=".txt";
    sprintf(num, "%d", id);
    int i,j;
    for(i=0;i<strlen(nome_file)+1;i++){
        risultato[i]=nome_file[i];
    }
    for(i=0;i<strlen(end1);i++){
        risultato[strlen(nome_file)+i]=end1[i];
    }
    for(i=0;i<strlen(num);i++){
        risultato[strlen(nome_file)+4+i]=num[i];
    }
    for(i=0;i<strlen(end2);i++){
        risultato[strlen(nome_file)+4+strlen(num)+i]=end2[i];
    }
    risultato[strlen(nome_file)+4+strlen(num)+strlen(end2)]='\0';
    if (debug) printf("File name is: %s",risultato);

    //************************************* SELECT THE STRIP OF THE LATTICE FOR THE WORKER ******************************
    /*Net count of strip dimension*/
    if(n_thread==1) chamb_y=(cham_dim_y-2);
    if(n_thread!=1) {
        chamb_y=((int)(floor((cham_dim_y-2)/n_thread)));
    }
    if((id==0 || id==n_thread-1) && n_thread!=1) chamb_y=(int)(floor((cham_dim_y-2)/n_thread+1)-1);
    if(id==n_thread-2 && n_thread!=1){
        if((cham_dim_y-2)%n_thread!=0) chamb_y=chamb_y+(cham_dim_y-2)%n_thread;
    }
    if (debug) printf("\n");

    if (debug)
    {
        printf("Worker id %i cham_dim_x %i, cham_dim_y %i\n", id, cham_dim_x, cham_dim_y);
        printf("Worker id %i chamb_x %i, chamb_y %i\n", id, chamb_x, chamb_y);
        printf("Worker id %i riga_start %i, riga_stop %i\n", id, riga_start, riga_stop);
    }
    if (debug) printf("\n");
    if (debug){
        printf("Central data before algorithm at  thread: %i\n",id);
        for(j=riga_start;j<riga_stop;j++){
            for(i=1;i<chamb_x;i++){//i è indice della x: indice esterno in ciclo, più interno in matrice   
                ParticleData* p = &data_array[i * cham_dim_y + j];
                double sx = p->pos_x;
                double sy = p->pos_y;
                double vx = p->speed_x;
                double vy = p->speed_y;
                printf("Posizione x %lf y %lf e velocità %lf %lf\n",sx,sy,vx,vy);
            }
        }
        if (debug) printf("\n");
    }
    //*****************************************************   CYCLE ON TIME   **********************************************
    printf("Cycle on time started\n");
    for (int step = 0; step < num_steps_ext; step++) {
        //----------------------- PREPARE THE FORCING AND THEN MAKE IT AVAILABLE FOR THE WORKERS ---------------------------------
        //This is redundant and only serves the lower boundary (and then remains the same) so only the first worker does it
        if(id==0){
            check=h_shape_f(chamb_x,force_pt_x,omega_insidex,len_omega_insidex,f0x,len_f0x,step*hinitial);//last is time, first is the border dimension
            if (check!=0) printf("h_shape f on x check failure with code: %i \n", check);

            check=h_shape_f(chamb_x,force_pt_y,omega_insidey,len_omega_insidey,f0y,len_f0y,step*hinitial);//last is time
            if (check!=0) printf("h_shape f on y check failure with code: %i \n", check);

            check=h_shape_f(chamb_x,force_pt_z,omega_insidez,len_omega_insidez,f0z,len_f0z,step*hinitial);//last is time
            if (check!=0) printf("h_shape f on z check failure with code: %i \n", check);
            if(debug){
                printf("force_pt_x of 0 at step %i: %lf \n", step,force_pt_x[0]);
                printf("force_pt_y of 0 at step %i: %lf \n", step,force_pt_y[0]);
                printf("force_pt_z of 0 at step %i: %lf \n", step,force_pt_z[0]);
            }
        }
    
        //--------------------------------  Wait until the first worker has calculated the forcing  ---------------------------------
        pthread_barrier_wait(&barrier);// You need to make sure everyone has read and is ready before you launch operations
        //----------------------- Explain and launch the work that each thread will perform in the time loop  ------------------
        if (debug) printf("Thread %d esegue incarico al passo %d\n", id, step + 1);
        solver_method((void *)data);
    
        //In the meantime the code has written the new positions on the kx and kp but has not yet put them on the lattice because the others
        //processors must wait until everyone is finished before copying onto the common sheet
        // In caso di debug: Simula un piccolo lavoro (ad esempio una pausa) usleep(500000);  // Pausa di 0.5 secondi per simulare un lavoro
        
        //--------------------------------  Synchronize threads at the end of each step   ---------------------------------
        pthread_barrier_wait(&barrier);// At the end of the barrier, all threads will copy the data

        //--------------------------------    Lattice update  ---------------------------------
                                /*As the loop updates, global calculations are performed.*/
        if (debug) printf("Aggiornamento reticolo interno\n");

        //CALCULATING OUTPUT VALUES FOR ANALYSIS
        //For data analysis: since the various threads already have the locations to work on, they can take care of them. However, since 
        //they must pass the data to shared memory, they can perform data analysis during the data read cycle, but they must do so at the 
        //right time interval to avoid clogging up the work. Therefore, there are two copying cycles: one without analysis and one with 
        //analysis. The boundaries are not analyzed. Furthermore, the analysis is performed on position, velocity, and force. For all three,
        // a double cycle must be performed to evaluate the deviation. The problematic one is force because the difference in velocity 
        //between the non-updated lattice and the updating lattice is evaluated, given its definition: delta v over delta t.

        //FORCE ANALYSIS
        if (step%output_time==0){
            //Before updating I create a temporary lattice to calculate the delta v, so I have the forces
            //It could have been done on the spot, however it appears to corrupt the data.
            count_average=0;
            deltaVxaverage=0;
            deltaVyaverage=0;

            force_field_matrix = (double**)malloc((riga_stop-riga_start) * sizeof(double*));
            if(force_field_matrix==NULL) printf("WARNING FORCE FIELD MATRIX ALLOCATION FAILURE");

            for(j=riga_start;j<riga_stop;j++){
                force_field_matrix[j - riga_start] = (double*)malloc(chamb_x * dimension * sizeof(double));
                if (force_field_matrix[j - riga_start] == NULL) printf("WARNING: row %d allocation failed\n", j);
                for(i=1;i<chamb_x;i++){//i is the index of x: outer index in loop, innermost in matrix  
                    ParticleData* p = &data_array[i * cham_dim_y + j];
                    double sx = p->pos_x;
                    double sy = p->pos_y;
                    double vx = p->speed_x;
                    double vy = p->speed_y;
                    if (debug) printf("Posizione x %lf y %lf e velocità %lf %lf",sx,sy,vx,vy);
                    speed_x=(kxp_pt+(i)*cham_dim_y*2*dimension+(j)*2*dimension+2);
                    speed_y=(kxp_pt+(i)*cham_dim_y*2*dimension+(j)*2*dimension+3);

                    deltaVx=speed_x[0]-*((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+2*sizeof(double)));
                    deltaVy=speed_y[0]-*((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+3*sizeof(double)));

                    if (debug) printf("deltaVx VALE: %.8lf\n",deltaVx);
                    force_field_matrix[j - riga_start][2*i] = deltaVx;
                    force_field_matrix[j - riga_start][2*i + 1] = deltaVy;

                    deltaVxaverage+=deltaVx;
                    deltaVyaverage+=deltaVy;
                    count_average++;
                }

            }
            deltaVxaverage=deltaVxaverage/count_average;
            deltaVyaverage=deltaVyaverage/count_average;
            msq_deltaVx=0;
            msq_deltaVy=0;
            if (debug){
                for(j=riga_start;j<riga_stop;j++){
                    for(i=1;i<chamb_x;i++){
                        printf("Forcefield matrix j:%i i:%i: %.10lf\n",j,i,force_field_matrix[j-riga_start][i]);
                    }
                }
            }
        //NEL CASO SI DOVESSERO METTERE DELLE BOUNDARY LATERALI DIVERSE ANDREBBE AGGIORNATO ANCHE QUI!
            for(j=riga_start;j<riga_stop;j++){
                for(i=1;i<chamb_x;i++){//i is the index of x: outer index in loop, innermost in matrix  

                    positionx=(kxp_pt+(i)*cham_dim_y*2*dimension+(j)*2*dimension);
                    positiony=(kxp_pt+(i)*cham_dim_y*2*dimension+(j)*2*dimension+1);
                    speed_x=(kxp_pt+(i)*cham_dim_y*2*dimension+(j)*2*dimension+2);
                    speed_y=(kxp_pt+(i)*cham_dim_y*2*dimension+(j)*2*dimension+3);

                    deltaVx=speed_x[0]-*((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+2*sizeof(double)));
                    deltaVy=speed_y[0]-*((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+3*sizeof(double)));
                    msq_deltaVx+=pow(deltaVx-deltaVxaverage,2);
                    msq_deltaVy+=pow(deltaVy-deltaVyaverage,2);

                    *((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)))=positionx[0];
                    *((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+sizeof(double)))=positiony[0];
                    *((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+2*sizeof(double)))=speed_x[0];
                    *((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+3*sizeof(double)))=speed_y[0];
                }
            }
            msq_deltaVx=msq_deltaVx/count_average;
            msq_deltaVy=msq_deltaVy/count_average;
            // printf("msq_deltaVx: %.15lf\n",msq_deltaVx);
            sprintf(num2, "%d", step);
            for(i=0;i<strlen(num2);i++){
                risultato[strlen(nome_file)+4+strlen(num)+2+i]=num2[i];
            }
            for(i=0;i<strlen(txt);i++){
                risultato[strlen(nome_file)+4+strlen(num)+2+strlen(num2)+i]=txt[i];
            }
            printf("Location del risultato %s\n",risultato);
            file = fopen(risultato, "w");

            if (file == NULL) {
                printf("Error in opening result file %s\n",risultato);
                return NULL;
            }
            averagex=0;
            averagey=0;
            average_speed_x=0;
            average_speed_y=0;
            count_average=0;
            for(j=riga_start;j<riga_stop;j++){
                for(int i=1;i<chamb_x;i++){//i is the index of x: outer index in loop, innermost in matrix
                    position_element_x=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)))[0];

                    initial_len_x=i;
                    averagex+=position_element_x-initial_len_x;

                    position_element_y=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+sizeof(double)))[0];
                    initial_len_y=j;
                    averagey+=position_element_y-initial_len_y;

                    speed_element_x=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+2*sizeof(double)))[0];
                    average_speed_x+=speed_element_x;

                    speed_element_y=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+3*sizeof(double)))[0];
                    average_speed_y+=speed_element_y;

                    count_average++;
                }
            } 
            averagex=averagex/count_average;
            averagey=averagey/count_average;
            average_speed_x=average_speed_x/count_average;
            average_speed_y=average_speed_y/count_average;
            msqx=0;
            msqy=0;
            v_msqx=0;
            v_msqy=0;
            count_average=0;

            for(j=riga_start;j<riga_stop;j++){
                for(int i=1;i<chamb_x;i++){//i is the index of x: outer index in loop, innermost in matrix
                    position_element_x=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)))[0];
                    position_element_y=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+sizeof(double)))[0];
                    initial_len_x=i;
                    initial_len_y=j;
                    speed_element_x=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+2*sizeof(double)))[0];
                    speed_element_y=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+3*sizeof(double)))[0];
                    msqx+=pow((position_element_x-initial_len_x-averagex),2);
                    msqy+=pow((position_element_y-initial_len_y-averagey),2);
                    v_msqx+=pow((speed_element_x-average_speed_x),2);
                    v_msqy+=pow((speed_element_y-average_speed_y),2);
                    count_average++;
                }
            }
            msqx=msqx/count_average;
            msqy=msqy/count_average;
            v_msqx=v_msqx/count_average;
            v_msqy=v_msqy/count_average;
            printf("Average x vale: %lf \n",averagex);
            printf("MSQx vale: %lf \n",msqx);
            
            fprintf(file, "%.15lf %.15lf\n", averagex, averagey);
            fprintf(file, "%.15lf %.15lf\n", msqx, msqy);
            fprintf(file, "%.15lf %.15lf\n", average_speed_x, average_speed_y);
            fprintf(file, "%.15lf %.15lf\n", v_msqx, v_msqy);
            fprintf(file, "%.15lf %.15lf\n", deltaVxaverage, deltaVyaverage);
            fprintf(file, "%.15lf %.15lf\n", msq_deltaVx, msq_deltaVy);

            //CALCULATING OUTPUT VALUES FOR THE PLOT
            //printf("CHECK SU RIGHE RIFERIMENTO %i %i %i %i\n",riga_riferimento+1,y_giu_plot, riga_riferimento+chamb_y+1, y_su_plot);
            for(j=riga_start;j<riga_stop;j++){
                if(j>=y_giu_plot && j<y_su_plot){
                    for(int i=x_sn_plot;i<x_dx_plot;i++){//i is the index of x: outer index in loop, innermost in matrix
                        position_element_x=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)))[0];
                        position_element_y=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+sizeof(double)))[0];
                        speed_element_x=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+2*sizeof(double)))[0];
                        speed_element_y=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+3*sizeof(double)))[0];
                        fprintf(file, "%i %i %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n", i,j,position_element_x, position_element_y, speed_element_x, speed_element_y, force_field_matrix[j - riga_start][2*i],force_field_matrix[j - riga_start][2*i+1]);//force_field_matrix[j-riga_start][2*i],force_field_matrix[j-riga_start][2*i+1]);
                        if (debug) printf("NEL FILE: %i %i %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n", i,j,position_element_x, position_element_y, speed_element_x, speed_element_y, force_field_matrix[j - riga_start][2*i],force_field_matrix[j - riga_start][2*i+1]);// force_field_matrix[j-riga_start][2*i],force_field_matrix[j-riga_start][2*i+1]);
                    }
                }
            }
            //PLOT FOR THE UPPER AND LOWER EXTREMES
            if(id==0){
                j=0;
                if(j>=y_giu_plot && j<=y_su_plot){
                    for(int i=x_sn_plot;i<x_dx_plot;i++){//i is the index of x: outer index in loop, innermost in matrix
                        position_element_x=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)))[0];
                        position_element_y=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+sizeof(double)))[0];
                        speed_element_x=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+2*sizeof(double)))[0];
                        speed_element_y=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+3*sizeof(double)))[0];
                        fprintf(file, "%i %i %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n", i,j,position_element_x, position_element_y, speed_element_x, speed_element_y, 0.0,0.0);
                    }
                }
            }
            if(id==n_thread-1){
                j=riga_stop;
                if(j>=y_giu_plot && j<=y_su_plot){
                    for(int i=x_sn_plot;i<x_dx_plot;i++){//i is the index of x: outer index in loop, innermost in matrix
                        position_element_x=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)))[0];
                        position_element_y=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+sizeof(double)))[0];
                        speed_element_x=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+2*sizeof(double)))[0];
                        speed_element_y=((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+3*sizeof(double)))[0];
                        fprintf(file, "%i %i %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n", i,j,position_element_x, position_element_y, speed_element_x, speed_element_y, 0.0,0.0);
                    }
                } 
            }
            for(j=riga_start;j<riga_stop;j++){
                free(force_field_matrix[j - riga_start]);
            }
            free(force_field_matrix);
            fclose(file);
        //NEL CASO SI DOVESSERO METTERE DELLE BOUNDARY LATERALI DIVERSE ANDREBBE AGGIORNATO ANCHE QUI!
        }
        else{
            for(j=riga_start;j<riga_stop;j++){
                for(i=1;i<chamb_x;i++){//i is the index of x: outer index in loop, innermost in matrix   

                    positionx=(kxp_pt+(i)*cham_dim_y*2*dimension+(j)*2*dimension);
                    positiony=(kxp_pt+(i)*cham_dim_y*2*dimension+(j)*2*dimension+1);
                    speed_x=(kxp_pt+(i)*cham_dim_y*2*dimension+(j)*2*dimension+2);
                    speed_y=(kxp_pt+(i)*cham_dim_y*2*dimension+(j)*2*dimension+3);

                    *((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)))=positionx[0];
                    *((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+sizeof(double)))=positiony[0];
                    *((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+2*sizeof(double)))=speed_x[0];
                    *((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+3*sizeof(double)))=speed_y[0];
                }
            }
        //NEL CASO SI DOVESSERO METTERE DELLE BOUNDARY LATERALI DIVERSE ANDREBBE AGGIORNATO ANCHE QUI!
        }
        //LOWER BOUNDARY UPDATE
        if(id==0){
            for(i=1;i<chamb_x;i++){//i is the index of x: outer index in loop, innermost in matrix
                j=0;//SET THE FIRST ROW  
                positionx=(kxp_pt+(i)*cham_dim_y*2*dimension+(j)*2*dimension);
                positiony=(kxp_pt+(i)*cham_dim_y*2*dimension+(j)*2*dimension+1);
                speed_x=(kxp_pt+(i)*cham_dim_y*2*dimension+(j)*2*dimension+2);
                speed_y=(kxp_pt+(i)*cham_dim_y*2*dimension+(j)*2*dimension+3);

                *((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)))=positionx[0];
                *((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+sizeof(double)))=positiony[0];
                *((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+2*sizeof(double)))=speed_x[0];
                *((double*)(address+(i)*cham_dim_y*2*dimension*sizeof(double)+(j)*2*dimension*sizeof(double)+3*sizeof(double)))=speed_y[0];
            }
        }
        
        
    //--------------------------------  Synchronize threads: all workers must copy before new cycle   ---------------------------------
        pthread_barrier_wait(&barrier);// Once the barrier is over, a new time cycle can be made because everyone will read the new data.
        printf("Syncro before copy\n");
    }

    return NULL;
}
/*#####################################################################################################################################*/
/*###############################################            MAIN CODE            ###############################################*/
/*Main code is the interface of C_Lattice with LSM. It must organize the rows for the assignements function (the workers are defined
through thread_data). First it must open the shared memory (both direct and dual lattice), then pass the keys and the main data.*/
int main(int argc, char **argv)
{
    //##########################################    Definition of variables   ##################################################
    //Reading the variables passed by the manager
    int chamb_x=atoi(argv[1]), chamb_y=atoi(argv[2]), chamb_z=atoi(argv[3]),D=atoi(argv[4]),t_steps=atoi(argv[6]),NUM_THREADS=atoi(argv[7]), output_time=atoi(argv[14]);
    int extr_start_plot_x=atoi(argv[16]), extr_start_plot_y=atoi(argv[18]), extr_start_plot_z=atoi(argv[20]);
    int extr_end_plot_x=atoi(argv[17]), extr_end_plot_y=atoi(argv[19]), extr_end_plot_z=atoi(argv[21]);
    int i,j;
    double f0x[MAX_NUM_FLOAT], f0y[MAX_NUM_FLOAT], f0z[MAX_NUM_FLOAT];
    double omega_insidex[MAX_NUM_FLOAT], omega_insidey[MAX_NUM_FLOAT], omega_insidez[MAX_NUM_FLOAT];
    int count_forcing_x, count_forcing_y, count_forcing_z, count_forcing_omega_x, count_forcing_omega_y, count_forcing_omega_z;
    
    //FORCING VARIABLES
    parse_vector_string(argv[8], f0x, &count_forcing_x);
    parse_vector_string(argv[9], omega_insidex, &count_forcing_omega_x);
    parse_vector_string(argv[10], f0y, &count_forcing_y);
    parse_vector_string(argv[11], omega_insidey, &count_forcing_omega_y);
    parse_vector_string(argv[12], f0z, &count_forcing_z);
    parse_vector_string(argv[13], omega_insidez, &count_forcing_omega_z);

    //Shared memory variables
    int sem_id,shm_id, check, force_points=chamb_x;    
    int foo = 16;
    unsigned char* address = (unsigned char*)&foo;
    unsigned char* dual_lattice_pt = (unsigned char*)&foo;

    //Pointers for dynamics
    double *force_pt_x, *force_pt_y, *force_pt_z;//The workers' leader defines the lattice constants and the forcings
    double *kxp_pt;//Matrix that workers will temporarily work on before updating the lattice again
    kxp_pt = (double *) malloc(sizeof(double) * chamb_x*chamb_y*2*D);
    if(kxp_pt == NULL)
    {
        printf("Out of memory in C_Lattice before worker's call\n");
        exit(1);
    }
    //Dynamic variables
    double t_sens=atof(argv[5]);
    double tau=(2*pi)*t_steps*t_sens;
    double h = (2*pi)*t_sens;  // Time step
    

    //Definitions for multithread
    pthread_t threads[NUM_THREADS];
    ThreadData thread_data[NUM_THREADS];  // Array of data structures
    pthread_barrier_init(&barrier, NULL, NUM_THREADS); // Initialize the barrier with the number of threads

    //Print of important info
    if (debug)
    {
        printf("\n***************************   Into C code   ***************************\n\n");
        printf("                               VARIABLES \n");
        printf("Cham X:                           %i \n", chamb_x);
        printf("Cham Y:                           %i \n", chamb_y);
        printf("Cham Z:                           %i \n", chamb_z);
        printf("Dimension:                        %i \n", D      );
        printf("Time sensitivity:                 %lf\n", t_sens );
        printf("Time steps:                       %i \n", t_steps);
        printf("Number of threds:                 %i \n\n", NUM_THREADS);
        for(i=0;i<count_forcing_x;i++){
        printf("Forcing amplitudes X %i:             %lf \n", i, f0x[i]);
        printf("Forcing frequencie X %i:             %lf \n", i, omega_insidex[i]);}
        printf("\n");
        for(i=0;i<count_forcing_y;i++){
        printf("Forcing amplitudes Y %i:             %lf \n", i, f0y[i]);
        printf("Forcing frequencie Y %i:             %lf \n",i,  omega_insidey[i]);}
        printf("\n");
        for(i=0;i<count_forcing_y;i++){
        printf("Forcing amplitudes Z %i:              %lf \n", i, f0z[i]);
        printf("Forcing frequencie Z %i:              %lf \n",i,  omega_insidez[i]);}
        printf("\n");
        printf("Output normalized time:           %i\n\n", output_time);
    }
    //Definition of forcing: defined on the bottom edge SO JUST PUT X AXIS SIZE
    if (debug) printf("                           MALLOC of FORCING \n");

    //*******x axes*******/
    force_pt_x = (double *) malloc(sizeof(double) * chamb_x);
    if(force_pt_x == NULL){
        printf("WARNING: Out of memory on malloc call for forces on x\n");
        exit(1);
    }
    check=h_shape_f(force_points,force_pt_x,omega_insidex, count_forcing_x,f0x,count_forcing_omega_x,0);//last is time
    if (debug) printf("h_shape f on x check: %i \n", check);

    //*******y axes*******/
    force_pt_y = (double *) malloc(sizeof(double) * chamb_x);//Here is write the term chambx
    if(force_pt_y == NULL){
        printf("WARNING: Out of memory on malloc call for forces on y\n");
        exit(1);
    }
    check=h_shape_f(force_points,force_pt_y,omega_insidey, count_forcing_y,f0y,count_forcing_omega_y,0);//last is time
    if (debug) printf("h_shape f on y check: %i \n", check);

    //*******z axes*******/
    force_pt_z = (double *) malloc(sizeof(double) * chamb_x);//Here is write the term chambx
    
    if(force_pt_z == NULL){
        printf("WARNING: Out of memory on malloc call for forces on z\n");
        exit(1);
    }
    check=h_shape_f(force_points,force_pt_z,omega_insidez, count_forcing_z,f0z,count_forcing_omega_z,0);//last is time
    if (debug) printf("h_shape f on z check: %i \n", check);

    if (debug){
        printf("force_pt_x of 0        at time 0: %lf \n",force_pt_x[0]);
        printf("force_pt_y of 0        at time 0: %lf \n",force_pt_y[0]);
        printf("force_pt_z of 0        at time 0: %lf \n",force_pt_z[0]);

        printf("force_pt_x at full len at time 0: %lf \n",force_pt_x[chamb_x-1]);
        printf("force_pt_y at full len at time 0: %lf \n",force_pt_y[chamb_x-1]);
        printf("force_pt_z at full len at time 0: %lf \n",force_pt_z[chamb_x-1]);
    }
    // --------------------------------------

    //######################################## OPEN THE FIRST LATTICE #####################################
    printf("\n                           Opening first lattice\n");
    sem_id = semget(1234567,0,1);//the first one is key, the other one is semaphore value
    sem_id=1;//Bypass the semaphore
    if (-1 == sem_id) 
        {
            sem_id = 0;
            printf("Getting a handle to the semaphore failed; errno is %d \n", errno);
        }
    else 
        {
            // get a handle to the shared memory
            shm_id = shmget(1234567, 16, 0);
            
            if (shm_id == -1) {
                shm_id = 0;
                printf("Couldn't get a handle to the shared memory; errno is %d \n", errno);
            }
            else {
                printf("Shared memory's id is %d \n", shm_id);

                // Attach the memory.
                address = shmat(shm_id, NULL, 0);

                if ((void *)-1 == address) {
                    address = NULL;
                    printf("Attaching the shared memory failed; errno is %d \n", errno);
                }
                else {
                    printf("Shared memory address = %p \n", address);
                    
                }
            }
        }
    //######################################## OPEN THE DUAL LATTICE #####################################
    printf("\n                           Opening second lattice\n");
    sem_id = semget(12345678,0,1);//the first one is key, the other one is semaphore value
    sem_id=1;//Bypass semaphore
    if (-1 == sem_id) 
        {
            sem_id = 0;
            printf("Getting a handle to the semaphore failed FOR DUAL LATTICE; errno is %d \n", errno);
        }
    else 
        {
            // get a handle to the shared memory
            shm_id = shmget(12345678, 16, 0);
            
            if (shm_id == -1) {
                shm_id = 0;
                printf("Couldn't get a handle to the shared memory FOR DUAL LATTICE; errno is %d \n", errno);
            }
            else {
                printf("Shared memory's id is %d FOR DUAL LATTICE\n", shm_id);

                // Attach the memory.
                dual_lattice_pt = shmat(shm_id, NULL, 0);

                if ((void *)-1 == dual_lattice_pt) {
                    dual_lattice_pt = NULL;
                    printf("Attaching the shared memory failed FOR DUAL LATTICE; errno is %d \n", errno);
                }
                else {
                    printf("Shared memory address FOR DUAL LATTICE = %p \n", dual_lattice_pt);
                }
            }
        }
    //##########################################    Assign the work    ##################################################
    //The job must get the correct memory locations. Run the RK through the correct time steps so that each
    //step is like the other jobs and at the end of each step the grid is updated.

    //----------------------------------------- Create threads and pass data to them -------------------------------------------
    //Strips are evaluated to divide the work
    int* extremes_vec_start=malloc(sizeof(int) * NUM_THREADS);
    int* extremes_vec_end=malloc(sizeof(int) * ((int) NUM_THREADS));
    for(i=0;i<NUM_THREADS;i++){//controlled initialization to zero and not to random values
        extremes_vec_start[i]=0;
        extremes_vec_end[i]=0;
    }
    int quoto,resto;
    extremes_vec_start[0]=1;
    if( NUM_THREADS==1) extremes_vec_end[0]=chamb_y-2+extremes_vec_start[0]; //Case of dimension 1, however, you do not work with less than 3*3. This extreme is not included: start<=work<stop
    if( NUM_THREADS!=1) {
        quoto=((chamb_y-2)/ NUM_THREADS);
        resto=((chamb_y-2)% NUM_THREADS);
        extremes_vec_start[0]=1;
        extremes_vec_end[0]=chamb_y-2+extremes_vec_start[0];
        for(i=0;i<NUM_THREADS-1;i++){
            extremes_vec_start[i+1]=extremes_vec_start[i]+quoto;
        }
        for(i=0;i<NUM_THREADS-1;i++){
            extremes_vec_end[i]=extremes_vec_start[i+1];
        }
        extremes_vec_end[NUM_THREADS-1]=chamb_y-2+extremes_vec_start[0];
    }


    if(debug){
        printf("\n                Position of strips sent to workers\n");
        for(i=0;i<NUM_THREADS;i++){
            printf("Worker %i, initial: %i, final: %i\n",i,extremes_vec_start[i],extremes_vec_end[i]);
        }
    }
    //----------------------------------------- Create the data for each worker -------------------------------------------
    for (int m = 0; m < NUM_THREADS; m++) {
        thread_data[m].id = m ;  // Each thread has a unique ID
        thread_data[m].tau = tau;  // Assign values to each thread
        thread_data[m].passo = h;
        thread_data[m].num_steps_ext = t_steps;
        thread_data[m].num_steps_internal = 1;
        thread_data[m].address=address;//general address because it can get its address with the id and the size
        thread_data[m].dual_lattice_pt=dual_lattice_pt;
        thread_data[m].cham_dim_x=chamb_x;
        thread_data[m].cham_dim_y=chamb_y;
        thread_data[m].n_thread_tot=NUM_THREADS;
        thread_data[m].start=extremes_vec_start[m];
        thread_data[m].stop=extremes_vec_end[m];
        thread_data[m].dimension=D;//The physical dimension of the system
        thread_data[m].kxp_position=kxp_pt;

        thread_data[m].forcing_x=force_pt_x;
        thread_data[m].forcing_y=force_pt_y;
        thread_data[m].forcing_z=force_pt_z;

        thread_data[m].omega_insidex=omega_insidex;
        thread_data[m].omega_insidey=omega_insidey;
        thread_data[m].omega_insidez=omega_insidez;

        thread_data[m].f0x=f0x;
        thread_data[m].f0y=f0y;
        thread_data[m].f0z=f0z;

        thread_data[m].len_omega_insidex=count_forcing_omega_x;
        thread_data[m].len_omega_insidey=count_forcing_omega_y;
        thread_data[m].len_omega_insidez=count_forcing_omega_z;

        thread_data[m].len_f0x=count_forcing_x;
        thread_data[m].len_f0y=count_forcing_y;
        thread_data[m].len_f0z=count_forcing_z;

        thread_data[m].output_time=output_time;
        thread_data[m].nome_file=argv[22];
        thread_data[m].thermal_coupling=atof(argv[15]);

        thread_data[m].start_plotX=extr_start_plot_x;
        thread_data[m].stop_plotX=extr_end_plot_x;
        thread_data[m].start_plotY=extr_start_plot_y;
        thread_data[m].stop_plotY=extr_end_plot_y;
        // printf("\n\n ***************extr_start_plot_x:%i**************** \n",extr_start_plot_x);
        // printf(" ***************extr_end_plot_x:%i**************** \n",extr_end_plot_x);
        // printf(" ***************extr_start_plot_y:%i**************** \n",extr_start_plot_y);
        // printf(" ***************extr_end_plot_y:%i**************** \n\n",extr_end_plot_y);

        if (pthread_create(&threads[m], NULL, assignment, (void *)&thread_data[m]) != 0) {
            perror("Error in thread creation");
            return 1;
        }
    }    
    //-------------------------------------  Wait for the workers and destroy them as soon as they are finished ---------------------------------
    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }
    // Destroy the barrier
    pthread_barrier_destroy(&barrier);
    free(kxp_pt);
    free(force_pt_x);
    free(force_pt_y);
    free(force_pt_z);
    free(extremes_vec_start);
    free(extremes_vec_end);
    return 0;
}
