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
double h_shape_g(double, double, double);
int h_shape_f(int,double*,double*, int, double *, int, double);

// Struttura che contiene i dati da passare al thread
typedef struct {
    int id;  // ID del thread
    double tau;  // Un valore numerico da passare al thread
    double passo;
    int num_steps_ext;
    int num_steps_internal;
    unsigned char* address;
    unsigned char* dual_lattice_pt;
    int cham_dim_x;
    int cham_dim_y;
    int n_thread_tot;
    int start;
    int stop;
    int dimension;
    double* kxp_position; // Questo è il vettore di posizioni e velocità mandate al worker
    double* forcing_x;
    double* forcing_y;
    double* forcing_z;

    double* omega_insidex;
    double* omega_insidey;
    double* omega_insidez;

    double* f0x;
    double* f0y;
    double* f0z;

    int len_omega_insidex;
    int len_omega_insidey;
    int len_omega_insidez;

    int len_f0x;
    int len_f0y;
    int len_f0z;

    int output_time;
    double thermal_coupling;
    char* nome_file; //Questo va messo all'ultimo perché, quando traduce le stringhe, scambia lo spazio per una stringa

    int start_plotX;
    int stop_plotX;
    int start_plotY;
    int stop_plotY;

} ThreadData;

typedef struct {
    double* sn;
    double* dx;
    double* sopra;
    double* sotto;
    double* v_sn;
    double* v_dx;
    double* v_sopra;
    double* v_sotto;
    int* indexes;
    double* gamma_zero;
    double* km_pt;
    double* mu;
    double transverse;
    double* forze;
    double* u_rand_vect; //Randomic speed caused by photons

} InternalData;

typedef struct {
    double pos_x;
    double pos_y;
    double speed_x;
    double speed_y;
} ParticleData;
//Old codes
// in ThreadData:
// int h_shape_k(int,int,int,double*,double);//3 dim and a pt*
    // int address_x;
    // int address_y;
    // double* km_pt;
    // double* kp_position;