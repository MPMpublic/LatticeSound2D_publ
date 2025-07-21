# LatticeSound 2D open-source
#
# Copyright (C) 2025  Giorgio Lo Presti (MPMpublic)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
import numpy as np
import sys
import glob
import os
import matplotlib.pyplot as plt

debug=False

chamX=int(sys.argv[1])
chamY=int(sys.argv[2])
chamZ=int(sys.argv[3])
D=int(sys.argv[4])
tot_step=int(sys.argv[5])
output_time=int(sys.argv[6])
plot=eval(sys.argv[7])
n_proc=int(sys.argv[8])
output_dir=sys.argv[9]
times=[]
for t in range(tot_step):
    if t%output_time==0:
        times.append(t)

# Pattern per cercare tutti i .txt nella cartella
pattern = os.path.join(output_dir, 'Trd*.txt')

# Trova tutti i file che corrispondono al pattern
file_txt_base = glob.glob(pattern)
file_txt=[]
for t in times: 
    for i in range(n_proc):
        for file in file_txt_base:
            if str(i)==file.split("Trd")[1].split("_t")[0]:
                if str(t)==file.split("_t")[1].split(".")[0]:
                    if debug: print(t,i,file)
                    file_txt.append(file)

displacements_time=[]
sigma_time=[]
speed_time=[]
sigma_speed_time=[]
forces_time=[]
points_time=[]
sigma_forces_time=[]
# Legge e stampa il contenuto di ciascun file
for t_step in times:
    displacements=[]
    sigmas=[]
    speeds=[]
    sigma_speeds=[]
    forces=[]
    sigma_forces=[]
    points=[]
    for i in range(n_proc):
        with open(os.path.join(output_dir, 'Trd')+str(i)+"_t"+str(t_step)+".txt", 'r', encoding='utf-8') as f:
            linee = f.readlines()
        displ=[]
        sigm=[]
        spd=[]
        sigma_spd=[]
        frc=[]
        sgm_frc=[]
        
        for number in linee[0].split(" "):
            displ.append(float(number))
        for number in linee[1].split(" "):
            sigm.append(float(number))
        for number in linee[2].split(" "):
            spd.append(float(number))
        for number in linee[3].split(" "):
            sigma_spd.append(float(number))
        for number in linee[4].split(" "):
            frc.append(float(number))
        for number in linee[5].split(" "):
            sgm_frc.append(float(number))
            
        for linea in linee[6:len(linee)]:
            pt=[]
            for number in linea.split(" "):
                pt.append(float(number))
            points.append(pt)
        
        displacements.append(displ)
        sigmas.append(sigm)
        speeds.append(spd)
        sigma_speeds.append(sigma_spd)
        forces.append(frc)
        sigma_forces.append(sgm_frc)

    displacements_time.append(displacements)
    sigma_time.append(sigmas)
    speed_time.append(speeds)
    sigma_speed_time.append(sigma_speeds)
    forces_time.append(forces)
    sigma_forces_time.append(sigma_forces)
    points_time.append(points)

if debug: 
    print(displacements_time)
    print(sigma_time)
    print(speed_time)
    print(sigma_speed_time)
    print(forces_time)
    print(sigma_forces_time)

displacements=[]
sigma=[]
speed=[]
sigma_speed=[]
forces=[]
sigma_forces=[]
for i in range(len(displacements_time)):
    displacements.append([np.mean(np.array(displacements_time[i])[:, 0]),np.mean(np.array(displacements_time[i])[:, 1])])
    sigma.append([np.mean(np.array(sigma_time[i])[:, 0]),np.mean(np.array(sigma_time[i])[:, 1])])
    speed.append([np.mean(np.array(speed_time[i])[:, 0]),np.mean(np.array(speed_time[i])[:, 1])])
    sigma_speed.append([np.mean(np.array(sigma_speed_time[i])[:, 0]),np.mean(np.array(sigma_speed_time[i])[:, 1])])
    forces.append([np.mean(np.array(forces_time[i])[:, 0]),np.mean(np.array(forces_time[i])[:, 1])])
    sigma_forces.append([np.mean(np.array(sigma_forces_time[i])[:, 0]),np.mean(np.array(sigma_forces_time[i])[:, 1])])
    
if plot:
    print("\ndisplacements\n")
    print(displacements)
    print("\nsigma")
    print(sigma)
    print("\nspeed")
    print(speed)
    print("\nsigma speed")
    print(sigma_speed)
    print("\nforces")
    print(forces)
    print("\nsigma forces")
    print(sigma_forces)
with open(os.path.join(output_dir,"LS_result.txt"),"a") as file:
    file.write("\nDisplacements\n")
    for i in range(len(displacements_time)):
        file.write(str(displacements[i])+"\n")
    file.write("\nSigma Displacement\n")
    for i in range(len(displacements_time)):
        file.write(str(sigma[i])+"\n")
    file.write("\nSpeed\n")
    for i in range(len(displacements_time)):
        file.write(str(speed[i])+"\n")
    file.write("\nSigma Speed\n")
    for i in range(len(displacements_time)):
        file.write(str(sigma_speed[i])+"\n")
    file.write("\nForces\n")
    for i in range(len(displacements_time)):
        file.write(str(forces[i])+"\n")
    file.write("\nSigma Forces\n")
    for i in range(len(displacements_time)):
        file.write(str(sigma_forces[i])+"\n")
print()
if debug: print("first element read: ",points_time[0][0])

if plot:
    for t,points in enumerate(points_time):
        mom_points=np.array(points)
        # print(mom_points)
        speeds = np.sqrt(mom_points[:,4]**2 + mom_points[:,5]**2)
        forces = np.sqrt(mom_points[:,6]**2 + mom_points[:,7]**2)
        # print(speeds)
        max_speed=np.max(speeds)
        max_force=np.max(forces)
        plt.figure(figsize=(8, 6))
        for i, point in enumerate(points):
            plt.scatter(point[2], point[3],c=forces[i],vmin=0, vmax=max_force,cmap='viridis',s=40)#,edgecolors="blue"
            plt.quiver(point[2], point[3], point[4], point[5], speeds[i], angles='xy', scale_units='xy', scale=max_speed*5, cmap='coolwarm')
        plt.title(f"Displacement plot with vector velocities and forces intensity at time plot {t}")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.grid(True)
        # plt.axis('equal')
        plt.show()

#OLD DEPRECATED
# for t_step in range(tot_step):
#     displacements=[]
#     sigmas=[]
#     speeds=[]
#     sigma_speeds=[]
#     forces=[]
#     sigma_forces=[]
#     points=[]
#     for file in file_txt:
#         if str(t_step)==file.split("_t")[1].split(".")[0]:
#             print(file)
#             with open(file, 'r', encoding='utf-8') as f:
#                 linee = f.readlines()
#             displ=[]
#             sigm=[]
#             spd=[]
#             sigma_spd=[]
#             frc=[]
#             sgm_frc=[]
            
#             for number in linee[0].split(" "):
#                 displ.append(float(number))
#             for number in linee[1].split(" "):
#                 sigm.append(float(number))
#             for number in linee[2].split(" "):
#                 spd.append(float(number))
#             for number in linee[3].split(" "):
#                 sigma_spd.append(float(number))
#             for number in linee[4].split(" "):
#                 frc.append(float(number))
#             for number in linee[5].split(" "):
#                 sgm_frc.append(float(number))
                
#             for linea in linee[6:len(linee)]:
#                 pt=[]
#                 for number in linea.split(" "):
#                     pt.append(float(number))
#                 points.append(pt)
            
#             displacements.append(displ)
#             sigmas.append(sigm)
#             speeds.append(spd)
#             sigma_speeds.append(sigma_spd)
#             forces.append(frc)
#             sigma_forces.append(sgm_frc)
#             print(displ)
            
#     displacements_time.append(displacements)
#     sigma_time.append(sigmas)
#     speed_time.append(speeds)
#     sigma_speed_time.append(sigma_speeds)
#     forces_time.append(forces)
#     sigma_forces_time.append(sigma_forces)
#     points_time.append(points)
#             # displacements_time.append(displ)
#             # sigma_time.append(sigm)
#             # speed_time.append(spd)
#             # sigma_speed_time.append(sigma_spd)
#             # forces_time.append(frc)
#             # sigma_forces_time.append(sgm_frc)
#             # points_time.append(points)