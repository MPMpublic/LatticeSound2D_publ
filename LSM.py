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
"""This script is Lattice Sound Manager. Its role is to initialize the C_Lattice, giving to it the system dimensions and the initial 
   conditions and preparing the shared memory slots for the dynamic c codes (through sysv_ipc). It launches also Py_Analyses to analyse
   the system giving it the shared locations. The info definitions are taken from the input Starter.flsm and the matrices are built using
   the library LSM_lib.
   Before launching the C code, it must be compiled and before compiling, the memory must be already allocated.
   The oscillator is defined in this way: [xi,yi,zi,vxi,vyi,vzi]. The dual lattice is of the type: [mu,km,gamm,1]. All of these quantities
   are referred to the quantities of the "base element" (the one nomined in build_lattice).
   MEMORY ALLOCATION: it prepares location and semaphores (semaphores not used because the real multithread is operated by C_Lattice
   and every worker must be able to enter in every position and moment but is C_Lattice to decide what to do). However the number of 
   threads must be specified here. The memory is dimension of the system * 8 (double) * 2 (initial position and speed). LSM must furnish 
   the initial lattice (with the position and speed of every oscillator). At the end, it will close the shared memory.
   BE CAREFUL: READ THE NORMALIZATION IN INITIALIZATION PART (SPACE and SPEED OF NORMALIZED TO 1). Remember that shared memory
   works with bytes number."""
########################################################   IMPORT AND NAMES  ########################################################
import subprocess
import os
import sysv_ipc # type: ignore for suppress warning
from multiprocessing import Process
import numpy as np
from LSM_Lib import build_lattice, insert_object_2D, TwoD_object, rectangle, ellipse, leng, elastic_const, mass, scale_back2D

C_name="C_Lattice"
Py_name="Py_Analyses"

########################################################   SYSTEM DEFINITION  ########################################################
"""It is read by the Starter.flsm.
   It must contain the GEOMETRY: number of oscillator per each axes of the chamber (which must be coherent with the dimension D). 
   reticular_distance is the space [m] between an oscillator and the other on the same axes; OSCILLATOR INFO: elastic_const, mass and 
   evaluates tau; TIME SYM INFO (sensitivity and step).
"""
#######################################################
#Reading from starter file flsm system info
with open("Starter.flsm","r") as file:
    lfsm=file.readlines()
file.close()
strings_for_commands=[]
for i,line in enumerate(lfsm):
    if "debug" in line:
        debug=eval(line.split("=")[1][0:len(line.split("=")[1])-1])
    if "thread_num" in line:
        thread_num=int(line.split("=")[1][0:len(line.split("=")[1])-1])
    if "chamX" in line:
        chamX=int(line.split("=")[1][0:len(line.split("=")[1])-1])
    if "chamY" in line:
        chamY=int(line.split("=")[1][0:len(line.split("=")[1])-1])
    if "chamZ" in line:
        chamZ=int(line.split("=")[1][0:len(line.split("=")[1])-1])
    if "D" in line and "=" in line:
        D=int(line.split("=")[1][0:len(line.split("=")[1])-1])
    if "Temperature" in line:
        temperature=float(line.split("=")[1][0:len(line.split("=")[1])-1])
    if "scaling" in line:
        scaling=float(line.split("=")[1][0:len(line.split("=")[1])-1])
    if "thermal_coupling" in line:
        thermal_coupling=float(line.split("=")[1][0:len(line.split("=")[1])-1])
    if "temporal_sensitivity" in line:
        temporal_sensitivity_numerical=float(line.split("=")[1][0:len(line.split("=")[1])-1])
    if "temporal_steps" in line:
        temporal_steps=int(line.split("=")[1][0:len(line.split("=")[1])-1])
        output_time=temporal_steps
    if "material" in line:
        substance=str(line.split("=")[1][0:len(line.split("=")[1])-1])
    if "state" in line:
        state=line.split("=")[1][0:len(line.split("=")[1])-1]
    if "FORCING" in line:
        forcing_amplitude_x=eval(lfsm[i+1])
        forcing_amplitude_str_x=lfsm[i+1][0:len(lfsm[i+1])-1]
        forcing_omega_x=eval(lfsm[i+2])
        forcing_omega_str_x=lfsm[i+2][0:len(lfsm[i+2])-1]
        
        forcing_amplitude_y=eval(lfsm[i+4])
        forcing_amplitude_str_y=lfsm[i+4][0:len(lfsm[i+4])-1]
        forcing_omega_y=eval(lfsm[i+5])
        forcing_omega_str_y=lfsm[i+5][0:len(lfsm[i+5])-1]
        
        forcing_amplitude_z=eval(lfsm[i+7])
        forcing_amplitude_str_z=lfsm[i+7][0:len(lfsm[i+7])-1]
        forcing_omega_z=eval(lfsm[i+8])
        forcing_omega_str_z=lfsm[i+8][0:len(lfsm[i+8])-1]
        
    if "BUILDINGS" in line:
        n_build=int(lfsm[i+1].split("=")[1][0:len(lfsm[i+1].split("=")[1])-1])
        for j in range(n_build):
            strings_for_commands.append(lfsm[i+2+j][0:len(lfsm[i+2+j])-1])
            
    if "OUTPUT" in line:
        output_time=int(lfsm[i+1].split("=")[1][0:len(lfsm[i+1].split("=")[1])-1])
        output_dir=lfsm[i+2][0:len(lfsm[i+2])-1]
        plot=lfsm[i+3].split("=")[1][0:len(lfsm[i+3].split("=")[1])-1]
        ChamX_plot=eval(lfsm[i+4].split("=")[1][0:len(lfsm[i+4].split("=")[1])-1])
        ChamY_plot=eval(lfsm[i+5].split("=")[1][0:len(lfsm[i+5].split("=")[1])-1])
        ChamZ_plot=eval(lfsm[i+6].split("=")[1][0:len(lfsm[i+6].split("=")[1])-1])

if((int(ChamX_plot[1])>int(chamX)) or (int(ChamY_plot[1])>int(chamY)) or (int(ChamZ_plot[1])>int(chamZ))):
    print("-------------------  WARNING: PLOT DIMENSION OUTSIDE THE LATTICE  -------------------")

#######################################################
#Geometry: number of oscillators per each dimension D 
# and reticular_distance
                                               # dimension
reticular_distance=leng[substance]/scaling                               # [m]
if debug:
    x_len=chamX*reticular_distance
    y_len=chamY*reticular_distance
    z_len=chamZ*reticular_distance
    print("Body dimensions:",x_len,y_len,z_len)


#####################################################  DEFINITION PROGRAM COMPILING  #####################################################
"""It contains the info to compile the C code, and to lauch both C and Python. To launch C, it needs: chambers, dimension, zero elastic 
   const, temporal info and reticular distance. BUT NOW IS NOT COMPILED.
"""
########################################################
#################################### C program Compiling
command_compile="cd "+str(os.getcwd())+"&& gcc "+C_name+".c h_shape.c C_Lattice_Worker.c C_Lattice_dif_sys.c -pthread -lgsl -lgslcblas -lm  -o "+C_name
subprocess.run(command_compile,shell=True) # gcc -o oscillatore oscillatore.c -lgsl -lgslcblas -lm
# -g -fsanitize=address
########################################################
################################# Define codes to launch
def Py_code(a):
    command_pylaunch="cd "+str(os.getcwd())+"&& python3 "+Py_name+".py "+str(chamX)+" "+str(chamY)+" "+str(chamZ)+" "+str(D)+" "+str(temporal_steps)+" "+str(output_time)+" "+str(plot)+" "+str(thread_num)+" "+output_dir
    subprocess.run(command_pylaunch,shell=True)
    return 0
def C_code(b):
    command_exe="cd "+str(os.getcwd())+"&& ./"+C_name+" "+str(chamX)+" "+str(chamY)+" "+str(chamZ)+" "+str(D)+" "+str(temporal_sensitivity_numerical)+" "+str(temporal_steps)\
        +" "+str(thread_num)+" "+forcing_amplitude_str_x+" "+forcing_omega_str_x+" "+forcing_amplitude_str_y+" "+forcing_omega_str_y\
        +" "+forcing_amplitude_str_z+" "+forcing_omega_str_z+" "+str(output_time)+" "+str(thermal_coupling)+" "+str(ChamX_plot[0])\
        +" "+str(ChamX_plot[1])+" "+str(ChamY_plot[0])+" "+str(ChamY_plot[1])+" "+str(ChamZ_plot[0])+" "+str(ChamZ_plot[1])+" "+output_dir#os.path.join(os.getcwd(),)
    subprocess.run(command_exe,shell=True)
    if debug: print("C Command launched \n",command_exe)
    return 0

########################################################    MEMORY ALLOCATION   ########################################################
"""Just on big memory is created because all the c worker should reach every point of this memory. The multithread is operated by 
   C_lattice. Here the "try" is to cancel the previous memory in case of anuspected shout down. 1234567 is the first key.
   Note: Size is in byte. In caso come size: size=sysv_ipc.PAGE_SIZE (che è la page size del sistema operativo da cercare su internet 
   oppure suoi multipli). SizeMem=8*SisDym*D*2 because I work with double per each dimension and it needs the speed.
   In lattice_positions and lattice_info 0 or 1 stand for the kind of lattice"""
########################################################
########################## Create the lattice by library
if n_build!=0:
    if D==2:
        lattice=build_lattice([chamX,chamY],substance,state,temperature, scaling)
    if D==3:
        lattice=build_lattice([chamX,chamY,chamZ],substance,state,temperature, scaling)

    for j in range(n_build):
        if ("lattice" in strings_for_commands[j]) and (not("duallattice" in strings_for_commands[j])):
            alfa=int(strings_for_commands[j][8])
            beta=int(strings_for_commands[j][11])
            latticeres=lattice
            latticeres[0][alfa][beta]=eval(strings_for_commands[j].split("=")[1])
            print(latticeres[0][alfa][beta])
        elif "duallattice" in strings_for_commands[j]:
            alfa=int(strings_for_commands[j][12])
            beta=int(strings_for_commands[j][15])
            latticeres=lattice
            latticeres[1][alfa][beta]=eval(strings_for_commands[j].split("=")[1])
            print(latticeres[1][alfa][beta])
        else:
            building=eval(strings_for_commands[j])
            if D==2:
                latticeres=insert_object_2D(lattice,building,temperature,substance,scaling)#Here is temperature of base material (inside: basetemperature)
else:
    if D==2:
        latticeres=build_lattice([chamX,chamY],substance,state,temperature, scaling)
    if D==3:
        latticeres=build_lattice([chamX,chamY,chamZ],substance,state,temperature, scaling)

if debug:
    print("Real lattice")
    for l in latticeres[0]:
        print(l)
    print("Dual lattice")
    for l in latticeres[1]:
        print(l)

lattice_positions=np.array(latticeres[0],dtype=np.float64)
lattice_info=np.array(latticeres[1],dtype=np.float64)

scaled_len, scaled_mass, scaled_k, scaled_gamm, scaled_omega=scale_back2D(substance,state,temperature,scaling)
with open(os.path.join(output_dir,"LS_result.txt"),"w") as file:
    file.write(f"Effective length          [m]:        {scaled_len}\n")
    file.write(f"Effective mass            [kg]:       {scaled_mass}\n")
    file.write(f"Effective elastic const   [kg/s^2]:   {scaled_k}\n")
    file.write(f"Effective viscosity const [1/s]:      {scaled_gamm}\n")
    file.write(f"Effective omega0          [rad/s]:    {scaled_omega}\n")
    file.write(f"thermal_coupling          [ ]:        {thermal_coupling}\n")

########################################################
############################### Create the shared memory
SisDym=chamX*chamY*chamZ
SizeMem=8*SisDym*D*2

try:
    memory = sysv_ipc.SharedMemory(1234567, sysv_ipc.IPC_CREX,size=SizeMem)
    semaphore = sysv_ipc.Semaphore(1234567, sysv_ipc.IPC_CREX)
except:
    memory = sysv_ipc.SharedMemory(1234567, sysv_ipc.IPC_CREAT,size=0)
    semaphore = sysv_ipc.Semaphore(1234567, sysv_ipc.IPC_CREAT)
    sysv_ipc.remove_shared_memory(memory.id)
    sysv_ipc.remove_semaphore(semaphore.id)
    memory = sysv_ipc.SharedMemory(1234567, sysv_ipc.IPC_CREAT,size=SizeMem)
    semaphore = sysv_ipc.Semaphore(1234567, sysv_ipc.IPC_CREAT)


try:
    memory2 = sysv_ipc.SharedMemory(12345678, sysv_ipc.IPC_CREX,size=SizeMem)
    semaphore2 = sysv_ipc.Semaphore(12345678, sysv_ipc.IPC_CREX)
except:
    memory2 = sysv_ipc.SharedMemory(12345678, sysv_ipc.IPC_CREAT,size=0)
    semaphore2 = sysv_ipc.Semaphore(12345678, sysv_ipc.IPC_CREAT)
    sysv_ipc.remove_shared_memory(memory2.id)
    sysv_ipc.remove_semaphore(semaphore2.id)
    memory2 = sysv_ipc.SharedMemory(12345678, sysv_ipc.IPC_CREAT,size=SizeMem)
    semaphore2 = sysv_ipc.Semaphore(12345678, sysv_ipc.IPC_CREAT)

########################################################
################ In Case you need to print in debug mode
if debug:
    print("Prepared lattice")
    for i in range(chamX):
        print(lattice_positions[i][:])

for a in range(chamX):
    memory.write(lattice_positions[a].tobytes('C'),a*chamY*8*4)#here the second argument is that of the starting point,obj_list=[]#1*chamY*8*4
if debug: 
    print("Sent lattice")
    print(np.frombuffer(memory.read(), dtype=np.float64))

if debug:
    print("Prepared info")
    for i in range(chamX):
        print(lattice_info[i][:])

for a in range(chamX):
    memory2.write(lattice_info[a].tobytes('C'),a*chamY*8*4)#here the second argument is that of the starting point,obj_list=[]#1*chamY*8*4
if debug: 
    print("Sent info")
    print(np.frombuffer(memory2.read(), dtype=np.float64))

########################################################    CODE LAUNCH  ########################################################
"""Uses process to launch the sub-processes. Remember to close the shared memory."""
########################################################
################################## Launch of codes
# creating processes
processes=[]
p1 = Process(target=C_code, args=[0])
processes.append(p1)
p1.start()

# completing process
for p in processes:
    p.join()

processes=[]
# print("RICORDA DI RIACCENDERE IL CODICE DI ANALISI!!!")
p2 = Process(target=Py_code, args=[0])
processes.append(p2)
p2.start()
# completing process
for p in processes:
    p.join()

# print(memory.read())
########################################################
################### Closing shared memory FUNDAMENTAL!!!
sysv_ipc.remove_shared_memory(memory.id)
sysv_ipc.remove_semaphore(semaphore.id)

sysv_ipc.remove_shared_memory(memory2.id)
sysv_ipc.remove_semaphore(semaphore2.id)

######################################################    IMPORTANT WEBSITES  ######################################################
#https://www.ibm.com/docs/es/aix/7.2?topic=memory-creating-shared-segment-shmat-subroutine
#https://semanchuk.com/philip/sysv_ipc/
#https://stackoverflow.com/questions/12749176/how-to-access-variables-in-shared-memory
#https://github.com/martinohanlon/c_python_ipc
#https://devblogs.microsoft.com/oldnewthing/20210510-00/?p=105200


######################################################    Old pieces of code     #################################################

#s=bytes(number, encoding='utf-8')
#number.to_bytes(4, byteorder='big',signed=False)#"Proviamoci"

#s.encode()) #encode perché prende le str come binario

# # Create shared memory object
# memory = sysv_ipc.SharedMemory(123456)
# # Read value from shared memory
# memory_value = memory.read()
# # Find the 'end' of the string and strip
# i = memory_value.find('\0')
# if i != -1:
#     memory_value = memory_value[:i]
# print (memory_value)

# print(np.zeros((chamX,chamY,chamZ), dtype=float, order='C').tobytes('C'))
# memory.write(np.zeros((chamX,chamY,chamZ), dtype=float, order='C').tobytes('C'))
# A=np.array([[1,2],[3,4]])
# print(A[0,1])

# print("verifica solutore: ",np.cos(2*2.5132741229e-03))
# print("verifica solutore: ",1.169046*np.exp(-0.382022*2*2.5132741229e-03)-0.169046*np.exp(-2.681022*2*2.5132741229e-03))


#OLD LATTICE DEFINITION
# for j in range(chamY):
#     prearray=np.array([[float(j), float(i), 0., 0.] for i in range(chamX)], dtype=np.float64)
#     # print(j*chamX*8*4)
#     # print(prearray.tobytes('C'))
#     memory.write(prearray.tobytes('C'),j*chamX*8*4)#,

# memory.write((np.array([[0.,0.,0.,0.],[0.,1.,0.,0.],[0.,2.,0.,0.],[1.,0.,0.,0.],[1.2,1.3,0.,0.],[1.,2.,0.,0.], 
#                        [2.,0.,0.,0.],[2.,1.,0.,0.],[2.,2.,0.,0.]])).tobytes('C'))
#                        #Positions ([xi,yi,vxi,vyi] per la lungh y) per lungh x e la y cresce da zero
#                        #DAGLIELI NORMALIZZATI!
# memory.write((np.array([[0.,0.,0.,0.],[0.,1.,0.,0.],[0.,2.,0.,0.],[0.,3.,0.,0.],[1.,0.,0.,0.],[1.3,1.5,0.,0.],[1.,2.,0.,0.],[1.,3.,0.,0.], 
#                        [2.,0.,0.,0.],[2.,1.,0.,0.],[2.,2.,0.,0.],[2.,3.,0.,0.],[3.,0.,0.,0.],[3.,1.,0.,0.],[3.,2.,0.,0.],[3.,3.,0.,0.]])).tobytes('C'))
#                        #Positions ([xi,yi,vxi,vyi] per la lungh y) per lungh x e la y cresce da zero
#                        #DAGLIELI NORMALIZZATI!
# print(np.frombuffer(memory.read(), dtype=np.float64))
