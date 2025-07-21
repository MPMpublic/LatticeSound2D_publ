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
"""This is the LSM library. The following functions are found: norm_Max_Boltz, propose_in_position, build_lattice, insert_object_2D. 
   Also, the definitions and methods on the following classes are found: TwoD_object, rectangle(TwoD_object), ellipse(TwoD_object).
   The task of this library is to build the lattice and the dual lattice based on the information that LSM has read from the 
   initialization file."""
######################################################    IMPORT & CONST   ######################################################
import numpy as np
import time
import matplotlib.pyplot as plt

debug=False

k_b=1.380649e-23                                                                       # J K^(-1)
Na=6.023e23
h_bar=1.054571817e-34                                                                   # J s
sigma_b=5.670374419e-8                                                                  # W/(m^2 K^4)
e_charge=1.602176634e-19                                                                # C
J_to_eV=1/e_charge                                                                      # eV/J
m_H2O=18.015/Na/1000                                                                    # kg
m_air=(2*14.0067*0.78084+2*15.9994*0.209476+2*39.948*0.00934+44.009*0.00041393)/Na/1000 # kg
m_p=1.67262313e-27                                                                       # kg
m_e=9.1093837015e-31                                                                    # kg
eps=8.8541878176e-12                                                                    # C^2/(N m^2)

#------------------------------------------   BUILDING FUNDAMENTAL BY GLOBAL QUANTITIES   ------------------------------------------
"""While Young's modulus is scale invariant and mass scaling is trivial because it preserves density, the definition of gamma is deduced 
   so that the molecule loses energy based on the dynamics induced by molecular bonds. The viscous friction must therefore be of the order
   of the force exerted by the particles."""
Tgl_air=78.80                                                                           # [K] T gas liquid air
omega_air_gl=k_b*Tgl_air/h_bar                                                          # [rad/s]
k_air_gl=(0.5)*m_air*omega_air_gl**2                                                    # [kg/(s^2)]
c_air=1005                                                                              # [J/(kg k)]
n_air=2.460977e25                                                                       # [at/m^3]

Tsl_H2O=273.15                                                                          # [K] T gas liquid air
omega_H2O_sl=k_b*Tsl_H2O/h_bar/3                                                        # [rad/s]; divided by 3 for Golstone bosons
k_H2O_sl=(0.5)*m_H2O*omega_H2O_sl**2                                                    # [kg/(s^2)]
c_ice=2090                                                                              # [J/(kg k)]

H2O_dens_0degree=997                                                                    # [kg/m^3]
H2Olen=2.75e-10                                                                         # [m] at solidification temperature https://www.its.caltech.edu/~atomic/snowcrystals/ice/ice.htm
air_dens_0degree=1.2922                                                                 # [kg/m^3]
airlen=(m_air/air_dens_0degree)**(1/3)                                                  # [m] average length of an effective air-air bond

deltaE_air=c_air*Tgl_air*air_dens_0degree                                               # [J] very approximate energy necessary to freeze 1 m^3 of air 
P_diss_rad_air=(airlen**3)*6*sigma_b*273.15**4                                          # [W] Power dissipated by the volume of an "air molecule" by radiation P=(V0/V)*S*I with unit V (so 6 faces of 1m^2). It is multiplied by the volume of the molecule because it is then divided
freez_time=deltaE_air/P_diss_rad_air*(airlen**3)                                        # [s] time to freeze 1 m^3 of air by radiation in void
E_osc_air=(1/2)*m_air*((m_e/m_p)*(airlen/2)**2)*(omega_air_gl**2)                       # [J] energy of an "air molecule" for oscillation
average_speed_air_gl=((m_e/m_p)**0.5)*(airlen/2)*omega_air_gl                           # [m/s]
gamma0_air=-np.log(P_diss_rad_air/(E_osc_air/freez_time))*(k_air_gl*((m_e/m_p)**0.5)*   # [1/s]
           (airlen/2)/average_speed_air_gl)*(1/np.e)                                    
# gamma0_H2O=1/(freez_time-(P_diss_rad_air/E_osc_air))*(1/omega_air_gl*freez_time)/3 Questa era dell'ordine giusto ma sbagliata

deltaE_H2O=c_ice*Tsl_H2O*H2O_dens_0degree                                               # [J] very approximate energy necessary to freeze 1 m^3 of water 
P_diss_rad_H2O=(H2Olen**3)*6*sigma_b*273.15**4                                          # [W] Power dissipated by the volume of an water molecule by radiation P=(V0/V)*S*I with unit V (so 6 faces of 1m^2).
freez_time_H2O=deltaE_H2O/P_diss_rad_H2O*(H2Olen**3)                                    # [s] time to freeze 1 m^3 of water by radiation in void
E_osc_H2O=(1/2)*m_H2O*((m_e/m_p)*(H2Olen/2)**2)*(omega_H2O_sl**2)+(1/(2*np.pi*eps))*(e_charge**2/(H2Olen/2))*3 # J energia di una molecola per oscillazione times Goldstone (solid?) + l'energia del legame solido
average_speed_H2O_sl=((m_e/m_p)**0.5)*(H2Olen/2)*omega_H2O_sl                           # m/s
# gamma0_H2O=-np.log(P_diss_rad_H2O/(E_osc_H2O/freez_time_H2O))*(k_H2O_sl*((m_e/m_p)**0.5)*(H2Olen/2)/average_speed_H2O_sl)/3*(1/np.e) #1/s (dissipa su 3 dimensioni)
gamma0_H2O=1/(freez_time_H2O-(P_diss_rad_H2O/E_osc_H2O))*(1/omega_H2O_sl*freez_time_H2O)/3

leng={"H2_O":H2Olen, "air":airlen}
mass={"H2_O":m_H2O, "air":m_air}
therm_dil={"H2_O":2.1e-4, "air":3.66e-3} #k^-1
elastic_const={"H2_O":k_H2O_sl, "air":k_air_gl}
gamma0={"H2_O":gamma0_H2O, "air":gamma0_air}

###################################################    GAUSS DISTRIBUTION FOR SPEED   ###################################################
"""These functions return Gaussian distributions which, when appropriately multiplied, provide Maxwellian distributions of velocities 
   randomly on multiple axes in [m/s]"""

def norm_Max_Boltz_1D(npart,T,element):
    sigma=np.sqrt(k_b * T / mass[element])
    vx = np.random.normal(0, sigma, npart)
    return vx                             #[m/s]

def norm_Max_Boltz_2D(npart,T,element):
    sigma=np.sqrt(k_b * T / mass[element])
    vx = np.random.normal(0, sigma, npart)
    vy = np.random.normal(0, sigma, npart)
    return np.stack((vx,vy), axis=1)      #[m/s]

def norm_Max_Boltz_3D(npart,T,element):
    sigma=np.sqrt(k_b * T / mass[element])
    vx = np.random.normal(0, sigma, npart)
    vy = np.random.normal(0, sigma, npart)
    vz =  np.random.normal(0, sigma, npart)
    return np.stack((vx,vy,vz), axis=1)   #[m/s]

###################################################    FUNCTION FOR DIRECT LATTICE   ###################################################
"""This function returns a set of 4-element (if in 2D) or 6-element (if in 3D) vectors corresponding to the 2 (or 3) spatial positions and
   the velocities obtained by normalizing with respect to the omega (already scaled) the velocity obtained by imposing an effective 
   temperature (dependent on the scaling)."""
def propose_in_position(initial_positions, final_positions, material, aggregation, temperature, scaling):#temperature in kelvin
    vectors=[]
    if (isinstance(temperature,int) or isinstance(temperature,float)):
        if len(initial_positions)==3:
            n_particles=(final_positions[0]-initial_positions[0])*(final_positions[1]-initial_positions[1])
            if debug:
                v0=np.power(3*k_b*temperature/mass[material],0.5)#m/s
                print("mass [kg]: ",mass[material])
                print("v0 [m/s] is: ",v0)
                print("Boltz speed: ",norm_Max_Boltz_3D(n_particles,temperature,material)*(10*v0)/n_particles)
                # import scipy
                # print("Integral of Max Boltz distribution",scipy.integrate.quad(norm_Max_Boltz_3D, 0, np.inf, args=(temperature,material), full_output=0, epsabs=1.49e-08, epsrel=1.49e-08))
                print("Number of particles", n_particles)
            
            if debug: t0=time.time()
            velocità=np.array([])
            if n_particles>100_000_000:
                for i in range(int(n_particles/100000000)+int(n_particles%100000000)):
                    if i==0:
                        velocità=norm_Max_Boltz_3D(100000000,temperature,material)[:]
                    else:
                        velocità=np.concatenate((velocità,norm_Max_Boltz_3D(100000000,temperature,material)[:]), axis=0)
                velocità=velocità[0:n_particles]
            else:
                velocità=norm_Max_Boltz_3D(n_particles,temperature,material)
            velocità=np.array(velocità)/((leng[material]*(1+therm_dil[material]*(temperature-273.15))*scaling*np.power(elastic_const[material]/(mass[material]*scaling**3),0.5))*scaling**3)
            if debug:

                print(velocità)
                t1=time.time()
                v_sq=[]
                print("tempo produzione",t1-t0)
                v_sq=np.power(velocità[:, 0]**2+velocità[:, 1]**2+velocità[:, 2]**2,0.5)
                t2=time.time()
                print("tempo totale",t2-t0)
                plt.hist(v_sq, bins=100, density=True, alpha=0.7, label='3D Maxwell')
                plt.xlabel('|v|')
                plt.ylabel('Distribuzione')
                plt.title('Distribuzione Maxwell-Boltzmann 3D')
                plt.legend()
                plt.grid(True)
                plt.show()
            point=0
            for i in range(initial_positions[0], final_positions[0],1):
                for j in range(initial_positions[1], final_positions[1],1):
                    for k in range(initial_positions[2], final_positions[2],1):
                        if((i==initial_positions[0]) or (j==initial_positions[1]) or (k==initial_positions[2])):
                            vectors.append([i,j,k,0.,0.,0.])
                        elif((i==final_positions[0]-1) or (j==final_positions[1]-1) or (k==final_positions[2]-1)):
                            vectors.append([i,j,k,0.,0.,0.])
                        else:
                            vectors.append([i,j,k,velocità[point][0],velocità[point][1],velocità[point][2]])
                        point+=1
        
        elif len(initial_positions)==2:
            n_particles=(final_positions[0]-initial_positions[0])*(final_positions[1]-initial_positions[1])
            if debug:
                v0=np.power(3*k_b*temperature/mass[material],0.5)#m/s
                print("mass [kg]: ",mass[material])
                print("v0 [m/s] is: ",v0)
                print("Boltz speed: ",norm_Max_Boltz_2D(n_particles,temperature,material)*(10*v0)/n_particles)
                # import scipy
                # print("Integral of Max Boltz distribution",scipy.integrate.quad(norm_Max_Boltz_3D, 0, np.inf, args=(temperature,material), full_output=0, epsabs=1.49e-08, epsrel=1.49e-08))
                print("Number of particles", n_particles)

            if debug: t0=time.time()
            velocità=np.array([])
            if n_particles>100_000_000:
                for i in range(int(n_particles/100000000)+int(n_particles%100000000)):
                    if i==0:
                        velocità=norm_Max_Boltz_2D(100000000,temperature,material)[:]
                    else:
                        velocità=np.concatenate((velocità,norm_Max_Boltz_2D(100000000,temperature,material)[:]), axis=0)
                velocità=velocità[0:n_particles]
            else:
                velocità=norm_Max_Boltz_2D(n_particles,temperature,material)
            eff_len=np.real(leng[material]*(1+therm_dil[material]*(temperature-273.15)*0.66))*scaling
            velocità=np.array(velocità)/((eff_len*np.power(elastic_const[material]/(mass[material]*scaling**2),0.5))*scaling**2)

            if debug: 
                print(velocità)
                t1=time.time()
                v_sq=[]
                print("tempo produzione",t1-t0)
                v_sq=np.power(velocità[:, 0]**2+velocità[:, 1]**2,0.5)
                t2=time.time()
                print("tempo totale",t2-t0)
                plt.hist(v_sq, bins=100, density=True, alpha=0.7, label='2D Maxwell')
                plt.xlabel('|v|')
                plt.ylabel('Distribuzione')
                plt.title('Distribuzione Maxwell-Boltzmann 2D')
                plt.legend()
                plt.grid(True)
                plt.show()
            point=0
            for i in range(initial_positions[0], final_positions[0],1):
                vectorsmom=[]
                for j in range(initial_positions[1], final_positions[1],1):
                    if((i==initial_positions[0]) or (j==initial_positions[1])):
                        vectorsmom.append(np.array([i,j,0.,0.],dtype=np.float64))
                    elif((i==final_positions[0]-1) or (j==final_positions[1]-1)):
                        vectorsmom.append(np.array([i,j,0.,0.],dtype=np.float64))
                    else:
                        vectorsmom.append(np.array([i,j,velocità[point][0],velocità[point][1]],dtype=np.float64))
                    point+=1
                vectors.append(np.array(vectorsmom,dtype=np.float64))
                
        else:
            print("Dimension not supported")
    else:
        if (len(initial_positions)==3):
            for i in range(initial_positions[0], final_positions[0],1):
                for j in range(initial_positions[1], final_positions[1],1):
                    for k in range(initial_positions[1], final_positions[1],1):
                        vectors.append([i,j,k,0,0,0])
        elif(len(initial_positions)==2):
            for i in range(initial_positions[0], final_positions[0],1):
                for j in range(initial_positions[1], final_positions[1],1):
                    vectors.append([i,j,0,0])
        else:
            print("Dimension not supported")
    return np.array(vectors)

###################################################    FUNCTION FOR DIRECT LATTICE   ###################################################
""" You can insert 3 different types of objects (at least in the 2D code) which are still classes or subclasses 2 objects that inherit 
    each other's properties. These are used to create groups of particles with desired shapes and positions of selected materials that 
    will then be inserted into the lattice via insert_object_2D.

    They are a rectangle, an ellipse and a random object of elliptical shape (that is: you select an element, and then, based on a random 
    choice, this inserts the intruder data into the dual lattice in place of the material, or leaves things as they were before).
    These, keeping track of the position [i,j] give unscaled mass, k, gamma and the transverse coupling factor (0 or 1)."""
class TwoD_object:
    def __init__(self, position, material, aggregation, scaling, inclusion_ratio, objtemperature, *args):
        self.position = position
        self.material = material
        self.args = args
        self.scaling=scaling
        self.aggregation=aggregation
        self.objtemperature=objtemperature
        self.inclusion_ratio=inclusion_ratio
        try:
            if isinstance(args[0][0],str):
                self.geometry = args[0][0]
                self.geometric_dimensions= args[0][1]
            else:
                self.geometry = "undefined"
                self.geometric_dimensions= args[0][0]
        except:
            pass

    def describe(self):
        return f"Two_object at ({self.position[0]}, {self.position[0]}) made of {self.material} ({self.aggregation}) with geometry {self.geometry}, and geometric dim ({self.geometric_dimensions[0]},{self.geometric_dimensions[1]}). Scaling factor: {self.scaling}"

    @staticmethod
    def create(position, material, aggregation, scaling, inclusion_ratio, objtemperature,*args):
        if isinstance(args[0][0],str):
            geometry = args[0][0]
            geometric_dimensions= args[0][1]
        else:
            geometry = "undefined"
            geometric_dimensions= args[0][0]

        if geometry=="rectangle":
            return rectangle(position, material, aggregation, scaling, inclusion_ratio, objtemperature, args)
        elif geometry=="ellipse":
            return ellipse(position, material, aggregation, scaling, inclusion_ratio, objtemperature, args)
        else:
            return TwoD_object(position, material, aggregation, scaling, inclusion_ratio, objtemperature, args)
    
    def in_position(self,esterni):
        print("Undefined twoD_obj defined... ")
        result=[]
        if isinstance(self.args[0][0],str):
            #geometry = self.args[0][0]
            geometric_dimensions= self.args[0][1]
            try:
                estrema=esterni[1]
                assex=int(self.geometric_dimensions[1][0])
                assey=int(self.geometric_dimensions[1][1])
                try:
                    assez=int(self.geometric_dimensions[1][2])
                except:
                    pass
                geodim=len(self.geometric_dimensions[1])
            except:
                estrema=esterni
                assex=int(self.geometric_dimensions[0])
                assey=int(self.geometric_dimensions[1])
                try:
                    assez=int(self.geometric_dimensions[2])
                except:
                    pass
                geodim=len(self.geometric_dimensions)
        else:
            geometric_dimensions= self.args[0][0]
            try:
                estrema=esterni[1]
                assex=int(self.geometric_dimensions[1][0])
                assey=int(self.geometric_dimensions[1][1])
                try:
                    assez=int(self.geometric_dimensions[1][2])
                except:
                    pass
                geodim=len(self.geometric_dimensions[1])
            except:
                estrema=esterni
                assex=int(self.geometric_dimensions[0])
                assey=int(self.geometric_dimensions[1])
                try:
                    assez=int(self.geometric_dimensions[2])
                except:
                    pass
                geodim=len(self.geometric_dimensions)
                
            if geodim==2:

                for i in range(self.position[0]-int(estrema[0]),self.position[0]+int(estrema[0]),1):
                    for j in range(self.position[1]-int(estrema[1]),int(estrema[1])+self.position[1],1):
                        if (((i-self.position[0])**2)/assex**2 +((j-self.position[1])**2)/assey**2)<=1:
                                if self.aggregation=="solid":
                                    km=elastic_const[self.material]# *self.scaling The function scales linearly in 3D, approximately uniformly in 2D.
                                    gamm=gamma0[self.material]*self.scaling**(3/2)# TO MODIFY SCALING FUNCTION FOR GAMMA HERE AND IN insert_object_2D!!!
                                    for i in range(self.position[0]-int(estrema[0]),self.position[0]+int(estrema[0]),1):
                                        for j in range(self.position[1]-int(estrema[1]),int(estrema[1])+self.position[1],1):
                                            if (((i-self.position[0])**2)/assex**2 +((j-self.position[1])**2)/assey**2)<=1:
                                                if(np.random.randint(0,2)==1): 
                                                    result.append([[i,j],[True,mass[self.material]*self.scaling**(2),km,gamm,1]])
                                elif self.aggregation=="fluid":
                                    km=elastic_const[self.material]# *self.scaling The function scales linearly in 3D, approximately uniformly in 2D.
                                    gamm=gamma0[self.material]*self.scaling**(3/2)# TO MODIFY SCALING FUNCTION FOR GAMMA HERE AND IN insert_object_2D!!!
                                    for i in range(self.position[0]-int(estrema[0]),self.position[0]+int(estrema[0]),1):
                                        for j in range(self.position[1]-int(estrema[1]),int(estrema[1])+self.position[1],1):
                                            if (((i-self.position[0])**2)/assex**2 +((j-self.position[1])**2)/assey**2)<=1:
                                                if(np.random.randint(0,2)==1): 
                                                    result.append([[i,j],[True,mass[self.material]*self.scaling**(2),km,gamm,0]])
                                else:
                                    print("WARNING: NOT AN AGGREGATION STATE!")
                                    result.append([[i,j],[True,1,0,0,0]])

            if geodim==3:
                print("Dimension not supported in TwoD_object in_position function")
            return result


class rectangle(TwoD_object):
    def __init__(self, position, material, aggregation, scaling, inclusion_ratio, objtemperature, *args):
        super().__init__(position, material, aggregation, scaling, inclusion_ratio, objtemperature)
        self.geometric_dimensions = args[0][0]
        
    def in_position(self,esterni):
        try:
            esterni[1][0]
            estrema=esterni[1]
        except:
            estrema=esterni
        print("Rectangle defined ...")
        if len(self.geometric_dimensions[1])==2:
            result=[]
            km=elastic_const[self.material]#   no scal of k in 2d
            gamm=gamma0[self.material]*(self.scaling**(3/2))# TO MODIFY SCALING FUNCTION HERE AND IN insert_object_2D!!!
            if self.aggregation=="solid":
                print("Kind Solid")
                for i in range(self.position[0]-int(estrema[0]),self.position[0]+int(estrema[0]),1): # I put the for inside to save on the if statements
                    for j in range(self.position[1]-int(estrema[1]),int(estrema[1])+self.position[1],1):
                        result.append([[i,j],[True,mass[self.material]*(self.scaling**2),km,gamm,1]])
            elif self.aggregation=="fluid":
                print("Kind fluid")
                for i in range(self.position[0]-int(estrema[0]),self.position[0]+int(estrema[0]),1): # I put the for inside to save on the if statements
                    for j in range(self.position[1]-int(estrema[1]),int(estrema[1])+self.position[1],1):
                        result.append([[i,j],[True,mass[self.material]*(self.scaling**2),km,gamm,0]])
            else:
                        print("WARNING: NOT AN AGGREGATION STATE!")
                        result.append([[0,0],[True,1,0,0,0]])
                    
            return result
        if len(self.geometric_dimensions[1])==3:
            print("Dimension not supported in rectangle in_position function")
            return False


    def describe(self):
        return (f"Rectangle at ({self.position[0]}, {self.position[0]}) made of {self.material} ({self.aggregation}), "
                f"with sides a = {self.geometric_dimensions[0]} and b = {self.geometric_dimensions[1]}. Scaling factor: {self.scaling}")


class ellipse(TwoD_object):
    def __init__(self, position, material, aggregation, scaling, inclusion_ratio, objtemperature, *args):
        super().__init__(position, material, aggregation, scaling, inclusion_ratio, objtemperature)
        self.geometric_dimensions = args[0][0]

    def describe(self):
        return (f"Ellipse at ({self.position[0]}, {self.position[0]}) made of {self.material} ({self.aggregation}), "
                f"with sides a = {self.geometric_dimensions[0]} and b = {self.geometric_dimensions[1]}. Scaling factor: {self.scaling}")
    
    def in_position(self,esterni):
        print("Ellipse defined ...")
        try:
            estrema=esterni[1]
            assex=estrema[0]
            assey=estrema[1]
            try:
                assez=int(self.geometric_dimensions[1][2])
            except:
                pass
            geodim=len(self.geometric_dimensions[1])
        except:
            estrema=esterni
            assex=estrema[0]
            assey=estrema[1]
            assez=estrema[2]
            geodim=len(self.geometric_dimensions)

        if geodim==2:
            
            result=[]
            km=elastic_const[self.material]#   no scale of k in 2d
            gamm=gamma0[self.material]*(self.scaling**(3/2))# TO MODIFY SCALING FUNCTION HERE AND IN insert_object_2D!!!
            if self.aggregation=="solid":
                print("Kind solid")
                for i in range(self.position[0]-int(estrema[0]),self.position[0]+int(estrema[0]),1):
                    for j in range(self.position[1]-int(estrema[1]),int(estrema[1])+self.position[1],1):
                        if (((i-self.position[0])**2)/assex**2 +((j-self.position[1])**2)/assey**2)<=1:
                            result.append([[i,j],[True,mass[self.material]*(self.scaling**2),km,gamm,1]]) 
            elif self.aggregation=="fluid":
                print("Kind fluid")
                for i in range(self.position[0]-int(estrema[0]),self.position[0]+int(estrema[0]),1):
                    for j in range(self.position[1]-int(estrema[1]),int(estrema[1])+self.position[1],1):
                        if (((i-self.position[0])**2)/assex**2 +((j-self.position[1])**2)/assey**2)<=1:
                            result.append([[i,j],[True,mass[self.material]*(self.scaling**2),km,gamm,0]]) 
            else:
                print("WARNING: NOT AN AGGREGATION STATE!")
                result.append([[0,0],[True,1,0,0,0]])
            return result
        if geodim==3:
            print("Dimension not supported in ellipse in_position function")
            return 0

###################################################    BUILD LATTICE FUNCTION   ###################################################
"""This puts the elements of the direct lattice in place (totally) and then builds the entire dual lattice of the basic substance but 
   does not insert the intruders because those have a scaling that has yet to be evaluated."""

def build_lattice(dims,material,state,temperature, scaling,obj_list=[]):
    if len(dims)==3:
        lattice_positions=propose_in_position([0,0,0],dims,material,state,temperature, scaling)
    if len(dims)==2:
        lattice_positions=propose_in_position([0,0],dims,material,state,temperature, scaling)
    lattice_info=[]
    km=1
    mu=1
    if len(dims)==2:
        gamm=(gamma0[material]*scaling**(3/2))/((mass[material]*(scaling**2))*np.power(elastic_const[material]/mass[material]*(scaling**2),0.5)*scaling)
        if obj_list==[]:
            for i in range(dims[0]):
                lattice_info_inf=[]
                for j in range(dims[1]):
                    if state=="solid":
                                    lattice_info_inf.append([mu,km,gamm,1])
                    if state=="liquid":
                                    lattice_info_inf.append([mu,km,gamm,0])
                lattice_info.append(lattice_info_inf)
        else:
            for i in range(dims[0]):
                lattice_info_inf=[]
                for j in range(dims[1]):
                    for geo in obj_list:
                        pos_geo=geo.in_position([i,j])
                        if pos_geo[0]:
                            lattice_info_inf.append(pos_geo[1:4])
                        else:
                            if state=="solid":
                                    lattice_info_inf.append([mu,km,gamm,1])
                            if state=="liquid":
                                    lattice_info_inf.append([mu,km,gamm,0])
                lattice_info.append(lattice_info_inf)
    if len(dims)==3:
        gamm=(gamma0[material]*scaling**(5/2))/((mass[material]*(scaling**2))*np.power(elastic_const[material]/mass[material]*(scaling**2),0.5)*scaling)
        if obj_list==[]:
            for i in range(dims[0]):
                lattice_info_inf1=[]
                for j in range(dims[1]):
                    lattice_info_inf2=[]
                    for k in range(dims[2]):
                        if state=="solid":
                                    lattice_info_inf2.append([mu,km,gamm,1])
                        if state=="liquid":
                                    lattice_info_inf2.append([mu,km,gamm,0])
                        lattice_info_inf1.append(lattice_info_inf2)
                lattice_info.append(lattice_info_inf1)
                    
        else:
            for i in range(dims[0]):
                lattice_info_inf1=[]
                for j in range(dims[1]):
                    lattice_info_inf2=[]
                    for k in range(dims[2]):
                        for geo in obj_list:
                            pos_geo=geo.in_position([i,j,k])
                            if pos_geo[0]:
                                lattice_info_inf2.append(pos_geo[1:4])
                            else:
                                if state=="solid":
                                    lattice_info_inf2.append([mu,km,gamm,1])
                                if state=="liquid":
                                    lattice_info_inf2.append([mu,km,gamm,0])
                    lattice_info_inf1.append(lattice_info_inf2)
                lattice_info.append(lattice_info_inf1)
    return lattice_positions, lattice_info

###############################################   FUNCTION TO INSERT THE EXTERNAL OBJs  ################################################
"""The external objects that are communicated by the LSM need a dedicated scaling also based on how many particles must be inserted in how 
   many places. Furthermore, these can have different temperatures and certainly different effective temperatures (and the scaling is not 
   necessarily invertible so you cannot put an effective temperature in a non-rigorous way).
   So it must act first on the direct lattice to normalize the equations and then on the dual to send the info to C_Lattice."""

def insert_object_2D(lattice,objec,basetemperature,basematerial,scaling):

    eff_len_base=np.real(leng[basematerial]*(1+(therm_dil[basematerial]*(basetemperature-273.15))**(2/3)))
    eff_len_obj=np.real(leng[objec.material]*(1+(therm_dil[objec.material]*(objec.objtemperature-273.15))**(2/3)))
    substitutional_scaling=(1/objec.inclusion_ratio)*(eff_len_base/eff_len_obj)**2
    modifiche=objec.in_position(objec.geometric_dimensions)
    massa1=mass[basematerial]*(scaling**2)
    k1=elastic_const[basematerial]
    omeg1=np.power(elastic_const[basematerial]/(mass[basematerial]*(scaling**2)),0.5)
    for m in modifiche:
        lattice[1][m[0][0]][m[0][1]][0]=m[1][1]/massa1*substitutional_scaling
        lattice[1][m[0][0]][m[0][1]][1]=m[1][2]/k1 #2d case: no scalig at all
        lattice[1][m[0][0]][m[0][1]][2]=m[1][3]/(massa1*omeg1)*(substitutional_scaling)**(3/2)
        lattice[1][m[0][0]][m[0][1]][3]=m[1][4]
    #The speed also needs to be adjusted because it has a different normalization
    total_scaling=scaling*substitutional_scaling
    
    temperature=objec.objtemperature

    velocità=np.array([])
    velocità=norm_Max_Boltz_2D(len(modifiche),temperature,objec.material)
    velocità=np.array(velocità)/((eff_len_obj*total_scaling*np.power(elastic_const[objec.material]/(mass[objec.material]*total_scaling**2),0.5))*scaling**2)

    for i,m in enumerate(modifiche):
        lattice[0][m[0][0]][m[0][1]][2]=velocità[i][0]
        lattice[0][m[0][0]][m[0][1]][3]=velocità[i][1]
    return lattice

def scale_back2D(substance,state,temperature,scaling):
    scaled_len=(scaling**2)*np.real(leng[substance]*(1+(therm_dil[substance]*(temperature-273.15))**(2/3)))
    scaled_mass=mass[substance]*(scaling**2)
    scaled_k=elastic_const[substance]
    scaled_gamm=gamma0[substance]*(scaling**0.5)
    scaled_omega=np.power(elastic_const[substance]/(mass[substance]*(scaling**2)),0.5)
    return scaled_len, scaled_mass, scaled_k, scaled_gamm, scaled_omega
