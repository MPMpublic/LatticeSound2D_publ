*This code is property of Giorgio Lo Presti (see MPMpublic) and distributed under GPL license as 
 it is intended to be open-source with the condition to use it for open-source 
 and research purposes*

Description:
LatticeSound 2D is a high-performances two dimensional multi-material meso-scalable 
multithread code for mixed oscillation modes in Lattice Thermal Quantum Field Theory 
of phonons.
It can simulate the sound, vibration or noise propagation inside a 2 dimensional 
lattice initialized with different materials of different shapes.

Theory ensures that scaling is performed in such a way as to maintain the results at
the correct resonance frequencies, so that the behavior of the aggregates reflects a
consistent stress distribution at the mesoscopic scale. It is possible to extend the
 size of the body by increasing the number of oscillation points per axis, or by 
increasing the scale factor (decreasing the sensitivity), or by increasing the 
sensitivity by performing the process in reverse.

As the scaling or number of centers varies, the system changes its effective mass 
and coupling to the electromagnetic field due to the running coupling.

The system propagates the vibration to each oscillation center by evaluating the 
position and velocity of its nearest neighbors at the previous instant. No approximation 
is imposed in the discretization of the wave equation. Resolution of the dynamics of
 each oscillator is guaranteed by tested libraries (GSL) written in C, which handle 
a limited portion of the lattice based on the number of threads selected. 
Transmission, calculations, and storage of steps by all cores are synchronized 
through barriers to prevent data corruption.

The code works at finite temperature (to evaluate the dynamics at zero thermal it is 
necessary to impose positive temperatures close to zero) and at ambient pressure, and
 each oscillation point can be customized both in terms of material and initial 
condition, even by imposing overpressures through the proportion factor of the 
intrusive impurities.

By selecting a forcing on the lower boundary, the dynamics of the forced system and 
the effective resonances and masses can be calculated. The forcing can be modulated 
in amplitude and frequency in a completely customizable way.
The lateral ends are fixed like the upper bound.


The purchase of thermal energy from the external environment can be selected using 
the thermal_coupling factor.

Warning: The number of oscillators must be greater than the number of threads 
selected. Do not go below a 3x3 system. The plot window must be smaller than the 
number of oscillators. Do not run plots with systems larger than 60x60 due to the 
computation time required for the images (results will still be provided).

The code automatically performs an initial data analysis, reporting the output data 
at a user-defined time interval. For each print step, for each thread, information 
is provided on the selected oscillators and the average displacement, velocity, 
force, and related offsets for all the oscillation centers assigned to that thread 
(even if manually excluded from the plot window, but in that case, the information 
on the individual oscillators will not be visible).

IMPORTANT:
This code is free-to-use and no support or responsibility is borne by the author.
This code is distributed WITHOUT ANY WARRANTY; without even the implied warranty 
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. The results of the code 
are not intended to be exhaustive and should be reviewed by experienced personnel 
in the field. The author declines all responsibility for the consequences resulting 
from improper use of the code. Thus:
The results produced by this software should not be considered authoritative. 
Decisions, analyses, or actions based solely on the output of this program are taken
entirely at the user's own risk.

You can ask for support at giorgio.lopresti@dfa.unict.it or 
giorgiolopresti@lns.infn.it

Installation:
Please, make sure you have installed the following packages before you continue
with install procedure:
• python3,
• GCC compiler,
• GSL-GNU,
• C Python,
• sysv_ipc,
• common python libraries as: subprocess, multiprocessing, matplotlib, os, sys, time, glob, numpy.

At this point, copy (or pull by git) the repository on the desired folder (possibly
not a system folder).
Advise: create a results folder. Remove previous file of results from results folder 
or create another one before relaunchthe code.

Launch:
-Organize the starter file (please: write the absolute path for results folder in the 
 folder path line in starter file)
-go to bash/prompt, go to the installation folder and launch as: 
                            python3 LSM.py

Use:
see manual

Licenses:
This project is licensed under the GNU General Public License v3.0.  
See the [LICENSE](./LICENSE) file for details.

See third party licenses and watch websites for more information
GCC compiler: https://gcc.gnu.org/git.html
GSL-GNU:https://www.gnu.org/software/gsl/#documentation
C Python: https://github.com/martinohanlon/c_python_ipc/blob/master/LICENSE
sysv_ipc: https://pypi.org/project/sysv-ipc/