# ---------------------------------------------------------------------#
# TUTORIAL to simulate a TriLMP membrane patch                         #
# Author: Maitane Muñoz-Basagoiti (maitane.munoz-basagoiti@ista.ac.at) #
#                                                                      #
# This code launches TriLMP for flat patch with fixed boundary         #
# conditions. The patch is initialized as a triangulated lattice       #
# which has been previously generated as is stored in                  #
# - mesh_coordinates.dat (constains positions of vertices)             #
# - mesh_faces.dat       (contains the oriented faces of the network)  #
# The network has 11500 vertices.                                      #
# ---------------------------------------------------------------------#

# .....................................................................#
#       PRELIMINARY STEPS NEEDED TO RUN THE TRILMP SIMULATION
# .....................................................................#

# Start by importing the necessary modules
import trimesh
import numpy as np
import pandas as pd
from pathlib import Path
from trimem.mc.trilmp import TriLmp

# Generate the 'checkpoints' directory 
# In this directory, we will save the (pickle) checkpoints
# to use in the second part of the simulation, i.e.,
# configurations from which we will let the system
# equilibrate
Path("checkpoints").mkdir(exist_ok=True)

# Initialize the network
# Remember that the coordinates of the network vertices are in mesh_coordinates.dat
# while mesh_faces contains the oriented faces (that is, faces with all normal vectors
# pointing in the same sense). 
mesh_coordinates = pd.read_csv('mesh_coordinates.dat', header = None, index_col = False, sep = ' ')
mesh_coordinates_array = mesh_coordinates[[0, 1, 2]].to_numpy()
mesh_faces = pd.read_csv('mesh_faces.dat', header = None, index_col = False, sep = ' ')
mesh_faces_array = mesh_faces[[0, 1, 2]].to_numpy()
mesh = trimesh.Trimesh(vertices=mesh_coordinates_array, faces = mesh_faces_array) 

# the number of vertices in the network
N = len(mesh.vertices)

# Rescale mesh distances
# As the vertices at the edge of the membrane patch are fixed
# in space, establishing the average edge distance in the system
# sets the projected area in the simulation. Note that the
# desired_average_distance has to be in [1.0, sqrt(3)],
# which correspond to the minimum and maximum values the edges
# of the network can adopt before incurring in energy penalties
old_average_distance = np.mean(mesh.edges_unique_length)
desired_average_distance = 1.05
scaling = desired_average_distance/old_average_distance
mesh.vertices *=scaling
current_average_distance = np.mean(mesh.edges_unique_length)
print('Old average edge length: ', current_average_distance)
print('Current average edge length: ', current_average_distance)

# You can print out the properties of the mesh
# you will use to initialize your simulation
print("Vertex properties: ")
print(f"MESH VERTICES : ", len(mesh.vertices))
print(f"MESH FACES    : ", len(mesh.faces))
print(f"MESH EDGES    : ", len(mesh.edges))

print("")
print("** STARTING TRILMP SIMULATION **")
print("")

# .....................................................................#
#                        TRILMP SIMULATION
# .....................................................................#
#
#
# Running a TriLMP simulation involves 3 stages:
#
# 1. First, one creates the TriLmp object. 
# This step initializes the simulation set-up. Here you define...
#       - The number of particles in the system, their positions, types and properties
#               + By default, the particles in fluid membrane are called 'vertices'
#               + The bead tethered to the membrane, that will pull on it, is called 'bead' here
#       - The dimensions of the simulation box and its boundary conditions (periodic or not)
#       - The mechanical properties of the membrane (i.e., the parameters in the TriMEM hamiltonian)
#       - The properties of the Hybrid Monte Carlo integration scheme (MD stages and MC stages)
#               + TriLMP will alternate randomly between MD (vertex motion) and MC (connectivity update) stages
#               + You can decide how many steps to integrate during the MD stages
#                 as well as what fraction of network edges to attempt flipping in the MC stages
#                 (see Siggel et al., 2022, for a discussion on the limitations of this parameter)
#       - The length of the equilibration stage
#       - The frequency at which performance, system info is outputted and checkpoints are created by TriLMP
#
# 2. Second, one adds LAMMPS commands to the simulation
# Here you can choose to pass as many LAMMPS commands as desired to your simulation.
# Note that there are some compulsory commands (see pair_styles or time-integration sections below)
# which must always be passed to ensure that TriLMP runs correctly.
# The LAMMPS commands can be passed in two different ways:
# - Pre-equilibration commands
#       + Pass them directly to the TriLmp object before the simulation has started running
#         (this is useful to equilibrate your system beforehand)
# - Post-equilibration commands
#       + Pass them by appending them to the postequilibration_commands list below
#          (they will be run once the pre-equilibration stage has ended, in order)
#
# 3. Third, run the simulation by calling the TriLmp.run(args) method with its corresponding arguments
#

# 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
# 1. CREATE THE TriLmp OBJECT

trilmp=TriLmp(initialize=True,                            # Use mesh above to initialize reference state of the system
            debug_mode=False,                             # For DEBUGGING purposes. Set it to True if you want the complete LAMMPS output
            periodic=True,                                # Boundary conditions of the system. Set to true if you want periodic boundaries.
            num_particle_types=2,                         # Total number of particle species in system (here two types: membrane type and bead type)
            mass_particle_type=[1.0, 1.0],                # Define the mass of the particles species in system
            group_particle_type=['vertices', 'bead'],     # Assign the different particle species a group, to be used below with LAMMPS commands

            mesh_points=mesh.vertices,                    # Introduce the mesh vertices
            mesh_faces=mesh.faces,                        # Introduce the mesh faces
            
            # Introduce the IDs of the vertices in the network at the edge of the patch
            # These vertices will remain fixed in space throughout the simulation
            vertices_at_edge=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 199, 200, 299, 300, 399, 400, 499, 500, 599, 600, 699, 700, 799, 800, 899, 900, 999, 1000, 1099, 1100, 1199, 1200, 1299, 1300, 1399, 1400, 1499, 1500, 1599, 1600, 1699, 1700, 1799, 1800, 1899, 1900, 1999, 2000, 2099, 2100, 2199, 2200, 2299, 2300, 2399, 2400, 2499, 2500, 2599, 2600, 2699, 2700, 2799, 2800, 2899, 2900, 2999, 3000, 3099, 3100, 3199, 3200, 3299, 3300, 3399, 3400, 3499, 3500, 3599, 3600, 3699, 3700, 3799, 3800, 3899, 3900, 3999, 4000, 4099, 4100, 4199, 4200, 4299, 4300, 4399, 4400, 4499, 4500, 4599, 4600, 4699, 4700, 4799, 4800, 4899, 4900, 4999, 5000, 5099, 5100, 5199, 5200, 5299, 5300, 5399, 5400, 5499, 5500, 5599, 5600, 5699, 5700, 5799, 5800, 5899, 5900, 5999, 6000, 6099, 6100, 6199, 6200, 6299, 6300, 6399, 6400, 6499, 6500, 6599, 6600, 6699, 6700, 6799, 6800, 6899, 6900, 6999, 7000, 7099, 7100, 7199, 7200, 7299, 7300, 7399, 7400, 7499, 7500, 7599, 7600, 7699, 7700, 7799, 7800, 7899, 7900, 7999, 8000, 8099, 8100, 8199, 8200, 8299, 8300, 8399, 8400, 8499, 8500, 8599, 8600, 8699, 8700, 8799, 8800, 8899, 8900, 8999, 9000, 9099, 9100, 9199, 9200, 9299, 9300, 9399, 9400, 9499, 9500, 9599, 9600, 9699, 9700, 9799, 9800, 9899, 9900, 9999, 10000, 10099, 10100, 10199, 10200, 10299, 10300, 10399, 10400, 10499, 10500, 10599, 10600, 10699, 10700, 10799, 10800, 10899, 10900, 10999, 11000, 11099, 11100, 11199, 11200, 11299, 11300, 11399, 11400, 11401, 11402, 11403, 11404, 11405, 11406, 11407, 11408, 11409, 11410, 11411, 11412, 11413, 11414, 11415, 11416, 11417, 11418, 11419, 11420, 11421, 11422, 11423, 11424, 11425, 11426, 11427, 11428, 11429, 11430, 11431, 11432, 11433, 11434, 11435, 11436, 11437, 11438, 11439, 11440, 11441, 11442, 11443, 11444, 11445, 11446, 11447, 11448, 11449, 11450, 11451, 11452, 11453, 11454, 11455, 11456, 11457, 11458, 11459, 11460, 11461, 11462, 11463, 11464, 11465, 11466, 11467, 11468, 11469, 11470, 11471, 11472, 11473, 11474, 11475, 11476, 11477, 11478, 11479, 11480, 11481, 11482, 11483, 11484, 11485, 11486, 11487, 11488, 11489, 11490, 11491, 11492, 11493, 11494, 11495, 11496, 11497, 11498, 11499],
            
            # Choose area constraint version
            # If area_frac == 1,  TriLMP uses a quadratic constraint on the area: kappa_a (1-A/A_0)^2
            # If area_frac == -1, TriLMP uses a linear constraint on the area   : kappa_a A
            area_frac=-1.0,                       

            # To get further details on the mechanical parameters below, please see Siggel et al., 2022
            kappa_b=20,                            # MEMBRANE MECHANICS: bending modulus (kB T)
            kappa_a=7.0,                           # MEMBRANE MECHANICS: surface tension (kB T/sigma^2)
            kappa_v=0,                             # MEMBRANE MECHANICS: constraint on volume change from target value (kB T)
            kappa_c=0.0,                           # MEMBRANE MECHANICS: constraint on area difference change (kB T)
            kappa_t=10000.0,                       # MEMBRANE MECHANICS: tethering potential to constrain edge length (kB T)
            kappa_r=1000.0,                        # MEMBRANE MECHANICS: repulsive potential to prevent surface intersection (kB T)

            # Define the number of bonds types in the system - Here 2 because:
            # - Type 1: Membrane edges
            # - Type 2: Membrane-bead tethering
            n_bond_types=2,  

            # Information regarding 'non-membrane' particles
            n_beads=1,                                          # Number of non-membrane beads
            n_bead_types=1,                                     # Number of non-membrane bead-types
            bead_pos=np.array([[0, 0, 6.050000000000001]]),     # Position of each of the non-membrane bead-types
            bead_types=np.array([2]),                           # Specify the type of each of the non-membrane beads (By default, start by assigning type '2'; type '1' is for the membrane)
            bead_sizes=np.array([10]),                          # Specify the diameter of each of the non-membrane beads                   

            # Simulation parameters
            step_size=0.001,                       # MD PART: Timestep for the time integration in the MD stages of the simulation
            traj_steps=100,                        # MD PART: Length of the MD stages in simulation steps
            flip_ratio=0.2,                        # MC PART: MC PART SIMULATION: fraction of edges to flip
            initial_temperature=1.0,               # MD/MC PART: (Initial) temperature of the system
            switch_mode="random",                  # MD/MC PART: Protocol to alternate sequence of MD and MC stages in the simulation. Options: 'random' or 'alternating' flip-or-move

            # Size of the simulation box
            box=(-59.74782608695652, 59.75217391304348, -59.36344801571301, 59.36344801571298, -20, 220), 

            # Number of equilibration rounds. Must be a multiple of traj_steps
            equilibration_rounds=100, 
        
            info=5000,                              # OUTPUT: frequency at which TriMEM information is printed in shell/terminal
            performance_increment=1000,             # OUTPUT: frequency at which performace stats in are printed in inp_performance.dat (MD+MC stage FREQUENCY)
            energy_increment=1000,                  # OUTPUT: frequency at which mesh properties in are printed in inp_system.dat (MD+MC stage FREQUENCY)
            checkpoint_every=100*5000,              # OUTPUT: frequency at which (pickle) checkpoints of the system are generated (MD+MC stage FREQUENCY)
   
            )

# 222222222222222222222222222222222222222222222222222222222222
# 2. INTRODUCE LAMMPS COMMANDS

# PRE-PRE-PRE-PRE-PRE-PRE-PRE-PRE-PRE-PRE-PRE-PRE-PRE
# Pre-equilibration stage: LAMMPS commands to pass
# before the simulation has stated
# PRE-PRE-PRE-PRE-PRE-PRE-PRE-PRE-PRE-PRE-PRE-PRE-PRE

# ..............................................................
# GROUPS
# [NECESSARY] 
# ..............................................................
# Aside from the 'vertices' and 'bead' groups defined in 
# the TriLmp object above, we will define two more groups:
# - 'vertex_edge' group : for vertices at the edge of the patch
# - BULK group : for vertices NOT at the edge of the patch
# ..............................................................

#                        --- START ---

# Vertices at the edge
trilmp.lmp.command('group vertex_edge id 1')
trilmp.lmp.command('group vertex_edge id 2')
trilmp.lmp.command('group vertex_edge id 3')
trilmp.lmp.command('group vertex_edge id 4')
trilmp.lmp.command('group vertex_edge id 5')
trilmp.lmp.command('group vertex_edge id 6')
trilmp.lmp.command('group vertex_edge id 7')
trilmp.lmp.command('group vertex_edge id 8')
trilmp.lmp.command('group vertex_edge id 9')
trilmp.lmp.command('group vertex_edge id 10')
trilmp.lmp.command('group vertex_edge id 11')
trilmp.lmp.command('group vertex_edge id 12')
trilmp.lmp.command('group vertex_edge id 13')
trilmp.lmp.command('group vertex_edge id 14')
trilmp.lmp.command('group vertex_edge id 15')
trilmp.lmp.command('group vertex_edge id 16')
trilmp.lmp.command('group vertex_edge id 17')
trilmp.lmp.command('group vertex_edge id 18')
trilmp.lmp.command('group vertex_edge id 19')
trilmp.lmp.command('group vertex_edge id 20')
trilmp.lmp.command('group vertex_edge id 21')
trilmp.lmp.command('group vertex_edge id 22')
trilmp.lmp.command('group vertex_edge id 23')
trilmp.lmp.command('group vertex_edge id 24')
trilmp.lmp.command('group vertex_edge id 25')
trilmp.lmp.command('group vertex_edge id 26')
trilmp.lmp.command('group vertex_edge id 27')
trilmp.lmp.command('group vertex_edge id 28')
trilmp.lmp.command('group vertex_edge id 29')
trilmp.lmp.command('group vertex_edge id 30')
trilmp.lmp.command('group vertex_edge id 31')
trilmp.lmp.command('group vertex_edge id 32')
trilmp.lmp.command('group vertex_edge id 33')
trilmp.lmp.command('group vertex_edge id 34')
trilmp.lmp.command('group vertex_edge id 35')
trilmp.lmp.command('group vertex_edge id 36')
trilmp.lmp.command('group vertex_edge id 37')
trilmp.lmp.command('group vertex_edge id 38')
trilmp.lmp.command('group vertex_edge id 39')
trilmp.lmp.command('group vertex_edge id 40')
trilmp.lmp.command('group vertex_edge id 41')
trilmp.lmp.command('group vertex_edge id 42')
trilmp.lmp.command('group vertex_edge id 43')
trilmp.lmp.command('group vertex_edge id 44')
trilmp.lmp.command('group vertex_edge id 45')
trilmp.lmp.command('group vertex_edge id 46')
trilmp.lmp.command('group vertex_edge id 47')
trilmp.lmp.command('group vertex_edge id 48')
trilmp.lmp.command('group vertex_edge id 49')
trilmp.lmp.command('group vertex_edge id 50')
trilmp.lmp.command('group vertex_edge id 51')
trilmp.lmp.command('group vertex_edge id 52')
trilmp.lmp.command('group vertex_edge id 53')
trilmp.lmp.command('group vertex_edge id 54')
trilmp.lmp.command('group vertex_edge id 55')
trilmp.lmp.command('group vertex_edge id 56')
trilmp.lmp.command('group vertex_edge id 57')
trilmp.lmp.command('group vertex_edge id 58')
trilmp.lmp.command('group vertex_edge id 59')
trilmp.lmp.command('group vertex_edge id 60')
trilmp.lmp.command('group vertex_edge id 61')
trilmp.lmp.command('group vertex_edge id 62')
trilmp.lmp.command('group vertex_edge id 63')
trilmp.lmp.command('group vertex_edge id 64')
trilmp.lmp.command('group vertex_edge id 65')
trilmp.lmp.command('group vertex_edge id 66')
trilmp.lmp.command('group vertex_edge id 67')
trilmp.lmp.command('group vertex_edge id 68')
trilmp.lmp.command('group vertex_edge id 69')
trilmp.lmp.command('group vertex_edge id 70')
trilmp.lmp.command('group vertex_edge id 71')
trilmp.lmp.command('group vertex_edge id 72')
trilmp.lmp.command('group vertex_edge id 73')
trilmp.lmp.command('group vertex_edge id 74')
trilmp.lmp.command('group vertex_edge id 75')
trilmp.lmp.command('group vertex_edge id 76')
trilmp.lmp.command('group vertex_edge id 77')
trilmp.lmp.command('group vertex_edge id 78')
trilmp.lmp.command('group vertex_edge id 79')
trilmp.lmp.command('group vertex_edge id 80')
trilmp.lmp.command('group vertex_edge id 81')
trilmp.lmp.command('group vertex_edge id 82')
trilmp.lmp.command('group vertex_edge id 83')
trilmp.lmp.command('group vertex_edge id 84')
trilmp.lmp.command('group vertex_edge id 85')
trilmp.lmp.command('group vertex_edge id 86')
trilmp.lmp.command('group vertex_edge id 87')
trilmp.lmp.command('group vertex_edge id 88')
trilmp.lmp.command('group vertex_edge id 89')
trilmp.lmp.command('group vertex_edge id 90')
trilmp.lmp.command('group vertex_edge id 91')
trilmp.lmp.command('group vertex_edge id 92')
trilmp.lmp.command('group vertex_edge id 93')
trilmp.lmp.command('group vertex_edge id 94')
trilmp.lmp.command('group vertex_edge id 95')
trilmp.lmp.command('group vertex_edge id 96')
trilmp.lmp.command('group vertex_edge id 97')
trilmp.lmp.command('group vertex_edge id 98')
trilmp.lmp.command('group vertex_edge id 99')
trilmp.lmp.command('group vertex_edge id 100')
trilmp.lmp.command('group vertex_edge id 101')
trilmp.lmp.command('group vertex_edge id 200')
trilmp.lmp.command('group vertex_edge id 201')
trilmp.lmp.command('group vertex_edge id 300')
trilmp.lmp.command('group vertex_edge id 301')
trilmp.lmp.command('group vertex_edge id 400')
trilmp.lmp.command('group vertex_edge id 401')
trilmp.lmp.command('group vertex_edge id 500')
trilmp.lmp.command('group vertex_edge id 501')
trilmp.lmp.command('group vertex_edge id 600')
trilmp.lmp.command('group vertex_edge id 601')
trilmp.lmp.command('group vertex_edge id 700')
trilmp.lmp.command('group vertex_edge id 701')
trilmp.lmp.command('group vertex_edge id 800')
trilmp.lmp.command('group vertex_edge id 801')
trilmp.lmp.command('group vertex_edge id 900')
trilmp.lmp.command('group vertex_edge id 901')
trilmp.lmp.command('group vertex_edge id 1000')
trilmp.lmp.command('group vertex_edge id 1001')
trilmp.lmp.command('group vertex_edge id 1100')
trilmp.lmp.command('group vertex_edge id 1101')
trilmp.lmp.command('group vertex_edge id 1200')
trilmp.lmp.command('group vertex_edge id 1201')
trilmp.lmp.command('group vertex_edge id 1300')
trilmp.lmp.command('group vertex_edge id 1301')
trilmp.lmp.command('group vertex_edge id 1400')
trilmp.lmp.command('group vertex_edge id 1401')
trilmp.lmp.command('group vertex_edge id 1500')
trilmp.lmp.command('group vertex_edge id 1501')
trilmp.lmp.command('group vertex_edge id 1600')
trilmp.lmp.command('group vertex_edge id 1601')
trilmp.lmp.command('group vertex_edge id 1700')
trilmp.lmp.command('group vertex_edge id 1701')
trilmp.lmp.command('group vertex_edge id 1800')
trilmp.lmp.command('group vertex_edge id 1801')
trilmp.lmp.command('group vertex_edge id 1900')
trilmp.lmp.command('group vertex_edge id 1901')
trilmp.lmp.command('group vertex_edge id 2000')
trilmp.lmp.command('group vertex_edge id 2001')
trilmp.lmp.command('group vertex_edge id 2100')
trilmp.lmp.command('group vertex_edge id 2101')
trilmp.lmp.command('group vertex_edge id 2200')
trilmp.lmp.command('group vertex_edge id 2201')
trilmp.lmp.command('group vertex_edge id 2300')
trilmp.lmp.command('group vertex_edge id 2301')
trilmp.lmp.command('group vertex_edge id 2400')
trilmp.lmp.command('group vertex_edge id 2401')
trilmp.lmp.command('group vertex_edge id 2500')
trilmp.lmp.command('group vertex_edge id 2501')
trilmp.lmp.command('group vertex_edge id 2600')
trilmp.lmp.command('group vertex_edge id 2601')
trilmp.lmp.command('group vertex_edge id 2700')
trilmp.lmp.command('group vertex_edge id 2701')
trilmp.lmp.command('group vertex_edge id 2800')
trilmp.lmp.command('group vertex_edge id 2801')
trilmp.lmp.command('group vertex_edge id 2900')
trilmp.lmp.command('group vertex_edge id 2901')
trilmp.lmp.command('group vertex_edge id 3000')
trilmp.lmp.command('group vertex_edge id 3001')
trilmp.lmp.command('group vertex_edge id 3100')
trilmp.lmp.command('group vertex_edge id 3101')
trilmp.lmp.command('group vertex_edge id 3200')
trilmp.lmp.command('group vertex_edge id 3201')
trilmp.lmp.command('group vertex_edge id 3300')
trilmp.lmp.command('group vertex_edge id 3301')
trilmp.lmp.command('group vertex_edge id 3400')
trilmp.lmp.command('group vertex_edge id 3401')
trilmp.lmp.command('group vertex_edge id 3500')
trilmp.lmp.command('group vertex_edge id 3501')
trilmp.lmp.command('group vertex_edge id 3600')
trilmp.lmp.command('group vertex_edge id 3601')
trilmp.lmp.command('group vertex_edge id 3700')
trilmp.lmp.command('group vertex_edge id 3701')
trilmp.lmp.command('group vertex_edge id 3800')
trilmp.lmp.command('group vertex_edge id 3801')
trilmp.lmp.command('group vertex_edge id 3900')
trilmp.lmp.command('group vertex_edge id 3901')
trilmp.lmp.command('group vertex_edge id 4000')
trilmp.lmp.command('group vertex_edge id 4001')
trilmp.lmp.command('group vertex_edge id 4100')
trilmp.lmp.command('group vertex_edge id 4101')
trilmp.lmp.command('group vertex_edge id 4200')
trilmp.lmp.command('group vertex_edge id 4201')
trilmp.lmp.command('group vertex_edge id 4300')
trilmp.lmp.command('group vertex_edge id 4301')
trilmp.lmp.command('group vertex_edge id 4400')
trilmp.lmp.command('group vertex_edge id 4401')
trilmp.lmp.command('group vertex_edge id 4500')
trilmp.lmp.command('group vertex_edge id 4501')
trilmp.lmp.command('group vertex_edge id 4600')
trilmp.lmp.command('group vertex_edge id 4601')
trilmp.lmp.command('group vertex_edge id 4700')
trilmp.lmp.command('group vertex_edge id 4701')
trilmp.lmp.command('group vertex_edge id 4800')
trilmp.lmp.command('group vertex_edge id 4801')
trilmp.lmp.command('group vertex_edge id 4900')
trilmp.lmp.command('group vertex_edge id 4901')
trilmp.lmp.command('group vertex_edge id 5000')
trilmp.lmp.command('group vertex_edge id 5001')
trilmp.lmp.command('group vertex_edge id 5100')
trilmp.lmp.command('group vertex_edge id 5101')
trilmp.lmp.command('group vertex_edge id 5200')
trilmp.lmp.command('group vertex_edge id 5201')
trilmp.lmp.command('group vertex_edge id 5300')
trilmp.lmp.command('group vertex_edge id 5301')
trilmp.lmp.command('group vertex_edge id 5400')
trilmp.lmp.command('group vertex_edge id 5401')
trilmp.lmp.command('group vertex_edge id 5500')
trilmp.lmp.command('group vertex_edge id 5501')
trilmp.lmp.command('group vertex_edge id 5600')
trilmp.lmp.command('group vertex_edge id 5601')
trilmp.lmp.command('group vertex_edge id 5700')
trilmp.lmp.command('group vertex_edge id 5701')
trilmp.lmp.command('group vertex_edge id 5800')
trilmp.lmp.command('group vertex_edge id 5801')
trilmp.lmp.command('group vertex_edge id 5900')
trilmp.lmp.command('group vertex_edge id 5901')
trilmp.lmp.command('group vertex_edge id 6000')
trilmp.lmp.command('group vertex_edge id 6001')
trilmp.lmp.command('group vertex_edge id 6100')
trilmp.lmp.command('group vertex_edge id 6101')
trilmp.lmp.command('group vertex_edge id 6200')
trilmp.lmp.command('group vertex_edge id 6201')
trilmp.lmp.command('group vertex_edge id 6300')
trilmp.lmp.command('group vertex_edge id 6301')
trilmp.lmp.command('group vertex_edge id 6400')
trilmp.lmp.command('group vertex_edge id 6401')
trilmp.lmp.command('group vertex_edge id 6500')
trilmp.lmp.command('group vertex_edge id 6501')
trilmp.lmp.command('group vertex_edge id 6600')
trilmp.lmp.command('group vertex_edge id 6601')
trilmp.lmp.command('group vertex_edge id 6700')
trilmp.lmp.command('group vertex_edge id 6701')
trilmp.lmp.command('group vertex_edge id 6800')
trilmp.lmp.command('group vertex_edge id 6801')
trilmp.lmp.command('group vertex_edge id 6900')
trilmp.lmp.command('group vertex_edge id 6901')
trilmp.lmp.command('group vertex_edge id 7000')
trilmp.lmp.command('group vertex_edge id 7001')
trilmp.lmp.command('group vertex_edge id 7100')
trilmp.lmp.command('group vertex_edge id 7101')
trilmp.lmp.command('group vertex_edge id 7200')
trilmp.lmp.command('group vertex_edge id 7201')
trilmp.lmp.command('group vertex_edge id 7300')
trilmp.lmp.command('group vertex_edge id 7301')
trilmp.lmp.command('group vertex_edge id 7400')
trilmp.lmp.command('group vertex_edge id 7401')
trilmp.lmp.command('group vertex_edge id 7500')
trilmp.lmp.command('group vertex_edge id 7501')
trilmp.lmp.command('group vertex_edge id 7600')
trilmp.lmp.command('group vertex_edge id 7601')
trilmp.lmp.command('group vertex_edge id 7700')
trilmp.lmp.command('group vertex_edge id 7701')
trilmp.lmp.command('group vertex_edge id 7800')
trilmp.lmp.command('group vertex_edge id 7801')
trilmp.lmp.command('group vertex_edge id 7900')
trilmp.lmp.command('group vertex_edge id 7901')
trilmp.lmp.command('group vertex_edge id 8000')
trilmp.lmp.command('group vertex_edge id 8001')
trilmp.lmp.command('group vertex_edge id 8100')
trilmp.lmp.command('group vertex_edge id 8101')
trilmp.lmp.command('group vertex_edge id 8200')
trilmp.lmp.command('group vertex_edge id 8201')
trilmp.lmp.command('group vertex_edge id 8300')
trilmp.lmp.command('group vertex_edge id 8301')
trilmp.lmp.command('group vertex_edge id 8400')
trilmp.lmp.command('group vertex_edge id 8401')
trilmp.lmp.command('group vertex_edge id 8500')
trilmp.lmp.command('group vertex_edge id 8501')
trilmp.lmp.command('group vertex_edge id 8600')
trilmp.lmp.command('group vertex_edge id 8601')
trilmp.lmp.command('group vertex_edge id 8700')
trilmp.lmp.command('group vertex_edge id 8701')
trilmp.lmp.command('group vertex_edge id 8800')
trilmp.lmp.command('group vertex_edge id 8801')
trilmp.lmp.command('group vertex_edge id 8900')
trilmp.lmp.command('group vertex_edge id 8901')
trilmp.lmp.command('group vertex_edge id 9000')
trilmp.lmp.command('group vertex_edge id 9001')
trilmp.lmp.command('group vertex_edge id 9100')
trilmp.lmp.command('group vertex_edge id 9101')
trilmp.lmp.command('group vertex_edge id 9200')
trilmp.lmp.command('group vertex_edge id 9201')
trilmp.lmp.command('group vertex_edge id 9300')
trilmp.lmp.command('group vertex_edge id 9301')
trilmp.lmp.command('group vertex_edge id 9400')
trilmp.lmp.command('group vertex_edge id 9401')
trilmp.lmp.command('group vertex_edge id 9500')
trilmp.lmp.command('group vertex_edge id 9501')
trilmp.lmp.command('group vertex_edge id 9600')
trilmp.lmp.command('group vertex_edge id 9601')
trilmp.lmp.command('group vertex_edge id 9700')
trilmp.lmp.command('group vertex_edge id 9701')
trilmp.lmp.command('group vertex_edge id 9800')
trilmp.lmp.command('group vertex_edge id 9801')
trilmp.lmp.command('group vertex_edge id 9900')
trilmp.lmp.command('group vertex_edge id 9901')
trilmp.lmp.command('group vertex_edge id 10000')
trilmp.lmp.command('group vertex_edge id 10001')
trilmp.lmp.command('group vertex_edge id 10100')
trilmp.lmp.command('group vertex_edge id 10101')
trilmp.lmp.command('group vertex_edge id 10200')
trilmp.lmp.command('group vertex_edge id 10201')
trilmp.lmp.command('group vertex_edge id 10300')
trilmp.lmp.command('group vertex_edge id 10301')
trilmp.lmp.command('group vertex_edge id 10400')
trilmp.lmp.command('group vertex_edge id 10401')
trilmp.lmp.command('group vertex_edge id 10500')
trilmp.lmp.command('group vertex_edge id 10501')
trilmp.lmp.command('group vertex_edge id 10600')
trilmp.lmp.command('group vertex_edge id 10601')
trilmp.lmp.command('group vertex_edge id 10700')
trilmp.lmp.command('group vertex_edge id 10701')
trilmp.lmp.command('group vertex_edge id 10800')
trilmp.lmp.command('group vertex_edge id 10801')
trilmp.lmp.command('group vertex_edge id 10900')
trilmp.lmp.command('group vertex_edge id 10901')
trilmp.lmp.command('group vertex_edge id 11000')
trilmp.lmp.command('group vertex_edge id 11001')
trilmp.lmp.command('group vertex_edge id 11100')
trilmp.lmp.command('group vertex_edge id 11101')
trilmp.lmp.command('group vertex_edge id 11200')
trilmp.lmp.command('group vertex_edge id 11201')
trilmp.lmp.command('group vertex_edge id 11300')
trilmp.lmp.command('group vertex_edge id 11301')
trilmp.lmp.command('group vertex_edge id 11400')
trilmp.lmp.command('group vertex_edge id 11401')
trilmp.lmp.command('group vertex_edge id 11402')
trilmp.lmp.command('group vertex_edge id 11403')
trilmp.lmp.command('group vertex_edge id 11404')
trilmp.lmp.command('group vertex_edge id 11405')
trilmp.lmp.command('group vertex_edge id 11406')
trilmp.lmp.command('group vertex_edge id 11407')
trilmp.lmp.command('group vertex_edge id 11408')
trilmp.lmp.command('group vertex_edge id 11409')
trilmp.lmp.command('group vertex_edge id 11410')
trilmp.lmp.command('group vertex_edge id 11411')
trilmp.lmp.command('group vertex_edge id 11412')
trilmp.lmp.command('group vertex_edge id 11413')
trilmp.lmp.command('group vertex_edge id 11414')
trilmp.lmp.command('group vertex_edge id 11415')
trilmp.lmp.command('group vertex_edge id 11416')
trilmp.lmp.command('group vertex_edge id 11417')
trilmp.lmp.command('group vertex_edge id 11418')
trilmp.lmp.command('group vertex_edge id 11419')
trilmp.lmp.command('group vertex_edge id 11420')
trilmp.lmp.command('group vertex_edge id 11421')
trilmp.lmp.command('group vertex_edge id 11422')
trilmp.lmp.command('group vertex_edge id 11423')
trilmp.lmp.command('group vertex_edge id 11424')
trilmp.lmp.command('group vertex_edge id 11425')
trilmp.lmp.command('group vertex_edge id 11426')
trilmp.lmp.command('group vertex_edge id 11427')
trilmp.lmp.command('group vertex_edge id 11428')
trilmp.lmp.command('group vertex_edge id 11429')
trilmp.lmp.command('group vertex_edge id 11430')
trilmp.lmp.command('group vertex_edge id 11431')
trilmp.lmp.command('group vertex_edge id 11432')
trilmp.lmp.command('group vertex_edge id 11433')
trilmp.lmp.command('group vertex_edge id 11434')
trilmp.lmp.command('group vertex_edge id 11435')
trilmp.lmp.command('group vertex_edge id 11436')
trilmp.lmp.command('group vertex_edge id 11437')
trilmp.lmp.command('group vertex_edge id 11438')
trilmp.lmp.command('group vertex_edge id 11439')
trilmp.lmp.command('group vertex_edge id 11440')
trilmp.lmp.command('group vertex_edge id 11441')
trilmp.lmp.command('group vertex_edge id 11442')
trilmp.lmp.command('group vertex_edge id 11443')
trilmp.lmp.command('group vertex_edge id 11444')
trilmp.lmp.command('group vertex_edge id 11445')
trilmp.lmp.command('group vertex_edge id 11446')
trilmp.lmp.command('group vertex_edge id 11447')
trilmp.lmp.command('group vertex_edge id 11448')
trilmp.lmp.command('group vertex_edge id 11449')
trilmp.lmp.command('group vertex_edge id 11450')
trilmp.lmp.command('group vertex_edge id 11451')
trilmp.lmp.command('group vertex_edge id 11452')
trilmp.lmp.command('group vertex_edge id 11453')
trilmp.lmp.command('group vertex_edge id 11454')
trilmp.lmp.command('group vertex_edge id 11455')
trilmp.lmp.command('group vertex_edge id 11456')
trilmp.lmp.command('group vertex_edge id 11457')
trilmp.lmp.command('group vertex_edge id 11458')
trilmp.lmp.command('group vertex_edge id 11459')
trilmp.lmp.command('group vertex_edge id 11460')
trilmp.lmp.command('group vertex_edge id 11461')
trilmp.lmp.command('group vertex_edge id 11462')
trilmp.lmp.command('group vertex_edge id 11463')
trilmp.lmp.command('group vertex_edge id 11464')
trilmp.lmp.command('group vertex_edge id 11465')
trilmp.lmp.command('group vertex_edge id 11466')
trilmp.lmp.command('group vertex_edge id 11467')
trilmp.lmp.command('group vertex_edge id 11468')
trilmp.lmp.command('group vertex_edge id 11469')
trilmp.lmp.command('group vertex_edge id 11470')
trilmp.lmp.command('group vertex_edge id 11471')
trilmp.lmp.command('group vertex_edge id 11472')
trilmp.lmp.command('group vertex_edge id 11473')
trilmp.lmp.command('group vertex_edge id 11474')
trilmp.lmp.command('group vertex_edge id 11475')
trilmp.lmp.command('group vertex_edge id 11476')
trilmp.lmp.command('group vertex_edge id 11477')
trilmp.lmp.command('group vertex_edge id 11478')
trilmp.lmp.command('group vertex_edge id 11479')
trilmp.lmp.command('group vertex_edge id 11480')
trilmp.lmp.command('group vertex_edge id 11481')
trilmp.lmp.command('group vertex_edge id 11482')
trilmp.lmp.command('group vertex_edge id 11483')
trilmp.lmp.command('group vertex_edge id 11484')
trilmp.lmp.command('group vertex_edge id 11485')
trilmp.lmp.command('group vertex_edge id 11486')
trilmp.lmp.command('group vertex_edge id 11487')
trilmp.lmp.command('group vertex_edge id 11488')
trilmp.lmp.command('group vertex_edge id 11489')
trilmp.lmp.command('group vertex_edge id 11490')
trilmp.lmp.command('group vertex_edge id 11491')
trilmp.lmp.command('group vertex_edge id 11492')
trilmp.lmp.command('group vertex_edge id 11493')
trilmp.lmp.command('group vertex_edge id 11494')
trilmp.lmp.command('group vertex_edge id 11495')
trilmp.lmp.command('group vertex_edge id 11496')
trilmp.lmp.command('group vertex_edge id 11497')
trilmp.lmp.command('group vertex_edge id 11498')
trilmp.lmp.command('group vertex_edge id 11499')
trilmp.lmp.command('group vertex_edge id 11500')

# BULK GROUP: All the 'vertices' not at the edge
trilmp.lmp.command("group BULK subtract vertices vertex_edge")

#                         -- END ---

# ..............................................................
# PAIR_STYLES
# [NECESSARY] 
# ..............................................................
# Here we define the LAMMPS interaction pair styles.
# The 'table' pair_styles MUST be one of them: this is how
# we enforce surface repulsion to prevent the self-intersection
# of the mesh. Rather than using a potential, we use a table
# to increase computational speed.
# ..............................................................

#                        --- START ---

# clean-up pair style as a safety measure 
# (TriLMP pre-defines some pair_styles when the object is created)
trilmp.lmp.command("pair_style none")

# Define the pair interactions for LAMMPS
# - The 'table' pair_style is compulsory (this is how TriLMP makes sure
#   that the mesh does not self-intersect
# - The 'harmonic/cut' pair_style prevents overlaps between membrane and bead
# - You can add more pair_styles by appending after 'harmonic/cut'
trilmp.lmp.command(f"pair_style hybrid/overlay table linear 2000 harmonic/cut")

# [Compulsory lines - do not remove] To avoid mesh self-intersection
trilmp.lmp.command("pair_modify pair table special lj/coul 0.0 0.0 0.0 tail no")
trilmp.lmp.command("pair_coeff 1 1 table trimem_srp.table trimem_srp")

# Set all interactions to zero just in case for added potentials
trilmp.lmp.command("pair_coeff * * harmonic/cut 0 0")

# Introduce the parameters of the harmonic/cut pair style
trilmp.lmp.command("pair_coeff 1 2 harmonic/cut 1 6.5")

# Increase communication cutoff to avoid LAMMPS warnings
trilmp.lmp.command(f"comm_modify cutoff 11.0")

#                         -- END ---

# ..............................................................
# ADD A TETHERED BEAD
# [NECESSARY] 
# ..............................................................
# In this section we define several commands associated with the
# pulling of the membrane by using a tethered bead. Those
# properties are:
# - What the bond that connects the bead to the membrane looks like
# - What membrane vertex we tether to the pulling bead
# - The (optional) calculation of the force that the membrane exerts
#   on the bead
# ..............................................................

#                        --- START ---

# Define the properties of the bond that tethers the bead to the membrane:
# We choose a harmonic spring with a very soft spring constant (k = 1 k_BT/sigma^2)
# so that the pulling deformation we impose on the membrane is very soft
# with the goal of mimicking optical-tweezer pulling experiments.
# The rest length of the bond if 6.05, as the diameter of the pulling bead
# is 10 and the 'diameter' of a vertex in the network is 1 (the smallest
# distance between two vertices before an energy penalty arises)
trilmp.lmp.command("bond_style hybrid zero nocoeff harmonic")
trilmp.lmp.command("bond_coeff 1 zero 0.0")
trilmp.lmp.command(f"bond_coeff 2 harmonic 0.5 6.050000000000001")

# Tether one of the central vertices in the membrane to the pulling bead
trilmp.lmp.command(f"create_bonds single/bond 2 5750 11501")
# Group the two connected particles (we assign them to the 'paired' group)
trilmp.lmp.command("group paired id 5750 11501")

# Compute the energy stored in the bond and the force the membrane exerts on the bead
trilmp.lmp.command("compute MutualForce paired bond/local dist dx dy dz engpot force fx fy fz")
trilmp.lmp.command(f"dump  aveForce all local 5000 mutual_force.dump index c_MutualForce[*]")

#                         -- END ---

# ..............................................................
# DUMPS, COMPUTES, FIXES...
# [OPTIONAL] 
# ..............................................................
# This section includes trajectory dumps, computes and other
# LAMMPS fixes that are useful to monitor the status of the
# simulation. This is totally optional for our purposes
# ..............................................................

#                        --- START ---

# Dump the simulation trajectory for visualization
trilmp.lmp.command(f"dump XYZ all custom 10000 trajectory.dump id type x y z")

# Dump out the network connectivity for visualization
trilmp.lmp.command("compute MEMBONDS vertices property/local batom1 batom2")
trilmp.lmp.command(f"dump DMEMBONDS vertices local 5000 mem.bonds index c_MEMBONDS[1] c_MEMBONDS[2]")
trilmp.lmp.command("dump_modify DMEMBONDS format line '%d %0.0f %0.0f'")

# Compute the position of the CM of the membrane patch
trilmp.lmp.command("compute MembraneCOM vertices com")

# Compute the temperature of the membrane patch bulk
trilmp.lmp.command("compute TempComputeMem BULK temp")

# Print out the computations above
trilmp.lmp.command(f"fix aveMEM all ave/time 5000 1 5000 c_TempComputeMem c_MembraneCOM[1] c_MembraneCOM[2] c_MembraneCOM[3] file 'membrane_properties.dat'")

#                         -- END ---

# ..............................................................
# TIME INTEGRATION AND THERMOSTATTING
# [NECESSARY] 
# ..............................................................
# In this section we define the time-integration scheme
# as well as the thermostat that will act on the membrane.
# Notice how neither the particles at the edge of the membrane
# patch, nor the bead are integrated here.
# ..............................................................

#                        --- START ---

# Integrate the equations of motion of particles in the membrane BULK
trilmp.lmp.command("fix NVEMEM BULK nve")
# Apply a thermostat to particles in the membrane BULK
trilmp.lmp.command(f"fix LGVMEM BULK langevin 1.0 1.0 1.0 123 zero yes")

#                         -- END ---

# POST-POST-POST-POST-POST-POST-POST-POST-POST-POST-POST-POST
# Post-equilibration commands
# POST-POST-POST-POST-POST-POST-POST-POST-POST-POST-POST-POST

# Post-equilibration commands are passed to TriLmp by
# appending them to the list below
postequilibration_commands = []

# Once the membrane has equilibrated, we can start moving the bead
# upwards with a prescribed velocity
postequilibration_commands.append("fix MOTIONBEAD bead move linear 0 0 0.1")


# 333333333333333333333333333333333333333333333333333333333333333333333333333333
# 3. RUN THE SIMULATION using the .run(args) method
trilmp.run(100000000, postequilibration_lammps_commands = postequilibration_commands)

print("")
print("*** End of the simulation ***")
print("")

        
