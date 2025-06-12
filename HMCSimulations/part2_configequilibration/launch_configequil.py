# ---------------------------------------------------------------------#
# TUTORIAL to simulate a TriLMP membrane patch                         #
# PART 2: Equilibration of the deformation                             #
# Author: Maitane Muñoz-Basagoiti (maitane.munoz-basagoiti@ista.ac.at) #
#                                                                      #
# This code reloads a single TriLMP pickle checkpoint.                 #
# Nonetheless, even though we do not need to initialize the system     #
# anymore, the same program structure as in the previous test          #
# (initialize Trilmp object, add LAMMPS commands and run) applies.     #
#                                                                      #
# Such checkpoint corresponds to a membrane patch with fixed boundary  #
# conditions in the z = 0 plane and a bead at some z = z0 position     #
# tethered to the membrane. In this simulation, such bead will remain  #
# fixed. As a result of the pulling from the bead, the membrane will   #
# relax into a tether/tube configuration.                              #
# The network has 11500 vertices.                                      #
# ---------------------------------------------------------------------#

# .....................................................................#
#       PRELIMINARY STEPS NEEDED AND CHECKPOINT RELOADING
# .....................................................................#

# Importe the necessary modules
from pathlib import Path
import trimem.core as m
from trimem.mc.trilmp import *

# Generate the 'checkpoints' directory (as for part 1) 
# In this directory, we will save the (pickle) checkpoints
# to use in the second part of the simulation, i.e.,
# configurations from which we will let the system
# equilibrate
Path("checkpoints").mkdir(exist_ok=True)

#11111111111111111111111111111111111111111111111111111111111111
# Load the checkpoint
# By doing so, 'loaded_checkpoint' is a new Trilmp object
# that remembers the network vertex positions and faces
# at some step of the simulation in the part_1, aside from the
# properties corresponding to the pulling bead
# This stage here would be the equivalent to stage 1 in
# the previous simulation
checkpoint_directory='./'
checkpoint_name='ckpt_example.pickle'
loaded_checkpoint = read_checkpoint(checkpoint_directory+checkpoint_name)

# 222222222222222222222222222222222222222222222222222222222222
# .....................................................................#
# LAMMPS ADDITIONS AND MODIFICATIONS
# [NECESSARY] 
# .....................................................................#
# The TriLMP configuration saved and just loaded only keeps information
# about the network vertices and faces. Therefore, some of the 
# previously defined information in the form of LAMMPS commands has to
# be defined again. 
# (More technically, when the checkpoint is unpickled, the LAMMPS lmp
# object is reinitialised, and we lose the previous one)
# .....................................................................#

#                        --- START ---

# ..........................
# REFRESH GROUP INFORMATION
# ..........................

# Vertices at the edge of the patch (they will remain fixed)
loaded_checkpoint.lmp.command('group vertex_edge id 1')
loaded_checkpoint.lmp.command('group vertex_edge id 2')
loaded_checkpoint.lmp.command('group vertex_edge id 3')
loaded_checkpoint.lmp.command('group vertex_edge id 4')
loaded_checkpoint.lmp.command('group vertex_edge id 5')
loaded_checkpoint.lmp.command('group vertex_edge id 6')
loaded_checkpoint.lmp.command('group vertex_edge id 7')
loaded_checkpoint.lmp.command('group vertex_edge id 8')
loaded_checkpoint.lmp.command('group vertex_edge id 9')
loaded_checkpoint.lmp.command('group vertex_edge id 10')
loaded_checkpoint.lmp.command('group vertex_edge id 11')
loaded_checkpoint.lmp.command('group vertex_edge id 12')
loaded_checkpoint.lmp.command('group vertex_edge id 13')
loaded_checkpoint.lmp.command('group vertex_edge id 14')
loaded_checkpoint.lmp.command('group vertex_edge id 15')
loaded_checkpoint.lmp.command('group vertex_edge id 16')
loaded_checkpoint.lmp.command('group vertex_edge id 17')
loaded_checkpoint.lmp.command('group vertex_edge id 18')
loaded_checkpoint.lmp.command('group vertex_edge id 19')
loaded_checkpoint.lmp.command('group vertex_edge id 20')
loaded_checkpoint.lmp.command('group vertex_edge id 21')
loaded_checkpoint.lmp.command('group vertex_edge id 22')
loaded_checkpoint.lmp.command('group vertex_edge id 23')
loaded_checkpoint.lmp.command('group vertex_edge id 24')
loaded_checkpoint.lmp.command('group vertex_edge id 25')
loaded_checkpoint.lmp.command('group vertex_edge id 26')
loaded_checkpoint.lmp.command('group vertex_edge id 27')
loaded_checkpoint.lmp.command('group vertex_edge id 28')
loaded_checkpoint.lmp.command('group vertex_edge id 29')
loaded_checkpoint.lmp.command('group vertex_edge id 30')
loaded_checkpoint.lmp.command('group vertex_edge id 31')
loaded_checkpoint.lmp.command('group vertex_edge id 32')
loaded_checkpoint.lmp.command('group vertex_edge id 33')
loaded_checkpoint.lmp.command('group vertex_edge id 34')
loaded_checkpoint.lmp.command('group vertex_edge id 35')
loaded_checkpoint.lmp.command('group vertex_edge id 36')
loaded_checkpoint.lmp.command('group vertex_edge id 37')
loaded_checkpoint.lmp.command('group vertex_edge id 38')
loaded_checkpoint.lmp.command('group vertex_edge id 39')
loaded_checkpoint.lmp.command('group vertex_edge id 40')
loaded_checkpoint.lmp.command('group vertex_edge id 41')
loaded_checkpoint.lmp.command('group vertex_edge id 42')
loaded_checkpoint.lmp.command('group vertex_edge id 43')
loaded_checkpoint.lmp.command('group vertex_edge id 44')
loaded_checkpoint.lmp.command('group vertex_edge id 45')
loaded_checkpoint.lmp.command('group vertex_edge id 46')
loaded_checkpoint.lmp.command('group vertex_edge id 47')
loaded_checkpoint.lmp.command('group vertex_edge id 48')
loaded_checkpoint.lmp.command('group vertex_edge id 49')
loaded_checkpoint.lmp.command('group vertex_edge id 50')
loaded_checkpoint.lmp.command('group vertex_edge id 51')
loaded_checkpoint.lmp.command('group vertex_edge id 52')
loaded_checkpoint.lmp.command('group vertex_edge id 53')
loaded_checkpoint.lmp.command('group vertex_edge id 54')
loaded_checkpoint.lmp.command('group vertex_edge id 55')
loaded_checkpoint.lmp.command('group vertex_edge id 56')
loaded_checkpoint.lmp.command('group vertex_edge id 57')
loaded_checkpoint.lmp.command('group vertex_edge id 58')
loaded_checkpoint.lmp.command('group vertex_edge id 59')
loaded_checkpoint.lmp.command('group vertex_edge id 60')
loaded_checkpoint.lmp.command('group vertex_edge id 61')
loaded_checkpoint.lmp.command('group vertex_edge id 62')
loaded_checkpoint.lmp.command('group vertex_edge id 63')
loaded_checkpoint.lmp.command('group vertex_edge id 64')
loaded_checkpoint.lmp.command('group vertex_edge id 65')
loaded_checkpoint.lmp.command('group vertex_edge id 66')
loaded_checkpoint.lmp.command('group vertex_edge id 67')
loaded_checkpoint.lmp.command('group vertex_edge id 68')
loaded_checkpoint.lmp.command('group vertex_edge id 69')
loaded_checkpoint.lmp.command('group vertex_edge id 70')
loaded_checkpoint.lmp.command('group vertex_edge id 71')
loaded_checkpoint.lmp.command('group vertex_edge id 72')
loaded_checkpoint.lmp.command('group vertex_edge id 73')
loaded_checkpoint.lmp.command('group vertex_edge id 74')
loaded_checkpoint.lmp.command('group vertex_edge id 75')
loaded_checkpoint.lmp.command('group vertex_edge id 76')
loaded_checkpoint.lmp.command('group vertex_edge id 77')
loaded_checkpoint.lmp.command('group vertex_edge id 78')
loaded_checkpoint.lmp.command('group vertex_edge id 79')
loaded_checkpoint.lmp.command('group vertex_edge id 80')
loaded_checkpoint.lmp.command('group vertex_edge id 81')
loaded_checkpoint.lmp.command('group vertex_edge id 82')
loaded_checkpoint.lmp.command('group vertex_edge id 83')
loaded_checkpoint.lmp.command('group vertex_edge id 84')
loaded_checkpoint.lmp.command('group vertex_edge id 85')
loaded_checkpoint.lmp.command('group vertex_edge id 86')
loaded_checkpoint.lmp.command('group vertex_edge id 87')
loaded_checkpoint.lmp.command('group vertex_edge id 88')
loaded_checkpoint.lmp.command('group vertex_edge id 89')
loaded_checkpoint.lmp.command('group vertex_edge id 90')
loaded_checkpoint.lmp.command('group vertex_edge id 91')
loaded_checkpoint.lmp.command('group vertex_edge id 92')
loaded_checkpoint.lmp.command('group vertex_edge id 93')
loaded_checkpoint.lmp.command('group vertex_edge id 94')
loaded_checkpoint.lmp.command('group vertex_edge id 95')
loaded_checkpoint.lmp.command('group vertex_edge id 96')
loaded_checkpoint.lmp.command('group vertex_edge id 97')
loaded_checkpoint.lmp.command('group vertex_edge id 98')
loaded_checkpoint.lmp.command('group vertex_edge id 99')
loaded_checkpoint.lmp.command('group vertex_edge id 100')
loaded_checkpoint.lmp.command('group vertex_edge id 101')
loaded_checkpoint.lmp.command('group vertex_edge id 200')
loaded_checkpoint.lmp.command('group vertex_edge id 201')
loaded_checkpoint.lmp.command('group vertex_edge id 300')
loaded_checkpoint.lmp.command('group vertex_edge id 301')
loaded_checkpoint.lmp.command('group vertex_edge id 400')
loaded_checkpoint.lmp.command('group vertex_edge id 401')
loaded_checkpoint.lmp.command('group vertex_edge id 500')
loaded_checkpoint.lmp.command('group vertex_edge id 501')
loaded_checkpoint.lmp.command('group vertex_edge id 600')
loaded_checkpoint.lmp.command('group vertex_edge id 601')
loaded_checkpoint.lmp.command('group vertex_edge id 700')
loaded_checkpoint.lmp.command('group vertex_edge id 701')
loaded_checkpoint.lmp.command('group vertex_edge id 800')
loaded_checkpoint.lmp.command('group vertex_edge id 801')
loaded_checkpoint.lmp.command('group vertex_edge id 900')
loaded_checkpoint.lmp.command('group vertex_edge id 901')
loaded_checkpoint.lmp.command('group vertex_edge id 1000')
loaded_checkpoint.lmp.command('group vertex_edge id 1001')
loaded_checkpoint.lmp.command('group vertex_edge id 1100')
loaded_checkpoint.lmp.command('group vertex_edge id 1101')
loaded_checkpoint.lmp.command('group vertex_edge id 1200')
loaded_checkpoint.lmp.command('group vertex_edge id 1201')
loaded_checkpoint.lmp.command('group vertex_edge id 1300')
loaded_checkpoint.lmp.command('group vertex_edge id 1301')
loaded_checkpoint.lmp.command('group vertex_edge id 1400')
loaded_checkpoint.lmp.command('group vertex_edge id 1401')
loaded_checkpoint.lmp.command('group vertex_edge id 1500')
loaded_checkpoint.lmp.command('group vertex_edge id 1501')
loaded_checkpoint.lmp.command('group vertex_edge id 1600')
loaded_checkpoint.lmp.command('group vertex_edge id 1601')
loaded_checkpoint.lmp.command('group vertex_edge id 1700')
loaded_checkpoint.lmp.command('group vertex_edge id 1701')
loaded_checkpoint.lmp.command('group vertex_edge id 1800')
loaded_checkpoint.lmp.command('group vertex_edge id 1801')
loaded_checkpoint.lmp.command('group vertex_edge id 1900')
loaded_checkpoint.lmp.command('group vertex_edge id 1901')
loaded_checkpoint.lmp.command('group vertex_edge id 2000')
loaded_checkpoint.lmp.command('group vertex_edge id 2001')
loaded_checkpoint.lmp.command('group vertex_edge id 2100')
loaded_checkpoint.lmp.command('group vertex_edge id 2101')
loaded_checkpoint.lmp.command('group vertex_edge id 2200')
loaded_checkpoint.lmp.command('group vertex_edge id 2201')
loaded_checkpoint.lmp.command('group vertex_edge id 2300')
loaded_checkpoint.lmp.command('group vertex_edge id 2301')
loaded_checkpoint.lmp.command('group vertex_edge id 2400')
loaded_checkpoint.lmp.command('group vertex_edge id 2401')
loaded_checkpoint.lmp.command('group vertex_edge id 2500')
loaded_checkpoint.lmp.command('group vertex_edge id 2501')
loaded_checkpoint.lmp.command('group vertex_edge id 2600')
loaded_checkpoint.lmp.command('group vertex_edge id 2601')
loaded_checkpoint.lmp.command('group vertex_edge id 2700')
loaded_checkpoint.lmp.command('group vertex_edge id 2701')
loaded_checkpoint.lmp.command('group vertex_edge id 2800')
loaded_checkpoint.lmp.command('group vertex_edge id 2801')
loaded_checkpoint.lmp.command('group vertex_edge id 2900')
loaded_checkpoint.lmp.command('group vertex_edge id 2901')
loaded_checkpoint.lmp.command('group vertex_edge id 3000')
loaded_checkpoint.lmp.command('group vertex_edge id 3001')
loaded_checkpoint.lmp.command('group vertex_edge id 3100')
loaded_checkpoint.lmp.command('group vertex_edge id 3101')
loaded_checkpoint.lmp.command('group vertex_edge id 3200')
loaded_checkpoint.lmp.command('group vertex_edge id 3201')
loaded_checkpoint.lmp.command('group vertex_edge id 3300')
loaded_checkpoint.lmp.command('group vertex_edge id 3301')
loaded_checkpoint.lmp.command('group vertex_edge id 3400')
loaded_checkpoint.lmp.command('group vertex_edge id 3401')
loaded_checkpoint.lmp.command('group vertex_edge id 3500')
loaded_checkpoint.lmp.command('group vertex_edge id 3501')
loaded_checkpoint.lmp.command('group vertex_edge id 3600')
loaded_checkpoint.lmp.command('group vertex_edge id 3601')
loaded_checkpoint.lmp.command('group vertex_edge id 3700')
loaded_checkpoint.lmp.command('group vertex_edge id 3701')
loaded_checkpoint.lmp.command('group vertex_edge id 3800')
loaded_checkpoint.lmp.command('group vertex_edge id 3801')
loaded_checkpoint.lmp.command('group vertex_edge id 3900')
loaded_checkpoint.lmp.command('group vertex_edge id 3901')
loaded_checkpoint.lmp.command('group vertex_edge id 4000')
loaded_checkpoint.lmp.command('group vertex_edge id 4001')
loaded_checkpoint.lmp.command('group vertex_edge id 4100')
loaded_checkpoint.lmp.command('group vertex_edge id 4101')
loaded_checkpoint.lmp.command('group vertex_edge id 4200')
loaded_checkpoint.lmp.command('group vertex_edge id 4201')
loaded_checkpoint.lmp.command('group vertex_edge id 4300')
loaded_checkpoint.lmp.command('group vertex_edge id 4301')
loaded_checkpoint.lmp.command('group vertex_edge id 4400')
loaded_checkpoint.lmp.command('group vertex_edge id 4401')
loaded_checkpoint.lmp.command('group vertex_edge id 4500')
loaded_checkpoint.lmp.command('group vertex_edge id 4501')
loaded_checkpoint.lmp.command('group vertex_edge id 4600')
loaded_checkpoint.lmp.command('group vertex_edge id 4601')
loaded_checkpoint.lmp.command('group vertex_edge id 4700')
loaded_checkpoint.lmp.command('group vertex_edge id 4701')
loaded_checkpoint.lmp.command('group vertex_edge id 4800')
loaded_checkpoint.lmp.command('group vertex_edge id 4801')
loaded_checkpoint.lmp.command('group vertex_edge id 4900')
loaded_checkpoint.lmp.command('group vertex_edge id 4901')
loaded_checkpoint.lmp.command('group vertex_edge id 5000')
loaded_checkpoint.lmp.command('group vertex_edge id 5001')
loaded_checkpoint.lmp.command('group vertex_edge id 5100')
loaded_checkpoint.lmp.command('group vertex_edge id 5101')
loaded_checkpoint.lmp.command('group vertex_edge id 5200')
loaded_checkpoint.lmp.command('group vertex_edge id 5201')
loaded_checkpoint.lmp.command('group vertex_edge id 5300')
loaded_checkpoint.lmp.command('group vertex_edge id 5301')
loaded_checkpoint.lmp.command('group vertex_edge id 5400')
loaded_checkpoint.lmp.command('group vertex_edge id 5401')
loaded_checkpoint.lmp.command('group vertex_edge id 5500')
loaded_checkpoint.lmp.command('group vertex_edge id 5501')
loaded_checkpoint.lmp.command('group vertex_edge id 5600')
loaded_checkpoint.lmp.command('group vertex_edge id 5601')
loaded_checkpoint.lmp.command('group vertex_edge id 5700')
loaded_checkpoint.lmp.command('group vertex_edge id 5701')
loaded_checkpoint.lmp.command('group vertex_edge id 5800')
loaded_checkpoint.lmp.command('group vertex_edge id 5801')
loaded_checkpoint.lmp.command('group vertex_edge id 5900')
loaded_checkpoint.lmp.command('group vertex_edge id 5901')
loaded_checkpoint.lmp.command('group vertex_edge id 6000')
loaded_checkpoint.lmp.command('group vertex_edge id 6001')
loaded_checkpoint.lmp.command('group vertex_edge id 6100')
loaded_checkpoint.lmp.command('group vertex_edge id 6101')
loaded_checkpoint.lmp.command('group vertex_edge id 6200')
loaded_checkpoint.lmp.command('group vertex_edge id 6201')
loaded_checkpoint.lmp.command('group vertex_edge id 6300')
loaded_checkpoint.lmp.command('group vertex_edge id 6301')
loaded_checkpoint.lmp.command('group vertex_edge id 6400')
loaded_checkpoint.lmp.command('group vertex_edge id 6401')
loaded_checkpoint.lmp.command('group vertex_edge id 6500')
loaded_checkpoint.lmp.command('group vertex_edge id 6501')
loaded_checkpoint.lmp.command('group vertex_edge id 6600')
loaded_checkpoint.lmp.command('group vertex_edge id 6601')
loaded_checkpoint.lmp.command('group vertex_edge id 6700')
loaded_checkpoint.lmp.command('group vertex_edge id 6701')
loaded_checkpoint.lmp.command('group vertex_edge id 6800')
loaded_checkpoint.lmp.command('group vertex_edge id 6801')
loaded_checkpoint.lmp.command('group vertex_edge id 6900')
loaded_checkpoint.lmp.command('group vertex_edge id 6901')
loaded_checkpoint.lmp.command('group vertex_edge id 7000')
loaded_checkpoint.lmp.command('group vertex_edge id 7001')
loaded_checkpoint.lmp.command('group vertex_edge id 7100')
loaded_checkpoint.lmp.command('group vertex_edge id 7101')
loaded_checkpoint.lmp.command('group vertex_edge id 7200')
loaded_checkpoint.lmp.command('group vertex_edge id 7201')
loaded_checkpoint.lmp.command('group vertex_edge id 7300')
loaded_checkpoint.lmp.command('group vertex_edge id 7301')
loaded_checkpoint.lmp.command('group vertex_edge id 7400')
loaded_checkpoint.lmp.command('group vertex_edge id 7401')
loaded_checkpoint.lmp.command('group vertex_edge id 7500')
loaded_checkpoint.lmp.command('group vertex_edge id 7501')
loaded_checkpoint.lmp.command('group vertex_edge id 7600')
loaded_checkpoint.lmp.command('group vertex_edge id 7601')
loaded_checkpoint.lmp.command('group vertex_edge id 7700')
loaded_checkpoint.lmp.command('group vertex_edge id 7701')
loaded_checkpoint.lmp.command('group vertex_edge id 7800')
loaded_checkpoint.lmp.command('group vertex_edge id 7801')
loaded_checkpoint.lmp.command('group vertex_edge id 7900')
loaded_checkpoint.lmp.command('group vertex_edge id 7901')
loaded_checkpoint.lmp.command('group vertex_edge id 8000')
loaded_checkpoint.lmp.command('group vertex_edge id 8001')
loaded_checkpoint.lmp.command('group vertex_edge id 8100')
loaded_checkpoint.lmp.command('group vertex_edge id 8101')
loaded_checkpoint.lmp.command('group vertex_edge id 8200')
loaded_checkpoint.lmp.command('group vertex_edge id 8201')
loaded_checkpoint.lmp.command('group vertex_edge id 8300')
loaded_checkpoint.lmp.command('group vertex_edge id 8301')
loaded_checkpoint.lmp.command('group vertex_edge id 8400')
loaded_checkpoint.lmp.command('group vertex_edge id 8401')
loaded_checkpoint.lmp.command('group vertex_edge id 8500')
loaded_checkpoint.lmp.command('group vertex_edge id 8501')
loaded_checkpoint.lmp.command('group vertex_edge id 8600')
loaded_checkpoint.lmp.command('group vertex_edge id 8601')
loaded_checkpoint.lmp.command('group vertex_edge id 8700')
loaded_checkpoint.lmp.command('group vertex_edge id 8701')
loaded_checkpoint.lmp.command('group vertex_edge id 8800')
loaded_checkpoint.lmp.command('group vertex_edge id 8801')
loaded_checkpoint.lmp.command('group vertex_edge id 8900')
loaded_checkpoint.lmp.command('group vertex_edge id 8901')
loaded_checkpoint.lmp.command('group vertex_edge id 9000')
loaded_checkpoint.lmp.command('group vertex_edge id 9001')
loaded_checkpoint.lmp.command('group vertex_edge id 9100')
loaded_checkpoint.lmp.command('group vertex_edge id 9101')
loaded_checkpoint.lmp.command('group vertex_edge id 9200')
loaded_checkpoint.lmp.command('group vertex_edge id 9201')
loaded_checkpoint.lmp.command('group vertex_edge id 9300')
loaded_checkpoint.lmp.command('group vertex_edge id 9301')
loaded_checkpoint.lmp.command('group vertex_edge id 9400')
loaded_checkpoint.lmp.command('group vertex_edge id 9401')
loaded_checkpoint.lmp.command('group vertex_edge id 9500')
loaded_checkpoint.lmp.command('group vertex_edge id 9501')
loaded_checkpoint.lmp.command('group vertex_edge id 9600')
loaded_checkpoint.lmp.command('group vertex_edge id 9601')
loaded_checkpoint.lmp.command('group vertex_edge id 9700')
loaded_checkpoint.lmp.command('group vertex_edge id 9701')
loaded_checkpoint.lmp.command('group vertex_edge id 9800')
loaded_checkpoint.lmp.command('group vertex_edge id 9801')
loaded_checkpoint.lmp.command('group vertex_edge id 9900')
loaded_checkpoint.lmp.command('group vertex_edge id 9901')
loaded_checkpoint.lmp.command('group vertex_edge id 10000')
loaded_checkpoint.lmp.command('group vertex_edge id 10001')
loaded_checkpoint.lmp.command('group vertex_edge id 10100')
loaded_checkpoint.lmp.command('group vertex_edge id 10101')
loaded_checkpoint.lmp.command('group vertex_edge id 10200')
loaded_checkpoint.lmp.command('group vertex_edge id 10201')
loaded_checkpoint.lmp.command('group vertex_edge id 10300')
loaded_checkpoint.lmp.command('group vertex_edge id 10301')
loaded_checkpoint.lmp.command('group vertex_edge id 10400')
loaded_checkpoint.lmp.command('group vertex_edge id 10401')
loaded_checkpoint.lmp.command('group vertex_edge id 10500')
loaded_checkpoint.lmp.command('group vertex_edge id 10501')
loaded_checkpoint.lmp.command('group vertex_edge id 10600')
loaded_checkpoint.lmp.command('group vertex_edge id 10601')
loaded_checkpoint.lmp.command('group vertex_edge id 10700')
loaded_checkpoint.lmp.command('group vertex_edge id 10701')
loaded_checkpoint.lmp.command('group vertex_edge id 10800')
loaded_checkpoint.lmp.command('group vertex_edge id 10801')
loaded_checkpoint.lmp.command('group vertex_edge id 10900')
loaded_checkpoint.lmp.command('group vertex_edge id 10901')
loaded_checkpoint.lmp.command('group vertex_edge id 11000')
loaded_checkpoint.lmp.command('group vertex_edge id 11001')
loaded_checkpoint.lmp.command('group vertex_edge id 11100')
loaded_checkpoint.lmp.command('group vertex_edge id 11101')
loaded_checkpoint.lmp.command('group vertex_edge id 11200')
loaded_checkpoint.lmp.command('group vertex_edge id 11201')
loaded_checkpoint.lmp.command('group vertex_edge id 11300')
loaded_checkpoint.lmp.command('group vertex_edge id 11301')
loaded_checkpoint.lmp.command('group vertex_edge id 11400')
loaded_checkpoint.lmp.command('group vertex_edge id 11401')
loaded_checkpoint.lmp.command('group vertex_edge id 11402')
loaded_checkpoint.lmp.command('group vertex_edge id 11403')
loaded_checkpoint.lmp.command('group vertex_edge id 11404')
loaded_checkpoint.lmp.command('group vertex_edge id 11405')
loaded_checkpoint.lmp.command('group vertex_edge id 11406')
loaded_checkpoint.lmp.command('group vertex_edge id 11407')
loaded_checkpoint.lmp.command('group vertex_edge id 11408')
loaded_checkpoint.lmp.command('group vertex_edge id 11409')
loaded_checkpoint.lmp.command('group vertex_edge id 11410')
loaded_checkpoint.lmp.command('group vertex_edge id 11411')
loaded_checkpoint.lmp.command('group vertex_edge id 11412')
loaded_checkpoint.lmp.command('group vertex_edge id 11413')
loaded_checkpoint.lmp.command('group vertex_edge id 11414')
loaded_checkpoint.lmp.command('group vertex_edge id 11415')
loaded_checkpoint.lmp.command('group vertex_edge id 11416')
loaded_checkpoint.lmp.command('group vertex_edge id 11417')
loaded_checkpoint.lmp.command('group vertex_edge id 11418')
loaded_checkpoint.lmp.command('group vertex_edge id 11419')
loaded_checkpoint.lmp.command('group vertex_edge id 11420')
loaded_checkpoint.lmp.command('group vertex_edge id 11421')
loaded_checkpoint.lmp.command('group vertex_edge id 11422')
loaded_checkpoint.lmp.command('group vertex_edge id 11423')
loaded_checkpoint.lmp.command('group vertex_edge id 11424')
loaded_checkpoint.lmp.command('group vertex_edge id 11425')
loaded_checkpoint.lmp.command('group vertex_edge id 11426')
loaded_checkpoint.lmp.command('group vertex_edge id 11427')
loaded_checkpoint.lmp.command('group vertex_edge id 11428')
loaded_checkpoint.lmp.command('group vertex_edge id 11429')
loaded_checkpoint.lmp.command('group vertex_edge id 11430')
loaded_checkpoint.lmp.command('group vertex_edge id 11431')
loaded_checkpoint.lmp.command('group vertex_edge id 11432')
loaded_checkpoint.lmp.command('group vertex_edge id 11433')
loaded_checkpoint.lmp.command('group vertex_edge id 11434')
loaded_checkpoint.lmp.command('group vertex_edge id 11435')
loaded_checkpoint.lmp.command('group vertex_edge id 11436')
loaded_checkpoint.lmp.command('group vertex_edge id 11437')
loaded_checkpoint.lmp.command('group vertex_edge id 11438')
loaded_checkpoint.lmp.command('group vertex_edge id 11439')
loaded_checkpoint.lmp.command('group vertex_edge id 11440')
loaded_checkpoint.lmp.command('group vertex_edge id 11441')
loaded_checkpoint.lmp.command('group vertex_edge id 11442')
loaded_checkpoint.lmp.command('group vertex_edge id 11443')
loaded_checkpoint.lmp.command('group vertex_edge id 11444')
loaded_checkpoint.lmp.command('group vertex_edge id 11445')
loaded_checkpoint.lmp.command('group vertex_edge id 11446')
loaded_checkpoint.lmp.command('group vertex_edge id 11447')
loaded_checkpoint.lmp.command('group vertex_edge id 11448')
loaded_checkpoint.lmp.command('group vertex_edge id 11449')
loaded_checkpoint.lmp.command('group vertex_edge id 11450')
loaded_checkpoint.lmp.command('group vertex_edge id 11451')
loaded_checkpoint.lmp.command('group vertex_edge id 11452')
loaded_checkpoint.lmp.command('group vertex_edge id 11453')
loaded_checkpoint.lmp.command('group vertex_edge id 11454')
loaded_checkpoint.lmp.command('group vertex_edge id 11455')
loaded_checkpoint.lmp.command('group vertex_edge id 11456')
loaded_checkpoint.lmp.command('group vertex_edge id 11457')
loaded_checkpoint.lmp.command('group vertex_edge id 11458')
loaded_checkpoint.lmp.command('group vertex_edge id 11459')
loaded_checkpoint.lmp.command('group vertex_edge id 11460')
loaded_checkpoint.lmp.command('group vertex_edge id 11461')
loaded_checkpoint.lmp.command('group vertex_edge id 11462')
loaded_checkpoint.lmp.command('group vertex_edge id 11463')
loaded_checkpoint.lmp.command('group vertex_edge id 11464')
loaded_checkpoint.lmp.command('group vertex_edge id 11465')
loaded_checkpoint.lmp.command('group vertex_edge id 11466')
loaded_checkpoint.lmp.command('group vertex_edge id 11467')
loaded_checkpoint.lmp.command('group vertex_edge id 11468')
loaded_checkpoint.lmp.command('group vertex_edge id 11469')
loaded_checkpoint.lmp.command('group vertex_edge id 11470')
loaded_checkpoint.lmp.command('group vertex_edge id 11471')
loaded_checkpoint.lmp.command('group vertex_edge id 11472')
loaded_checkpoint.lmp.command('group vertex_edge id 11473')
loaded_checkpoint.lmp.command('group vertex_edge id 11474')
loaded_checkpoint.lmp.command('group vertex_edge id 11475')
loaded_checkpoint.lmp.command('group vertex_edge id 11476')
loaded_checkpoint.lmp.command('group vertex_edge id 11477')
loaded_checkpoint.lmp.command('group vertex_edge id 11478')
loaded_checkpoint.lmp.command('group vertex_edge id 11479')
loaded_checkpoint.lmp.command('group vertex_edge id 11480')
loaded_checkpoint.lmp.command('group vertex_edge id 11481')
loaded_checkpoint.lmp.command('group vertex_edge id 11482')
loaded_checkpoint.lmp.command('group vertex_edge id 11483')
loaded_checkpoint.lmp.command('group vertex_edge id 11484')
loaded_checkpoint.lmp.command('group vertex_edge id 11485')
loaded_checkpoint.lmp.command('group vertex_edge id 11486')
loaded_checkpoint.lmp.command('group vertex_edge id 11487')
loaded_checkpoint.lmp.command('group vertex_edge id 11488')
loaded_checkpoint.lmp.command('group vertex_edge id 11489')
loaded_checkpoint.lmp.command('group vertex_edge id 11490')
loaded_checkpoint.lmp.command('group vertex_edge id 11491')
loaded_checkpoint.lmp.command('group vertex_edge id 11492')
loaded_checkpoint.lmp.command('group vertex_edge id 11493')
loaded_checkpoint.lmp.command('group vertex_edge id 11494')
loaded_checkpoint.lmp.command('group vertex_edge id 11495')
loaded_checkpoint.lmp.command('group vertex_edge id 11496')
loaded_checkpoint.lmp.command('group vertex_edge id 11497')
loaded_checkpoint.lmp.command('group vertex_edge id 11498')
loaded_checkpoint.lmp.command('group vertex_edge id 11499')
loaded_checkpoint.lmp.command('group vertex_edge id 11500')

# BULK GROUP: All the 'vertices' not at the edge
loaded_checkpoint.lmp.command("group BULK subtract vertices vertex_edge")

# ...............................
# REFRESH PAIR STYLE INFORMATION
# ...............................

# clean-up pair style as a safety measure 
loaded_checkpoint.lmp.command("pair_style none")

# Define the pair interactions for LAMMPS
# - The 'table' pair_style is compulsory (this is how TriLMP makes sure
#   that the mesh does not self-intersect
# - The 'harmonic/cut' pair_style prevents overlaps between membrane and bead
# - You can add more pair_styles by appending after 'harmonic/cut'loaded_checkpoint.lmp.command(f"pair_style hybrid/overlay table linear 2000 harmonic/cut")
loaded_checkpoint.lmp.command(f"pair_style hybrid/overlay table linear 2000 harmonic/cut")

# [Compulsory lines - do not remove] To avoid mesh self-intersection
loaded_checkpoint.lmp.command("pair_modify pair table special lj/coul 0.0 0.0 0.0 tail no")
loaded_checkpoint.lmp.command("pair_coeff 1 1 table trimem_srp.table trimem_srp")

# Set all interactions to zero just in case for added potentials
loaded_checkpoint.lmp.command("pair_coeff * * harmonic/cut 0 0")

# Introduce the parameters of the harmonic/cut pair style
loaded_checkpoint.lmp.command("pair_coeff 1 2 harmonic/cut 1 6.5")

# Increase communication cutoff to avoid LAMMPS warnings
loaded_checkpoint.lmp.command(f"comm_modify cutoff 11.0")

# ..................................
# REFRESH TETHERED BEAD INFORMATION
# ..................................

# Define the properties of the bond that tethers the bead to the membrane:
# We choose a harmonic spring with a very soft spring constant (k = 1 k_BT/sigma^2)
# so that the pulling deformation we impose on the membrane is very soft
# with the goal of mimicking optical-tweezer pulling experiments.
# The rest length of the bond if 6.05, as the diameter of the pulling bead
# is 10 and the 'diameter' of a vertex in the network is 1 (the smallest
# distance between two vertices before an energy penalty arises)
loaded_checkpoint.lmp.command("bond_style hybrid zero nocoeff harmonic")
loaded_checkpoint.lmp.command("bond_coeff 1 zero 0.0")
loaded_checkpoint.lmp.command(f"bond_coeff 2 harmonic 0.5 6.050000000000001")

# Tether one of the central vertices in the membrane to the pulling bead
loaded_checkpoint.lmp.command(f"create_bonds single/bond 2 5750 11501")
# Group the two connected particles (we assign them to the 'paired' group)
loaded_checkpoint.lmp.command("group paired id 5750 11501")

# Compute the energy stored in the bond and the force the membrane exerts on the bead
# Here mutual_force.dump is the key file from which we will be able to extract the
# force that it takes to extrude a membrane tether
loaded_checkpoint.lmp.command("compute MutualForce paired bond/local dist dx dy dz engpot force fx fy fz")
loaded_checkpoint.lmp.command(f"dump  aveForce all local 50000 mutual_force.dump index c_MutualForce[*]")

# ..................................
# REFRESH DUMPS, COMPUTES, FIXES
# ..................................

# Dump the simulation trajectory for visualization
loaded_checkpoint.lmp.command(f"dump XYZ all custom 50000 trajectory.dump id x y z type")

# Dump out the network connectivity for visualization
# (Make sure that the frequency at which you print bonds agrees with the frequency
# at which you print configurations in the trajectory for a correct visualization)
loaded_checkpoint.lmp.command("compute MEMBONDS vertices property/local batom1 batom2")
loaded_checkpoint.lmp.command(f"dump DMEMBONDS vertices local 50000 mem.bonds index c_MEMBONDS[1] c_MEMBONDS[2]")
loaded_checkpoint.lmp.command("dump_modify DMEMBONDS format line '%d %0.0f %0.0f'")

# (Lazy) extraction of the bead position for the given checkpoint
# As it wont move, printing frequency is not important here
loaded_checkpoint.lmp.command(f"dump BEADPOS bead custom 10000000 beadpos.dump id x y z type")

# Extraction of the properties acting on the membrane vertex on which we are pulling 
# This vertex is first included in the 'PULLING' group and then its properties dumped
loaded_checkpoint.lmp.command(f"group PULLING id 5750")
loaded_checkpoint.lmp.command(f"dump PULLPOS PULLING custom 50000 properties_pulled.dump id x y z fx fy fz type")

# Compute the position of the CM of the membrane patch
loaded_checkpoint.lmp.command("compute MembraneCOM vertices com")

# Compute the temperature of the membrane patch bulk
loaded_checkpoint.lmp.command("compute TempComputeMem BULK temp")

# Print out the computations above
loaded_checkpoint.lmp.command(f"fix  aveMEM all ave/time 50000 1 50000 c_TempComputeMem file 'membrane_properties.dat'")

# ..................................
# REFRESH INTEGRATORS AND THERMOSTATS
# ..................................

# include the integrators (pre-equilibration)
loaded_checkpoint.lmp.command("fix NVEMEM BULK nve")
loaded_checkpoint.lmp.command(f"fix LGVMEM BULK langevin 1.0 1.0 1.0 123 zero yes")

# 333333333333333333333333333333333333333333333333333333333333333333333333333333
# RUN THE SIMULATION
# Note how on this occasion we do not do anything about the position of the
# bead in the system. It will remain fixed at the position z = z0 saved in the
# checkpoint
loaded_checkpoint.run(100000000)

print("")
print("*** End of the simulation ***")
print("")

        
