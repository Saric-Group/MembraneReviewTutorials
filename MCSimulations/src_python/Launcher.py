import numpy as np
import os, glob
import ProduceInitialConfig

class MCSimulation:
    def __init__(self,kappa=20.0, surftension=7.0, network_columns = 0,
    N = 0, Ncoll = 0, Nedges =0, 
    Nbound =0, sigma_colloid = 0, 
    seed = 123, is_mem_fluid = 1,  
    start_from_config = 0, speed_pulling = 0.001, 
    position_bead = 0, k_harmonic = 1.0, 
    r_eq = 0.5, LX=100, LY = 100, LZ = 100):

        print("MC Simulation class initialized")

    def write_in_ves(self):

        file_content=f"""{self.kappa} {self.surftension} {self.N} {self.Ncoll} {self.Nedges} {self.Nbound} {self.sigma_colloid} {self.seed} {self.is_mem_fluid} {self.start_from_config} {self.speed_pulling} {self.position_bead} {self.k_harmonic} {self.r_eq} {self.LX} {self.LY} {self.LZ} 
kappa
surftension
N
Ncoll
Nedges
Nbound
sigma_colloid
seed
is_mem_fluid
start_from_config
speed_pulling
position_bead
k_harmonic
r_eq
LX 
LY 
LZ    
"""

        # write in the file
        with open(self.dirname+"/in.ves", 'w') as f:
            f.write(file_content)

def create_directory_structure(exp_codename, name, CLASS, varied_parameters, launch_directories):

    """
    Creates directories that contain computer experiments + saves them
    in a datafile so that it is easier to launch them using bash.
    """
    
    def include_attribute(obj, written, varied_parameters, name):
        index = 0
        for v in varied_parameters:
            if(hasattr(obj, v)) and written[index] == False:
                temp_name = v.replace("_", "")
                name += "_"+temp_name+"_{}".format(getattr(obj, v))
                written[index] = True
            index+=1

        return written, name
    counter = 0

    written = np.full(len(varied_parameters), False)
    written, name = include_attribute(CLASS, written, varied_parameters, name)
 
    new_path = name
    isExist = os.path.exists(new_path)
    if not isExist:
       os.makedirs(new_path)
    
    bash_path = new_path.replace(exp_codename+"/", "")
    launch_directories.write(bash_path)
    launch_directories.write('\n')
    return new_path

def prepare_bash_launch(experiment_codename):

    """
    Considers how many simulations have been launched already and tells you
    what directory to pass to slurm script in order to launch new simulations.
    """
    
    # create director_* files for launching directory
    bash_path = experiment_codename+"/bashdirs"
    isExist = os.path.exists(bash_path)
    if not isExist:
       os.makedirs(bash_path)

    # create outputs directory
    output_path = experiment_codename+"/output"
    isExist = os.path.exists(output_path)
    if not isExist:
       os.makedirs(output_path)

    # check how many directory files there are
    counter_directories = 1
    for files in glob.glob(bash_path+'/directories*'):
        counter_directories += 1

    print("LAUNCH ---> bashdirs/directories_"+str(counter_directories))
    launch_directories = open(bash_path+'/directories_'+str(counter_directories)+'.dat', 'w')
    return launch_directories

if __name__=="__main__":

    # initialize simulation class
    MCSim = MCSimulation()

    # set the name of your executable
    EXENAME='EXEtutorial'
    path_executable = f'../src_c/{EXENAME}'

    # set the experiment codename 
    experiment_codename = "../example_simulation"

    # prepare launching information for Slurm
    launch_directories = prepare_bash_launch(experiment_codename)

    # parameters that we will want to change (edit this list as desired)
    varied_parameters = ["kappa", "surftension", "network_columns", "speed_pulling", "k_harmonic", "r_eq", "position_bead", "seed", "scale_lattice"]

    # set mechanical properties of the membrane
    # kappa is the prefactor of the bending energy in the Hamiltonian
    # surface tension is the prefactor that penalizes area increases in the Hamiltonian
    MCSim.kappa = 20.0
    MCSim.surftension = 1.0

    # generate flat patch configuration of desired extension
    MCSim.network_columns = 20
    first_seed            = 123
    radius                = 0.95 # will define the bulk/mobile vertices
    scale_lattice         = 1.3 # will dictate the initial average bond length
    displace_z            = 70  # will place the flat patch in the z = -displace_z plane for the simulation
    MCSim.N, MCSim.Nedges, MCSim.Nbound, MCSim.LX = ProduceInitialConfig.prepare_config_hexagonal_flatpatch(MCSim.network_columns, first_seed, displace_z , scale_lattice = scale_lattice, radius = radius)

    # set parameters regarding pulling bead
    MCSim.Ncoll = 1
    MCSim.sigma_colloid = 10.0
    MCSim.speed_pulling = 0.0001
    MCSim.k_harmonic = 1.0
    MCSim.r_eq = 5.5
    MCSim.scale_lattice = scale_lattice

    # simulation box properties
    MCSim.LX +=2
    MCSim.LX = int(MCSim.LX)
    MCSim.LY = int(MCSim.LX)
    MCSim.LZ = int(195)

    # basic running parameters
    MCSim.is_mem_fluid = 1
    MCSim.start_from_config = 0

    base_position = 50
    positions_bead = np.arange(base_position, base_position*10, base_position*10+1)
    
    # launch the simulations
    simulations_launched = 0

    # add as many replica as you want
    for seed in [111]:
        
        # update the seed
        MCSim.seed = seed 

        for position in positions_bead:

            # count how many simulations are being launched
            simulations_launched+=1

            # introduce the new position
            MCSim.position_bead = position 

            # generate the directory
            dirname = experiment_codename+"/runs/"
            path = create_directory_structure(experiment_codename, dirname, MCSim, varied_parameters, launch_directories)
            MCSim.dirname = path

            MCSim.write_in_ves()

            # write the name of the file
            filename = "flatpatch_N_"+str(MCSim.N)+"_Sigma_"+str(int(MCSim.sigma_colloid))+"_.dat"
            os.system(f'mv {filename} {path}')
            os.system(f'cp {path_executable} {path}')

    print("Total simulations launched: ", simulations_launched)




