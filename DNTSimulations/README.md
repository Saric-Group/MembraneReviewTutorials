# Tether extrusion with a Dynamically Triangulated Network (DTN)

## 1. Generate initial configurations

Running a TriLMP simulation involves 3 stages:

1. First, one creates the TriLmp object. 
This step initializes the simulation set-up. Here you define...
       - The number of particles in the system, their positions, types and properties
               + By default, the particles in fluid membrane are called 'vertices'
               + The bead tethered to the membrane, that will pull on it, is called 'bead' here
       - The dimensions of the simulation box and its boundary conditions (periodic or not)
       - The mechanical properties of the membrane (i.e., the parameters in the TriMEM hamiltonian)
       - The properties of the Hybrid Monte Carlo integration scheme (MD stages and MC stages)
               + TriLMP will alternate randomly between MD (vertex motion) and MC (connectivity update) stages
               + You can decide how many steps to integrate during the MD stages
                 as well as what fraction of network edges to attempt flipping in the MC stages
                 (see Siggel et al., 2022, for a discussion on the limitations of this parameter)
       - The length of the equilibration stage
       - The frequency at which performance, system info is outputted and checkpoints are created by TriLMP
2. Second, one adds LAMMPS commands to the simulation
 Here you can choose to pass as many LAMMPS commands as desired to your simulation.
 Note that there are some compulsory commands (see pair_styles or time-integration sections below)
 which must always be passed to ensure that TriLMP runs correctly.
 The LAMMPS commands can be passed in two different ways:
 - Pre-equilibration commands
       + Pass them directly to the TriLmp object before the simulation has started running
         (this is useful to equilibrate your system beforehand)
 - Post-equilibration commands
       + Pass them by appending them to the postequilibration_commands list below
          (they will be run once the pre-equilibration stage has ended, in order)
3. Third, run the simulation by calling the TriLmp.run(args) method with its corresponding arguments

## 2. Equilibrate the membrane
