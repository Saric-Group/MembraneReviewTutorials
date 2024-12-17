Instructions for running a simulation of a membrane tube equilibration using the YLZ potential
==============================================================================================

Installation of LAMMPS with the required in-house changes
---------------------------------------------------------
- Download the LAMMPS source code (tested with version 29Aug2024)
- Replace files `compute_pressure.h` and `compute_pressure.cpp` in the src directory
- Compile and install LAMMPS (as outlined in LAMMPS documentation)
    - Installing with an MPI library for parallelisation is highly recommended

Installation of the required Python libraries
---------------------------------------------
- You need a python3 distribution with the following libraries
    -   numpy
    -   scipy
    -   ovito

Running of the LAMMPS simulation
--------------------------------
- Make sure that you have the `lmpinp*` and `lmpdat*` files in the same directory
- Run the simulation as `lmp -in lmpinp*`
- The ouput files are:
    -   `log.lammps` = standard LAMMPS log file
    -   `lmpout_xyz*.dat` = the trajectory, i.e. positions of particles
    -   `lmpout_boxsize*.dat` = size of the box
    -   `lmpout_Pprime.dat` = the value of P' as outlined in the article

Extracting the radii of the tubes from the simulation
-----------------------------------------------------
- Run `ovito_find_radius.py lmpout_xyz*`

Extracting the average radius for the given simulation
------------------------------------------------------
- Run `get_equilibrium_radius.py radii*`
- Put in the first and the last timestep that you wish to consider (For more detail on how to choose these rigorously, see the Cooke model tutorial)
