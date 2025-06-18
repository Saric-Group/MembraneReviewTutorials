# Tether extrusion with a Dynamically Triangulated Network I: Monte Carlo simulations with C code

In this tutorial, we will measure the force required to extrude a tube from a membrane patch using a dynamically triangulated mesh model and Monte Carlo (MC) simulations. By the end of the tutorial, you should be able to produce the force-elongation profile shown in the review. 

## Table of contents
1. [Description of the MC code](#1-description-of-the-mc-code)
    - [Editable variables](#a-editable-variables)
    - [Backbone of the program](#b-backbone-of-the-program)
    - [Specifics of tube pulling](#c-specifics-of-tube-pulling)
2. [Compiling the code](#2-compiling-the-code)
3. [Running the code](#3-running-the-code)
    - [Python codes to initialize the system](#a-python-codes-to-initialize-the-system)
    - [The in.ves file in detail](#b-the-inves-file-in-detail)
    - [Data produced](#c-data-produced)
4. [Reporting bugs, feedback and credit](#4-reporting-bugs-beedback-and-credit)

## 1. Description of the MC code
The MC code can be found in the [`src_c`](src_c) directory. The main files are `main.c`, `Functions.h`, `DataStructures.h` and `PreprocessorDeclarations.h`. The rest of the files are associated to [`uthash`](https://troydhanson.github.io/uthash/), which is a library for hash tables used in the code to efficiently identify bulk and edge vertices (see review for details on the simulation set-up). 

### A. Editable variables

The main variables of the simulation are provided in the `in.ves` file (see below). However, the `PreprocessorDeclarations.h` contains 3 global variables that you can edit to run your simulation. These are:
- `MCEQUILIBRATE`, which sets the number of MC sweeps during the equilibration stage
- `MCSTEPS`, which sets the total number of MC sweeps in the program
- `WRITE_CONF`, which sets the frequency at which the program outputs data


### B. Backbone of the program

To generate fluid membrane configurations, the program loops over two steps:

1. Random mesh vertex moves, which are attempted by the `MC_xyz` function
2. Random connectivity update moves, which are attempted by the `switch_bond` function

### C. Specifics of tube pulling

To extrude a membrane tube, we tether a bead via a harmonic bond to a single vertex in the mesh. Once the membrane has equilibrated after `MCEQUILIBRATE` steps, we gradually update the position of the pulling bead through increments `speed_pulling` until it reaches a final prescribed position; it will remain fixed in place once this happens.  By design, the pulling bead can only move along the z-direction. Additionally, to keep the membrane tension constant, the simulation method implements a specific set-up based on bulk and edge vertices; see the review for further details and discussion on the set-up.

## 2. Compiling the code
To compile the code,
1. Go to directory where `*.c` and `*.h` files are
2. Run `gcc -o EXENAME *.c` in the terminal

Choose any executable name that what you want, i.e., substitute `EXENAME`. Please note that to compile in a HPC cluster, you might have to run 

```gcc *.c -lm -o EXENAME``` 

to also link the math libraries.

## 3. Running the code

To run the code, you have to place the executable `EXENAME` in a directory, together with the `in.ves` and `flatpatch*.dat` files. Then, simply execute `./EXENAME` in the command line.

### A. Python codes to initialize the system

In [`src_python`](src_python) we provide python codes that help you establish your working directory as well as all the files needed for the simulation. All you need to do is edit the `Launcher.py` file and run it: it will produce a directory for the simulation, together with the `in.ves` and `flatpath*.dat` files that correspond to the simulations you want to run. Check the `Launcher.py` file to see the parameters that you can change. 

### B. The `in.ves` file in detail

The `in.ves` file contains the parameters needed to run the simulation. These are provided in the first line of the file in order, and they are:
 - `kappa`: prefactor of the discretized bending energy, units of $k_BT = 1$
 - `surftension`: prefactor for the area term in the Hamiltonian, units of $\sigma^{2} k_B T$, where $\sigma = 1$ is the simulation unit of length
 - `N`: total number of mesh vertices
 - `Ncoll`: total number of pulling beads, by default 1.
 - `Nedges`: number of edge mesh vertices (fixed vertices)
 - `Nbound`: ID of the mesh vertex which is tethered to the bead
 - `Sigma_colloid`: diameter of the pulling bead, hard coded to $10~\sigma$
 - `seed`: simulation seed for random number generation
 - `is_mem_fluid`: set to 1 for fluid membranes, 0 for elastic ones (no connectivity update)
 - `start_from_config`: flag dictating whether we start from an already existing configuration; in this tutorial, it is set to 0
 - `speed_pulling`: distance that the pulling bead travels each time its position is updated
 - `position_bead`: final position of the pulling bead (measured from reference membrane plane)
 - `k_harmonic`: elastic constant of the harmonic bond that tethers the mesh vertex to the pulling bead
 - `r_eq`: equilibrium position of the harmonic bond tethering the pulling bead to the mesh
 - `LX`: length of the simulation box in the x-direction
 - `LY`: length of the simulation box in the y-direction
 - `LZ`: length of the simulation box in the z-direction

 ### C. Data produced

 The folder [example_simulation](../example_simulation) contains the typical output of a simulation. The data therein corresponds to a small membrane patch ran for 550k MC sweeps (only equilibration stage, the pulling bead has not started moving). The relevant files are:
 - `out.dump` contains the coordinates of all particles in the system in time
 - `bond.dump` contains the connectivity of the mesh in time
 - `energy.dump` contains the energy of the system (by columns:   (1) MC sweep number, (2) bending modulus (3) surface tension (4) bending energy (5) membrane area (6) energy stored in the harmonic bond (7) change in bending energy (8) change in the energy stored in the harmonic bond (9) change in the area energy).
 - `conf*` file corresponds to the initial configuration of the system
 - `position_pulling_bead.dat` contains the coordinates of the pulling bead and the membrane vertex attached to it.
 - `setnb_init.dat` contains the initial mesh connectivity
 - `settriangles_init.dat` contains the initial mesh triangles

 Both `out.dump` and `bond.dump` can be fed to software like [OVITO](https://www.ovito.org/) for visualization.


## 4. Reporting bugs, beedback and credit

Please report any bugs. Feedback on how to improve this tutorial is welcome.

The C code shared in this tutorial has been edited and adapted for the tutorial by members of the [Šarić group](https://github.com/Saric-Group/) during 2023-2025. The code was inherited from previous members of the group, and therefore some of the functions in the C code had already been developed. Please contact us if you worked in some of the functions presented in this code in the past and would like to be credited.