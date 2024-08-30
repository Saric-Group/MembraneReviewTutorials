# Tether extrusion with a Dynamically Triangulated Network

Generating the force-extrusion profile shown in OurReviewEtAl., (2024) requires two steps in the analysis, which are detailed below. All the codes necessary to generate the data are provided. These have been tested and thoroughly commented with the goal of teaching the simulator that follows the tutorial the basic functioning of the [TriLMP software](https://github.com/Saric-Group/trimem_sbeady). At the end of the tutorial, the simulator will be able to simulate Dynamically Triangulated Networks (DTNs) using TriLMP.

## 1. Generate initial configurations

To measure the force required to extrude a membrane tether, we will first run a single, long simulation. In this simulation, the membrane is initialized as a flat patch with $L \approx 100\sigma$ with fixed boundary conditions in the $z =0$ plane. Additionally, we attach a bead to one of the vertices in our DTN using a harmonic bond. We then move the bead with constant velocity along $z>0$ while membrane dynamics are integrated, which gives rise to a membrane tether. This simulation procedure generates a series of system checkpoints or configurations that will be used in the second stage of our analysis, where for a fixed bead position, we let the membrane relax.

### 1.1. Running the simulation
To run the simulation, go to the [`part1_generatetrajectory`](https://github.com/Saric-Group/MembraneReviewTutorials/tree/main/DNTSimulations/part1_generatetrajectory) directory. You can run the simulation directly by using:

```python launch_trajectory.py```

We additionally provide a bash script that allows you to submit the simulation into an HPC cluster.

### 1.2. A first contact with TriLMP

This first simulation allows you to understand the basic structure of TriLMP. Indeed, running a TriLMP simulation involves 3 stages, which we detail below.

**A. Create the TriLmp object**
This step initializes the simulation set-up. 

Here you define:
- The number of particles in the system, their positions, types and properties
  - By default, the particles in fluid membrane are called 'vertices'
  - The bead tethered to the membrane, that will pull on it, is called 'bead' here
- The dimensions of the simulation box and its boundary conditions (periodic or not)
- The mechanical properties of the membrane (i.e., the parameters in the TriMEM hamiltonian)
- The properties of the Hybrid Monte Carlo integration scheme (MD stages and MC stages)
  - TriLMP will alternate randomly between MD (vertex motion) and MC (connectivity update) stages
- The length of the equilibration stage
- The frequency at which performance, system info is outputted and checkpoints are created by TriLMP
  
**B. Add LAMMPS commands to the simulation**
 Here you can choose to pass as many LAMMPS commands as desired to your simulation.
 
 Note that there are some compulsory commands (see pair_styles or time-integration sections below)
 which must always be passed to ensure that TriLMP runs correctly.
 
 The LAMMPS commands can be passed in two different ways:
- *Pre-equilibration commands*: Pass them directly to the TriLmp object before the simulation has started running (this is useful to equilibrate your system beforehand)
- *Post-equilibration commands*: Pass them by appending them to the postequilibration_commands list below (they will be run once the pre-equilibration stage has ended, in order)
  
**C. Run**
Run the simulation by calling the TriLmp.run(args) method with its corresponding arguments.

## 2. Equilibrate the membrane from generated configurations
