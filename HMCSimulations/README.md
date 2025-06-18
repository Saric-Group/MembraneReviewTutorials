# Tether extrusion with a Dynamically Triangulated Network II: Parallelized Hybrid Monte Carlo Simulations with TriLMP

It is possible to simulate the configurations of a dynamically triangulated mesh in parallel, by using software such as TriMem (see [Siggel *et al.*, 2023](https://pubs.aip.org/aip/jcp/article/157/17/174801/2842373/TriMem-A-parallelized-hybrid-Monte-Carlo-software)) or [PyMembrane](https://github.com/fdmatoz/PyMembrane) (see [Matoz-Fernandez *et al.*, 2023](https://arxiv.org/abs/2308.12754)). In this tutorial, we show how to measure the force required to extrude a tube from a membrane patch using the [TriLMP software](https://github.com/Saric-Group/trimem_sbeady)

TriLMP is a modified version of the [TriMEM](https://github.com/bio-phys/trimem) package[^Siggel2022] developed within the Šarić Group at ISTA. To be able to run the tutorials of this section, you will have to install TriLMP. For that, please go to the [TriLMP Github repository](https://github.com/Saric-Group/trimem_sbeady) and follow the installation guide. Once you have verified that your installation is working fine, you are ready to run the tutorials of this section.

To efficiently extrude a tube from a membrane patch, the simulation procedure is separated in two steps: (1) the generation of an initial configuration and (2) the relaxation of the tube to its equilibrium position. Note that this tutorial is based on a membrane patch with fixed boundaries, which is a different set-up from that detailed in the [MCSimulations tutorial](../MCSimulations) and included the review, where a bulk membrane patch can draw area from its boundary. As the boundary conditions will affect the force-elongation profile measured; we recommend the simulation set-up and boundary conditions described in the review to extrude the tube. This tutorial can nonetheless be useful for the the reader to get acquainted with TriLMP as a simulation method.

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
Once the first stage of our analysis has concluded and we have generated enough configurations, we use those configurations as new starting points. Our goal now is to let the system relax in order to measure the equilibrium force that the membrane exerts on the pulling bead, for different bead positions $z = z_i$. For each initial configuration, the membrane will reach a certain maximum elongation and exert a certain force on the pulling bead. It is these observables that we use to generate the force-elongation profile in OurReviewEtAl., (2024). 

From the computational point of view, the only challenge in this section is to load the (pickled) checkpoints we have created. This can be easily done through the `read_checkpoint(args)` function in TriLMP. While the basic structure of a TriLMP program (see above) is not explicit anymore, it still holds, and the three stages can be easily recognised. Additionally, the LAMMPS commands that we have used for the first part are also recycled here. The difference now is that the pulling bead, which before was mobile, remains static.

### 2.1. Running the simulation
Given a checkpoint file, to run its equilibration go to the [part2_configequilibration](https://github.com/Saric-Group/MembraneReviewTutorials/tree/main/DNTSimulations/part2_configequilibration) directory. You can run the simulation directly by using:

```python launch_configequil.py```

Note that to obtain the full force-elongation profile you will need to equilibrate the system at many different configurations (many different elongations).

