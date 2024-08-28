# A guide to modeling mesoscale membrane deformations with coarse-grained computer simulations

This repository contains ready-to-run tutorials on the simulation tests presented in the membrane review by B. Meadowcroft, F. Frey, M. Muñoz-Basagoiti, A. Prada and M. Amaral.

## Three-beads-per-lipid membrane: Cooke model
Please add info needed to understand and run your tutorial(s).

## One-particle thick fluid membrane: YLZ model
Please add info needed to understand and run your tutorial(s).

## Dynamically triangulated membrane model: TriLMP
We use TriLMP to simulate a fluid membrane using a Dynamically Triangulated Network (DTN). TriLMP is a modified version of the Trimem python package(https://github.com/bio-phys/trimem). It couples Trimem to LAMMPS. It allows the direct use for MD simulations in connection with LAMMPS via the python interface of the latter. Hereby the calculation of the surface repulsion is dealt with by LAMMPS instead of Trimem. For details on how to use the package, please refer to the HowTo_TriLMP.md file. 
