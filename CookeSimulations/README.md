# Cooke Membrane Simulation & Fluctuation analysis example

<!-- This tutorial represents an opinionated view on how to do science with simulations.
Divert from its choices at your own risk. -->

## Pre-requisites

<!-- ###### Run Locally -->
We recommend first testing this software on a local computer, as opposed to running the code remotely on a cluster.

<!-- ##### Operating system -->
This software was tested on OSX and Linux.
If you are on a Windows machine please look into linux virtual machines or WSL (Windows Subsystem for Linux).

<!-- ##### Command line -->
This tutorial assumes familiarity with the command line and the shell `bash`.

## Software dependencies
For running simulations we use the simulator LAMMPS and for analyzing them we use Python and specific Python packages.
We install these in an isolated container using the `micromamba` package manager.
In a shell session, the container can then be activated, exposing the software within to the shell.

In short:
- [install micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)
- create environment with: `micromamba env create --name mr_cooke -f environment.yml`
- activate environment with `micromamba activate mr_cooke`

For visualization, install the basic GUI version of [Ovito](https://www.ovito.org/#download).

### Running a simulation

Check the simulation by running a small membrane patch for a short time with `lmp -var seed 1 -var halfL 12 -var duration 1000. -in in.lmp`; here `halfL` is half the simulation box size and `duration` is the simulation duration in time units $\tau$.
In a laptop computer this should take about 5mn.
This will allow a quick check of the software but will not provide good enough measurements.
To reproduce the results in the review paper, use `lmp -var seed 1 -var halfL 30 -var duration 60000. -in in.lmp`; this should complete in around 8h.
For these longer simulations, use a cluster, and for large box sizes parallelize with `mpi`.

To use `mpi` in a cluster, you must load the cluster's `mpi` library, usually via the `module` command, and also adjust the `environment.yml` file to include the [pseudo-package](https://conda-forge.org/docs/user/tipsandtricks/#using-external-message-passing-interface-mpi-libraries) that tells the package manager that `mpi` will be provided by the external environment.
The `mpi` version and implementation (`openmpi` or `mpich`) for both cluster and pseudo-package must be compatible.

Optionally, inspect trajectory with [OVITO Basic](https://www.ovito.org/#download) by opening first the topology file `topo.data` and then using a `Load Trajectory Modifier` to load the `traj.nc` file.
This can be done before the simulation is complete, and also from a remote computer connected by SSH to the cluster.

Run analysis with `python analyze.py .`; expect a plot of box size $L$ versus time, and a fitted spectrum to saved as PDFs in the same folder.
Only equilibrated modes with at least 20 uncorrelated samples are taken in consideration and plotted.

The output files for the short simulation and its analysis were included in this repository, for comparison.