# Cooke Membrane Simulation & Fluctuation analysis example

Setup:

- [install micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)
- create environment with: `micromamba env create --name mr_cooke -f environment.yml`
- activate environment with `micromamba activate mr_cooke`

Check the simulation by running a small membrane patch for a short time with `lmp -var seed 1 -var halfL 12 -var duration 1000. -in in.lmp`; here `halfL` is half the simulation box size and `duration` is the simulation duration in time units $\tau$.
In a laptop computer this should take about 5mn.
This will allow a quick check of the software but will not provide good enough measurements.
To reproduce the results in the review paper, use `lmp -var seed 1 -var halfL 30 -var duration 60000. -in in.lmp`; this should complete in around 8h.
For these longer simulations, use a cluster, and for large box sizes parallelize with `mpi` (check your cluster documentation).

Optionally, inspect trajectory with [OVITO Basic](https://www.ovito.org/#download) by opening first the topology file `topo.data` and then using a `Load Trajectory Modifier` to load the `traj.nc` file.
This can be done before the simulation is complete, and also from a remote computer connected by SSH to the cluster.

Run analysis with `python analyze.py .`; expect a plot of box size $L$ versus time, and a fitted spectrum to saved as PDFs in the same folder.
Only equilibrated modes with at least 20 uncorrelated samples are taken in consideration and plotted.

The output files for the short simulation and its analysis were included in this repository, for comparison.
