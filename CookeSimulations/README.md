# Cooke Membrane Simulation & Fluctuation analysis example

Setup LAMMPS and python for analysis:

- install `micromamba`
- create environment with: `micromamba env create --name mr_cooke -f environment.yml`
- activate environment with `micromamba activate mr_cooke`

Edit `in.lmp` to adjust membrane size and duration; small values will allow a quick check of the software but will not provide good enough measurements.
Run the simulation `lmp -in in.lmp`. 

For long simulations, use a cluster and parallelize LAMMPS with `mpi` according to the cluster documentation.

Optionally, inspect trajectory with ovito by opening first the topology file `topo.data` and then using a `Load Trajectory Modifier` to load the `traj.nc` file.

Run analysis with:

 `python analyze.py .`
