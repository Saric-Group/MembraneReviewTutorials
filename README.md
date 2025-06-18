# A guide to modeling mesoscale membrane deformations with coarse-grained computer simulations

This repository contains ready-to-run tutorials on the simulation tests presented in the membrane review by M. Muñoz-Basagoiti, F. Frey,  B. Meadowcroft, M. Amaral, A. Prada, and A. Saric[^RevCit] ([2502.09798](https://arxiv.org/abs/2502.09798); under review).

## Table of contents
1. [Tutorials](#tutorials)
    1. [Three-beads-per-lipid membrane: Cooke model](#cooke)
    2. [One-particle thick fluid membrane: YLZ model](#ylz)
    3. [Dynamically triangulated membrane model I: Monte Carlo simulations with C code](#mcsims)
    4. [Dynamically triangulated membrane model II: Parallelized Hybrid Monte Carlo simulations with TriLMP](#trilmp) 
2. [Reporting bugs and feedback](#bugs)

## Tutorials <a name="tutorials"></a>

### 1. Three-beads-per-lipid membrane: Cooke model <a name="cooke"></a>
[See respective folder](CookeSimulations).

### 2. One-particle thick fluid membrane: YLZ model <a name="ylz"></a>
[See respective folder](YLZSimulations).

### 3. Dynamically triangulated membrane model I: Monte Carlo simulations with C code  <a name="mcsims"></a>
[See respective folder](MCSimulations).

### 4. Dynamically triangulated membrane model II: Parallelized HMC simulations with TriLMP  <a name="trilmp"></a>
Simulate dynamically triangulated network in parallel by using TriLMP. TriLMP is a modified version of the [TriMEM](https://github.com/bio-phys/trimem) package[^Siggel2022] developed within the Šarić Group at ISTA. To be able to run the tutorials of this section, you will have to install TriLMP. For that, please go to the [TriLMP Github repository](https://github.com/Saric-Group/trimem_sbeady) and follow the installation guide for your operating system. Once you have verified that your installation is working fine, go to [HMCSimulations](https://github.com/Saric-Group/MembraneReviewTutorials/tree/main/HMCSimulations) to run the tutorials of this section.

## Reporting bugs and feedback  <a name="bugs"></a>

Please report any bugs. Feedback on how to improve these tutorials is welcome.

[^RevCit]: Munoz-Basagoiti, Frey, Meadowcroft, Amaral, Prada and Saric, A guide to modeling mesoscale membrane deformations with coarse-grained computer simulations (submitted, 2025)
[^Siggel2022]: Siggel, M. et al. (2022) TriMem: A Parallelized Hybrid Monte
  Carlo Software for Efficient Simulations of Lipid Membranes.
  J. Chem. Phys. (in press) (2022); https://doi.org/10.1063/5.0101118
