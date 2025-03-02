# A guide to modeling mesoscale membrane deformations with coarse-grained computer simulations

This repository contains ready-to-run tutorials on the simulation tests presented in the membrane review by B. Meadowcroft, F. Frey, M. Muñoz-Basagoiti, A. Prada, M. Amaral and A. Saric[^RevCit].

## Table of contents
1. [Tutorials](#tutorials)
    1. [Three-beads-per-lipid membrane: Cooke model](#cooke)
    2. [One-particle thick fluid membrane: YLZ model](#ylz)
    3. [Dynamically triangulated membrane model: TriLMP](#trilmp)
2. [Reporting bugs and feedback](#bugs)

## Tutorials <a name="tutorials"></a>

### 1. Three-beads-per-lipid membrane: Cooke model <a name="cooke"></a>
[See respective folder](CookeSimulations).

### 2. One-particle thick fluid membrane: YLZ model <a name="ylz"></a>
[See respective folder](YLZSimulations).

### 3. Dynamically triangulated membrane model: TriLMP  <a name="trilmp"></a>
We choose TriLMP to simulate a fluid membrane using a Dynamically Triangulated Network (DTN). TriLMP is a modified version of the [TriMEM](https://github.com/bio-phys/trimem) package[^Siggel2022] developed within the Saric Group at ISTA. To be able to run the tutorials of this section, you will have to install TriLMP. For that, please go to the [TriLMP Github repository](https://github.com/Saric-Group/trimem_sbeady) and follow the installation guide for your operating system. Once you have verified that your installation is working fine, go to [DTNSimulations](https://github.com/Saric-Group/MembraneReviewTutorials/tree/main/DNTSimulations) to run the tutorials of this section.

## Reporting bugs and feedback  <a name="bugs"></a>

Please report any bugs. Feedback on how to improve this tutorials is welcome.

[^RevCit]: Munoz-Basagoiti, Frey, Meadowcroft, Amaral, Prada and Saric, A guide to modeling mesoscale membrane deformations with coarse-grained computer simulations (submitted, 2025)
[^Siggel2022]: Siggel, M. et al. (2022) TriMem: A Parallelized Hybrid Monte
  Carlo Software for Efficient Simulations of Lipid Membranes.
  J. Chem. Phys. (in press) (2022); https://doi.org/10.1063/5.0101118
