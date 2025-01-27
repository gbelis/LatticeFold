# Protein Folding Simulation on a 2D Lattice

Lattice protein models are highly simplified models of protein chains on lattice conformational space. Each whole amino acid is modeled as a single point and is restricted to be placed on vertices of a lattice (square or cube). The simplest model for modeling protein is the hydrophobic-polar (HP) protein model. The HP protein model is a highly simplified model for studying protein folds. It stems from the observation that hydrophobic interactions between amino acid residues are the driving force for proteins folding into their native state.

This project simulates the folding process of proteins on a 2D lattice using a Monte Carlo model. The simulation employs methods like move generation, energy calculation, and simulated annealing to explore possible configurations of the protein.

### Figure 1: Example of the folding process using single residue updates.
![protein_folding_small_updates](https://github.com/gbelis/LatticeFold/blob/main/visualisations/protein_folding_small_updates.gif)


### Figure 2: Example of the folding process using rotations.
![protein_folding_rot](https://github.com/gbelis/LatticeFold/blob/main/visualisations/protein_folding_rot.gif)

* For speed of rendering reasons, the gifs were generated using only 1 of every 5 states.


## Bibliography
- [A replica exchange Monte Carlo algorithm for protein folding in the HP model](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-342)
- [Protein Folding Prediction in a Cubic Lattice in Hydrophobic-Polar Model](https://www.liebertpub.com/doi/full/10.1089/cmb.2016.0181)
