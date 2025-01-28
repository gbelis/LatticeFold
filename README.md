# Protein Folding Simulation on a 2D Lattice

Lattice protein models are highly simplified representations of protein chains in a constrained lattice space, such as a 2D square or 3D cubic lattice. Each amino acid in the chain is modeled as a single point, and its placement is restricted to the vertices of the lattice. These models are crucial for studying the fundamental principles of protein folding without the complexity of continuous molecular dynamics simulations.

The simplest and most commonly used lattice protein model is the hydrophobic-polar (HP) model. This model stems from the observation that hydrophobic interactions among amino acid residues are one of the primary driving forces for proteins folding into their native states. The HP model simplifies protein sequences by categorizing amino acids into two types: **hydrophobic (H)** and **polar (P)**. 

This project simulates the folding process of proteins on a 2D lattice using a **Monte Carlo algorithm** combined with **simulated annealing** to explore low-energy configurations. The simulation incorporates move generation, energy calculation, and acceptance criteria to evaluate folding pathways and identify near-native conformations.

---

### Protein

Proteins are linear chains of amino acids that fold into specific three-dimensional shapes to perform essential biological functions, such as catalysis, structural support, and signaling. The sequence of amino acids, referred to as the **primary structure**, encodes all the information required for the protein to fold into its unique and functional structure. 

In this simulation, the HP model is employed to represent proteins. Each amino acid in the chain is classified as either **hydrophobic (H)** or **polar (P)**, depending on its interaction with water:
- **Hydrophobic residues (H)** tend to cluster together to minimize exposure to the surrounding polar environment, mimicking the hydrophobic effect observed in real proteins.
- **Polar residues (P)**, on the other hand, interact freely with the environment and do not contribute to stabilizing the folded structure.

The simulation uses a **2D lattice model**, where each amino acid occupies a vertex of a two-dimensional square lattice. Consecutive residues in the chain are constrained to be adjacent on the lattice. This simplification captures the key aspect of protein folding: **the trade-off between minimizing hydrophobic exposure and maintaining chain connectivity**.

The sequence of H and P residues and their spatial arrangement on the lattice determine the stability and energy of a given conformation. **Lower-energy configurations** are more stable and closer to the protein's native structure. Therefore, finding such configurations is the primary goal of this simulation.

---

### Move Generation

The simulation explores conformational space by generating moves within the **local neighborhood** of the current structure. The local neighborhood consists of all structures that can be reached by a single move from the current configuration. These moves ensure the chain remains connected and occupies valid lattice positions.

#### Types of Moves:
1. **End Move**:  
   This involves pivoting the residue at either end of the chain (residue 1 or residue \( n \)) relative to its connected neighbor. The pivot must place the residue at a free, adjacent position on the lattice.

2. **Corner Move**:  
   A corner move can be performed on residues (except the ends) when the two connected neighbors of a residue are mutually adjacent to an unoccupied lattice position. The residue can "snap" into this position, mimicking a local rearrangement.

3. **Crankshaft Move**:  
   This move is applicable when a residue is part of a U-shaped bend in the chain. The crankshaft move involves a 180° rotation of the U-shaped structure, provided the new positions are unoccupied. This move allows the chain to explore compact configurations efficiently.

4. **Chain Rotation**:  
   The chain can be rotated at any internal residue, where the rotation pivots the structure by \( \pm\frac{\pi}{2} \) or \( \pm\pi \). The rotation must maintain the chain's connectivity and prevent steric clashes.

#### Move Acceptance Criteria:
- **Steric Clash**: No two residues can occupy the same lattice position.
- **Chain Connectivity**: All residues must remain adjacent to their immediate neighbors in the chain.

These moves, combined with Monte Carlo sampling, allow the simulation to explore the conformational space effectively.

---
### Energy

The energy of a given conformation is determined by the **number of hydrophobic-hydrophobic (H-H) interactions** within the structure. An **H-H interaction** occurs when two hydrophobic residues that are not neighbors in the chain are adjacent on the lattice. Each H-H interaction contributes an energy of **-1**, promoting compact configurations that cluster hydrophobic residues together.

#### Total Energy:
$E = -\sum \text{(H-H interactions)}$


The energy calculation captures the essential physics of protein folding: hydrophobic residues seek to minimize their exposure to the polar environment by clustering together, forming stable, low-energy states. Simulated annealing is used to guide the simulation toward these states by balancing exploration (high temperatures) and exploitation (low temperatures).

---

### Monte Carlo Move Selection

The Monte Carlo simulation is a stochastic approach to exploring the conformational space of the protein. At each iteration, a move is proposed from the set of possible moves (end moves, corner moves, crankshaft moves, or chain rotations). The decision to accept or reject the proposed move is determined using the **Metropolis criterion**, which balances exploration and optimization.

#### Metropolis Criterion:
A move is accepted with the probability:

![Equation](https://latex.codecogs.com/svg.image?\color{white}P_{\text{accept}}=\begin{cases}1&\text{if}\Delta&space;E\leq&space;0\\\exp\left(-\frac{\Delta&space;E}{T}\right)&\text{if}\Delta&space;E>0\end{cases})

Where:
- $\Delta E$ is the change in energy caused by the proposed move $E_{\text{new}} - E_{\text{current}}$.
- $T$ is the system's temperature at the current iteration.

#### Steps in Monte Carlo Move Selection:
1. **Propose a Move**:  
   A random move is chosen from the available set (end, corner, crankshaft, or rotation).
   
2. **Check Validity**:  
   The proposed move is validated by ensuring:
   - **Steric clash avoidance**: No two residues occupy the same lattice position.
   - **Chain connectivity**: All residues remain adjacent to their neighbors.

3. **Calculate Energy Change**:  
   The energy of the new conformation is computed, and the change in energy (\( \Delta E \)) is determined.

4. **Acceptance Decision**:  
   - If $\Delta E \leq 0$, the move is always accepted, as it leads to a lower-energy state.
   - If $\Delta E > 0$, the move is accepted with a probability proportional to $\exp(-\Delta E / T)$, allowing uphill moves to escape local minima.

5. **Update State**:  
   If the move is accepted, the conformation is updated to the new state. Otherwise, the system retains its current state.

#### Importance of Temperature $T$:
The temperature parameter $T$ controls the balance between exploration and exploitation:
- At **high temperatures**, the system is more likely to accept uphill moves, enabling broad exploration of the conformational space.
- At **low temperatures**, the system predominantly accepts downhill moves, focusing on refining low-energy states.

#### Simulated Annealing:
The simulation uses a cooling schedule to gradually lower the temperature $T$, reducing the acceptance of high-energy states over time. This process mimics the natural cooling process during protein folding, allowing the simulation to converge to near-native conformations.

This Monte Carlo move selection process, combined with the energy model, enables the simulation to explore a vast number of potential folding pathways efficiently and identify low-energy configurations representative of the protein's folded state.

---

### Figures

#### Figure 1: Example of the Folding Process Using Single Residue Updates
![protein_folding_small_updates](https://github.com/gbelis/LatticeFold/blob/main/visualisations/protein_folding_small_updates.gif)  
This figure demonstrates the gradual folding process, where individual residues are repositioned to achieve lower-energy configurations.

#### Figure 2: Example of the Folding Process Using Rotations
![protein_folding_rot](https://github.com/gbelis/LatticeFold/blob/main/visualisations/protein_folding_rot.gif)  
This figure illustrates larger-scale conformational changes facilitated by chain rotations, enabling the chain to explore distant regions of the conformational space.

*Note: For rendering efficiency, only 1 in every 5 states is shown in the visualizations.*

---

### Methodological Justifications

1. **HP Model**:  
   While highly simplified, the HP model effectively captures the essence of protein folding, where hydrophobic interactions dominate. It provides a computationally efficient framework for exploring folding dynamics.

2. **2D Lattice**:  
   Restricting the simulation to a 2D lattice reduces computational complexity, allowing for the study of folding mechanisms without the overhead of 3D modeling.

3. **Monte Carlo Sampling**:  
   Monte Carlo methods are ideal for sampling vast conformational spaces, ensuring diverse configurations are explored.

4. **Simulated Annealing**:  
   This technique mimics the natural folding process by gradually lowering the "temperature" of the system, focusing on lower-energy states while avoiding local minima.

---

### Bibliography
1. [A replica exchange Monte Carlo algorithm for protein folding in the HP model](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-342)  
   This paper discusses advanced Monte Carlo algorithms and their application to HP models, offering insights into efficient sampling techniques.

2. [Protein Folding Prediction in a Cubic Lattice in Hydrophobic-Polar Model](https://www.liebertpub.com/doi/full/10.1089/cmb.2016.0181)  
   This article provides a comprehensive review of protein folding predictions using lattice-based HP models, emphasizing algorithmic and computational considerations.

--- 
