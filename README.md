# Light PBH Spikes

<img width="1366" height="320" alt="Light PBH spikes-2" src="https://github.com/user-attachments/assets/67cc5468-cba6-4619-ab92-de54e353fc8c" />

Based on arXiv:26XX.XXXXX.

**Authors:** Agnese Tolino, Bradley Kavanagh, Francesca Scarcella, Valentina De Romeri, and Daniele Gaggero.

## Overview

This repository accompanies arXiv:26XX.XXXXX and investigates whether primordial black holes (PBHs) in the mass range of approximately `10^17 – 10^28 g` can form long-lived density spikes around solar-mass PBHs (`M_\odot`).

For such spikes to survive, infalling PBHs must acquire sufficient angular momentum to avoid direct capture by the central black hole. 

Since PBHs may form with negligible peculiar velocities, additional mechanisms are needed to generate the required torque. This project studies those mechanisms through two complementary approaches:

1. **Numerical simulations (GADGET4) of mutual gravitational interactions**  
   Dedicated simulations in GADGET4 of PBH gravitational interactions in the presence of a fixed central potential, designed to capture small-scale encounters and their role in angular momentum growth around a heavier PBH.

2. **Analytical modeling of other contributions**  
   Notebook implementations of the formulas derived in the paper, describing additional sources of torque such as large-scale dark matter perturbations and shell capture. These notebooks also reproduce the main figures of the study.


## Repository Structure

- **`gadget_fixed_potential`**  
  A customized GADGET4 setup for simulating light-PBH interactions in the presence of a fixed central gravitational potential.

- **`paper_notebooks`**  
  Jupyter notebooks containing the analytical calculations presented in the paper, together with scripts to reproduce its main results and figures.
