# backmapping
Backmapping from Multiresolution Coarse-Grained Models to Atomic  Structures of Large Biomolecules by Restrained Molecular Dynamics  Simulations Using Bayesian Inference

Scripts and packages in this folder
1. "gmx455-blogr"
    - We added log harmonic energy term to GROMACS-4.5.5. Please install this modified package before running backmapping simulations.
    - Installation of the modified package is the same as that of GROMACS-4.5.5

2. "cgmap.py"
    - it's used for mapping aa atoms to CG sites and generating inputs for backmapping simulations
    - the script is only tested for python 2.7 & 3.7. There might have some issues for other python versions. Please contact us if there would have some bugs.

3. "blogr.sh"
    - bash script to run backmapping simulations
    - modify python and gromacs path and backmapping simulation parameters accordingly

4. "mdtemp.mdp"
    - MD parameter file. Modify it accordingly.
    
Backmapping simulations
1. Setup your simulation system. As an initial AA structure is needed in our method. We need to setup the simulation. We recommend using Amber force fields, since currently cgmap.py is only compatible with Amber force field.

2. Edit the main script blogr.sh.
    - Enviorements: Python path and modified Gromacs path.
    - Input files: Initial AA structure, Initial <name>.GRO (with velocities), Target CG structure, CG type (CA, COM, and others, run python cgmap.py -h for detail) and topology files.

3. run blogr.sh

4. After simulations, if neccessary, use clustering method to extract representative conformations as the final solution. This can be done by g_cluster implemented in Gromacs or any other clustering algorithms.

Reference:
Backmapping from Multiresolution Coarse-Grained Models to Atomic Structures of Large Biomolecules by Restrained Molecular Dynamics Simulations Using Bayesian Inference  
Junhui Peng, Chuang Yuan, Rongsheng Ma, and Zhiyong Zhang  
Journal of Chemical Theory and Computation 2019 15 (5), 3344-3353  
DOI: 10.1021/acs.jctc.9b00062
