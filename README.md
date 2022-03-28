<!-- ![OptDesign](https://myoctocat.com/assets/images/base-octocat.svg)
![OptDesign](https://github.com/chang88ye/OptDesign/blob/main/ToC.png) -->

<p float="centre">
  <img src="https://myoctocat.com/assets/images/base-octocat.svg" width="20%" />
  <img src="https://github.com/chang88ye/OptDesign/blob/main/ToC.png" width="75%" /> 
</p>

# OptDesign

OptDesign (Optimal Network Design) is a software platform for identifying genetic modification strategies, including knockout and up/down-regulation, for metabolic engineering design. It is based on genome-scale metabolic models (GSMM), and considers a two-step
procedure to identify the optimal combination of manipulations. The first step essentially reduces the search space by finding the most probable regulation targets. The second step considers the cell and a metabolic engineer as two different agenets playing a metabolic game: 
the metabolic engineer attempts to maximally violate the host cell's intention in avoiding overproducing biochemicals for homeostasis. 

# Getting Started
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

Prerequisites
  - MATLAB2018b and later
  - cobratoolbox (latest version: 3.0)
  - Gurobi8.0 above (or other software packages containing LP and ILP solvers)
Note: you need to use the new solveCobraMILP.m to replace the one in cobratoolbox if your solver is gurobi.

# Installation
1, Download the whole OptDesgin package directly from this repository or clone using

$ git clone --depth=1 https://github.com/Chang88ye/OptDesign.git OptDesign

2, Launch MATLAB and change the working directory to OptDesign as the current folder. Then run optdesign_setup
>> optdesign_setup

Note: optdesign_setup will help you add all the OptDesign files to matlab directory. It will also test if you have all prerequisites available for OptDesign.

3, Run a toy example that maximises succinate production from the E. coli core model
>> run_OptDesign

For the core model, this produces two flux distribution maps, one for the wild type and one for the mutant strain. You should be able to use the maps to inspect flux changes.

The solution structure of OptDesign is:

    - KOSet: set of knockcout reactions.
    
    - UPset: set of upregulation reactions whose absolute flux increases. Relative flux change direction is indicated by '+' (flux increase) or '-' (flux decrease) right after a reaction name in the set.
    
    - DOWNset: set of downregulation reactions whose absolute flux decreases. Relative flux change direction is indicated by '+' (flux increase) or '-' (flux decrease) right after a reaction name in the set.

4, If you want to apply OptDesign to other genome-scale models for different products, just change target product in run_OptDesign. 

5, This is an optional step, however, you need to run creat_envelopes.m in the envelopes folder if you want to create production envelopes for the design strategies identified by OptDesign.

# License
Released under the MIT license. All included network models provided under their respective licenses.

# Cite Us
Jiang S, Otero-Muras I, Banga JR, Kaiser M, and Krasnogor N (2021) OptDesign: Identifying optimal  design strategies in strain engineering for biochemical production. doi: https://doi.org/10.1101/2021.12.10.472123
