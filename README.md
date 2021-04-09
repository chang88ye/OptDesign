# OptDesign

OptDesign (Optimal Network Design) is a software platform for identifying genetic modification strategies, including knockout and up/down-regulation, for metabolic engineering design. It is based on genome-scale metabolic models (GSMM), and considers a two-step
procedure to identify the optimal combination of manipulations. The first step essentially reduces the search space by finding the most probable regulation targets. The second step considers the cell and a metabolic engineer as two different agenets playing a metabolic game: 
the metabolic engineer attempts to maximally violate the host cell's intention in avoiding overproducing biochemicals for homeostasis. 

# Getting Started
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

Prerequisites
MATLAB2016b and later
cobratoolbox (latest version: 3.0)
Gurobi8.0 above (or other software packages containing LP and ILP solvers)

# Installation
1, Download the whole OptDesgin package directly from this repository or clone using

$ git clone --depth=1 https://github.com/Chang88ye/OptDesign.git OptDesign

2, Launch MATLAB and change the working directory to NIHBA as the current folder. Then run optdesign_setup
>> optdesign_setup

Note: optdesign_setup will help you add all the NIHBA files to matlab directory. It will also test if you have all prerequisites available for OptDesign.

3, Run a toy example that maximises succinate production from the E. coli core model
>> run_OptDesign

For the core model, this produces two flux distribution maps, one for the wild type and one for the mutant strain. You should be able to use the maps to inspect flux changes.

4, If you want to apply OptDesign to other genome-scale models for different products, just change target product in run_OptDesign.

# License
Released under the MIT license. All included network models provided under their respective licenses.

# Cite Us
Jiang S, Otero-Muras I, Banga JR, Kaiser M, and Krasnogor N (2021) OptDesign: Identifying optimal  design strategies in strain engineering for biochemical production. 
