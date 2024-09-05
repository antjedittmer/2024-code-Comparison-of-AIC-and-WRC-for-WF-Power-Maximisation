# 2024-code-Comparison-and-Combination-of-Axial-Induction-and-Wake-Redirection-Control-for-Wind-Farm-Power-Output-Maximisation-and-Grid-Power-Reference-Tracking

## General

This folder contains supplementary documentation and code to reproduce the results and figures presented in the paper 

> A. Dittmer, B. Sharan and H. Werner, "Comparison and Combination of Axial Induction and Wake Redirection Control for Wind Farm Power Output Maximisation and Grid Power Reference Tracking"

accepted at the 23rd Wind & Solar Integration Workshop, 2024.

It may be used to recreate the simulation results and figures from the paper. To do so, run the script `mainGeneratePlotsWISO.m`.

Running the simulation for the first time takes roughly 20 minutes using a laptop with an 11th Gen Intel(R) Core(TM) i7-11850H @ 2.50GHz processor. 
This is due to data in WFSim beiing generated to design the Koopman matrix used in the Koopman MPC. Once the Koopman matrix is available, the figures are generated in approximately 20 seconds.
The output of the FLORIS optimization and the FAST.FArm simulation is provided as *.mat files.  

This paper uses the same slightly modified WFSim environment and MPC design which we describe in  

[Koopman Model Predictive Control for Wind Farm Yield Optimization with Combined Thrust and Yaw Control](https://www.sciencedirect.com/science/article/pii/S2405896323014209)

## Simulation Environment WFSim

The simulations WindFarmSimulator (WFSim) developed by Doekemeijer and Boersma [on GitHub](https://github.com/TUDelft-DataDrivenControl/WFSim) is utilized.
WFSim is used both for data generation in open loop to design the Koopman matrix and for closed loop testing.

## Optimization Environment
FLORIS (FLOw Redirection and Induction in Steady State) from NREL (National Renewable Energy Laboratory) is available [on GitHub](https://github.com/NREL/floris) is utilized for yaw optimization in steady-state. 


## Simulation Environment FAST.Farm
The simulation FAST.Farm is part of OpenFAST from NREL and available here [on GitHub](https://github.com/OpenFAST/openfast). It is used for load analysis.

## Evaluation 

The code in this repository was tested in the following environment:

* *Windows* 10 Enterprise 21H2
* *Matlab* 2021b (Toolboxes used: Control System Toolbox, Optimization Toolbox, and Global Optimization Toolbox)


