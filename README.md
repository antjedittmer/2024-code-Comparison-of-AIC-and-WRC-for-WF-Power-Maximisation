# 2023-code-IFAC-Koopman Model Predictive Control for Wind Farm Yield Optimization with Combined Thrust and Yaw Control

## General

This folder contains supplementary documentation and code to reproduce the results and figures presented in the paper 


> A. Dittmer, B. Sharan and H. Werner, "Comparison and Combination of Axial Induction and Wake Redirection Control for Wind Farm Power Output Maximisation and Grid Power Reference Tracking"

accepted at the 23rd Wind & Solar Integration Workshop, 2024.

It may be used to recreate the simulation results and figures from the paper. To do so, run the script `mainGeneratePlotsWISO.m`.

Running the simulation for the first time takes roughly 20 minutes using a laptop with an 11th Gen Intel(R) Core(TM) i7-11850H @ 2.50GHz processor.

## Optimization Environment
FLORIS (FLOw Redirection and Induction in Steady State)[on GitHub](https://github.com/NREL/floris) is utilized for yaw optimization in steady-state. 

## Simulation Environment WFSim

The simulations WindFarmSimulator (WFSim) developed by Doekemeijer and Boersma [on GitHub](https://github.com/TUDelft-DataDrivenControl/WFSim) is utilized.
WFSim is used both for data generation in open loop and for closed loop testing.

## Simulation Environment FAST.Farm
The simulation FAST.Farm is part of OpenFAST and available here [on GitHub](https://github.com/OpenFAST/openfast). It is used for load analysis.

## Evaluation 

The code in this repository was tested in the following environment:

* *Windows* 10 Enterprise 21H2
* *Matlab* 2021b (Toolboxes used: Control System Toolbox, Optimization Toolbox, and Global Optimization Toolbox)


