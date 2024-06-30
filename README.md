# Detecting 3D Hand Pointing Direction from RGB-D Data in Wide-Ranging HRI Scenarios

This repository contains data and codes for the paper ..... The presented method is an adaptation of the computationally efficient robust optimization method presented in __Efficient robust optimization for robust control
with constraints__ (available [here](https://link.springer.com/article/10.1007/s10107-007-0096-6))for a system of the form $` h_{i+1}= Ah_i + B_{1}u_i + B_{2} d_{i} + w_i`$ where $`d_{i}`$ is a known function of time and $`w_i`$ is the unknown disturbance. If you use this code or parts of it, please cite the paper.

## Abstract

## Installation
The optimization problems are implemented using [JuMP](https://github.com/JuliaOpt/JuMP.jl) modeling language for mathematical optimization.

 [Ipopt](https://ipoptjl.readthedocs.io/en/latest/ipopt.html) solver is used in the code, however, it may be changed with a QP solver when feasible.


