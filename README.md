# Detecting 3D Hand Pointing Direction from RGB-D Data in Wide-Ranging HRI Scenarios

This repository contains data and codes for the paper ..... The presented method is an adaptation of the method presented in __Efficient robust optimization for robust control
with constraints__ for a system of the form $` h_{i+1}= A_dh_i + B_{d_1}u_i + B_{d_2} d_{a_i} + w_i`$ where $`d_{a_i}`$ is some known function of time. The paper is available [here](https://link.springer.com/article/10.1007/s10107-007-0096-6). If you use this code or parts of it, please cite the paper.

## Abstract

## Installation
The optimization problems are implemented using [JuMP](https://github.com/JuliaOpt/JuMP.jl) modeling language for mathematical optimization.

 [Ipopt](https://ipoptjl.readthedocs.io/en/latest/ipopt.html) solver is used in the code, however, it may be changed with a QP solver when feasible.


