# A Robust Predictive Control Method  for Pump Scheduling in Water Distribution Networks

This repository contains the code accompanying the paper **" A Robust Predictive Control Method  for Pump Scheduling in Water Distribution Networks"**. The robust control method implemented here is based on the computationally efficient robust optimization technique presented in the paper ["Efficient robust optimization for robust control with constraints"](https://link.springer.com/article/10.1007/s10107-007-0096-6).

## Overview

The method in this repository is an extension of the robust optimization method developed in the aforementioned paper. While the core steps and approach remain the same, our method includes an additional term in the system model to better represent Water Distribution Networks (WDNs). Despite this addition, the fundamental procedure of the method is largely unchanged and retains its computational efficiency and robustness.

### System Model

The robust control method is designed for systems that can be described by the following state-space model:
```math
 h_{i+1} = Ah_i + B_{1}u_i + B_{2}d_i + w_i 
```
where:
-   $`h_i `$  is the state vector at time step  $`i`$.
-   $`A`$ is the state transition matrix.
-   $`u_i`$ is the control input vector at time step  $`i`$.
-   $`B_1`$ and  $`B_2`$  are input matrices.
-   $`d_i`$ is a known function of time.
-   $`w_i`$ is a bounded unknown disturbance.

### Applicability

While this implementation is specifically applied to a Water Distribution Network (WDN), the robust control method can be applied to any system that can be modeled in the above form.


## Files

- `dfStdStruct.jl`: The standard implementation of the presented robust controller.
- `dFEffStruct.jl`: The computationally efficient implementation of the presented robust controller.
- `dfEFFmain.jl`: An example usage of the presented methods in a receding horizon fashion.
- `parameters.jl`: Parameters for the WDN in the paper.

## How to cite this material?

https://doi.org/10.5281/zenodo.13908277

