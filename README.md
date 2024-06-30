# ...

This repository contains the codes for the paper ..... The presented method is an adaptation of the computationally efficient robust optimization method presented in __Efficient robust optimization for robust control
with constraints__ (available [here](https://link.springer.com/article/10.1007/s10107-007-0096-6)). Note that although the method is applied to a Water Distribution Network (WDN), it can be applied to any system of the form $` h_{i+1}= Ah_i + B_{1}u_i + B_{2} d_{i} + w_i`$ where $`d_{i}`$ is a known function of time and $`w_i`$ is bounded the unknown disturbance.
## Abstract

## Files
-`dfStdStruct.jl`: The standard implementation of the presented robust controller <br /> 
-`dFEffStruct.jl`: The computationally efficient implementation of the presented robust controller <br />
-`dfEFFmain.jl`: An example usage of the presented methods in a receding horizon fashion <br />
-`parameters.jl`: Parameters for the WDN in the paper.



