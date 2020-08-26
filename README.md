# Time_Delay_DMD
Code to supplement the paper:

"Principal Component Trajectories (PCT): Nonlinear dynamics as a superposition of time-delayed periodic orbits"

Daniel Dylewsky, Eurika Kaiser, Steven L. Brunton, and J. Nathan Kutz

https://arxiv.org/abs/2005.14321

## About:
This repository contains MATLAB code for a few simple examples of the methods outlined in this paper.

* Example 1 illustrates how delay embedding of time series data can be used to construct a basis of Principal Component Trajectories (PCTs) for sparse time-frequency representation of nonlinear dynamics

* Example 2 shows how nonlinear dynamics projected onto these PCTs can be linearly modeled with Dynamic Mode Decomposition (DMD)

* Example 3 demonstrates how the time_delay_dmd() function can be used to iterate through different combinations of delay embedding parameters

* Example 4 illustrates the procedure for unsupervised discovery of a forcing signal which can be used to construct a linear control model from the DMD results

More detailed documentation is included in the function files.

All code is standalone, but it is highly recommended that the user download the Optimized DMD package for MATLAB: https://github.com/duqbo/optdmd. To use optDMD with this code, follow its setup instructions and then execute the time_delay_dmd function with the flag dmd_type=='opt' instead of dmd_type=='exact'

