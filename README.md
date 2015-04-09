# KernelFit

C++ classes for single and multidimensional non-parametric Gaussian kernel
regression. These objects allow for the fitting of smooth profiles through
noisy data. The functions that solve for the profile are 
*embarrassingly parallel* and use OpenMP to gain large speedups.

**dependencies:**
You should be able to include 
![example](Figures/KernelFit1D.png "Results of KernelFit1D")
