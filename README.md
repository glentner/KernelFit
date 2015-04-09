# KernelFit

C++ classes for single and multidimensional non-parametric Gaussian kernel
regression. These objects allow for the fitting of smooth profiles through
noisy data. The functions that solve for the profile are 
*embarrassingly parallel* and use OpenMP to gain large speedups.

**Dependencies:**
STL and OpenMP

![example](Figures/KernelFit1D.png "Results of KernelFit1D")

**Figure** **1:** The above figure was plotting using Python and showcases the 
results of the TestKernelFit1D.cc program. A noisy sinc function was produced 
with both *red* and *white* noise. The blue dashed line demonstrates the smooth 
profile fit through the data. The red dashed line shows the analytical function.
