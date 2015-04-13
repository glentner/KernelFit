// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// KernelFit/Source/Test/TestKernelFit1D.cc

// Test of the KernelFit1D object. This test involves creating a `noisy`
// normalized sinc function and fitting a smooth profile through it.
// The output is are data files, `raw-1D.dat` and `fit-1D.dat`.

#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <fstream>
#include <omp.h>

#include "../KernelFit.hh"

#define N     1000
#define pi    3.141592653589793
#define range 3.0
#define scale 0.025

int main(){
	
	// create a uniform random line-space
	std::random_device device;
	std::mt19937_64 random_number( device( ) );
	std::vector<double> x(N, 0.0);
	for (std::size_t i = 0; i < x.size(); i++)
		x[i] = -(range * pi / 2.) + (range * pi * random_number() / 
			random_number.max());
	
	// solve the normalized sinc(x)
	std::vector<double> y(N, 0.0);
	for (std::size_t i = 0; i < y.size(); i++)
		y[i] = sin(pi * x[i]) / (pi * x[i]);
	
	// random `red` noise proportional to x
	for (std::size_t i = 0; i < y.size(); i++)
		y[i] += ( (2.0 * random_number() / random_number.max() - 1.0) 
			* scale * x[i] );
	
	// white noise
	for (std::size_t i = 0; i < y.size(); i++)
		y[i] += ( (2.0 * random_number() / random_number.max() - 1.0)
			* 2.5 * scale );
	
	// uniform line-space over range in x
	std::vector<double> xx(N, 0.0);
	xx[0] = -(range * pi / 2.0);
	for (std::size_t i = 1; i < xx.size(); i++)
		xx[i] = xx[i-1] + range * pi / (double(N) - 1);
	
	// set parallelism
	omp_set_num_threads(2);
	
	// solve KernelFit
	KernelFit1D<double> kernel(x, y, 3*pi * pi * range / double(N));
	std::vector<double> f = kernel.Solve(xx);
	
    // solve for standard deviation
    std::vector<double> s = kernel.StdDev(xx);
    
	// output raw data
	std::ofstream rawfile("Test/raw-1D.dat");
	if (rawfile) {
		
		for (std::size_t i = 0; i < x.size(); i++)
			rawfile << x[i] << " " << y[i] << std::endl;
		
	} else {
		
		std::cerr << "Failed to open output file, raw-1D.dat!\n";
		return 1;
	}
	
	// output profile fit
	std::ofstream fitfile("Test/fit-1D.dat");
	if (fitfile) {
		
		for (std::size_t i = 0; i < xx.size(); i++)
			fitfile << xx[i] << " " << f[i] << std::endl;
		
	} else {
		
		std::cerr << "Failed to open output file, fit-1D.dat!\n";
		return 1;
	}
    
    // output profile fit
    std::ofstream stdev("Test/stdev-1D.dat");
    if (stdev) {
        
        for (std::size_t i = 0; i < xx.size(); i++)
            stdev << xx[i] << " " << s[i] << std::endl;
        
    } else {
        
        std::cerr << "Failed to open output file, stdev-1D.dat!\n";
        return 1;
    }
    
	return 0;
}
