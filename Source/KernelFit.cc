// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// KernelFit/Source/KernelFit.cc

// This source file contains the definitions for the KernelFit objects.

#include <vector>
#include <omp.h>

#include "KernelFit.hh"

template<class T>
KernelFit1D<T>::KernelFit1D(const std::vector<T> &x, const std::vector<T> &y,
	const double bandwidth){
		
	// Constructor for the KernelFit1D object. Save the `x`, `y` data
	// and set an initial value for the `bandwidth`.
		
	_x = x;
	_y = y;
	
	_bandwidth = bandwidth;
	
}

template<class T>
T KernelFit1D<T>::Dispersion(){
	
	// Solve for the mean distance between the data points.
	// It is assumed that the `x` vector is in ascending order!
	
	T dispersion = 0;
	
	for (std::size_t i = 1; i < _x.size(); i++)
		dispersion += _x[i] - _x[i-1];
	
	return dispersion / _x.size();
}

template<class T>
void KernelFit1D<T>::SetBandwidth(const double multiple){
	
	// set the bandwidth to be a multible of the Dispersion()
	_bandwidth = multiple * Dispersion();
}

template<class T>
std::vector<T> KernelFit1D<T>::Solve(const std::vector<T> &x){
	//
	// solve for the smooth profile through the data at all `x`
	//
	
	std::vector<T> f( x.size(), 0.0);
	
	// omp_set_num_threads() should be called prior to here!
	#pragma omp parallel for shared(f)
	for (std::size_t i = 0; i <  x.size(); i++){
	
		T sum = 0.0;
		
		for (std::size_t j = 0; j < _x.size(); j++){
			
			T k   = Kernel(_x[j] - x[i]);
			f[i] += k * _y[j];
			sum  += k;
		}
		
		f[i] /= sum;
	}
	
	return f;
}

template<class T>
KernelFit2D<T>::KernelFit2D(const std::vector<T> &x, const std::vector<T> &y,
	const std::vector<T> &z, const double bandwidth){
		
	// Constructor for the KernelFit1D object. Save the `x`, `y` data
	// and set an initial value for the `bandwidth`.
		
	_x = x;
	_y = y;
	_z = z;
	
	_bandwidth = bandwidth;
	
}

template<class T>
T KernelFit2D<T>::Dispersion(){
	
	// Solve for the mean distance between the data points.
	
	T dispersion = 0;
	
	for (std::size_t i = 1; i < _x.size(); i++)
		dispersion += _x[i] - _x[i-1];
	
	return dispersion / _x.size();
}

template<class T>
void KernelFit2D<T>::SetBandwidth(const double multiple){
	
	// set the bandwidth to be a multible of the Dispersion()
	_bandwidth = multiple * Dispersion();
}

template<class T>
std::vector<T> KernelFit2D<T>::Solve(const std::vector<T> &x,
	const std::vector<T> &y){
	//
	// solve for the smooth profile through the data at all (x, y)
	//
	
	std::vector<T> f( x.size(), 0.0);
	
	// omp_set_num_threads() should be called prior to here!
	#pragma omp parallel for shared(f)
	for (std::size_t i = 0; i <  x.size(); i++){
	
		T sum = 0.0;
		
		for (std::size_t j = 0; j < _x.size(); j++){
			
			T k   = Kernel(_x[j] - x[i]);
			f[i] += k * _y[j];
			sum  += k;
		}
		
		f[i] /= sum;
	}
	
	return f;
}

template class KernelFit1D<double>;