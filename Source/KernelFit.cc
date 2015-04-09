// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// KernelFit/Source/KernelFit.cc

// This source file contains the definitions for the KernelFit objects.

#include <cmath>
#include <vector>
#include <algorithm>
#include <omp.h>

#include "KernelFit.hh"

template<class T>
KernelFit1D<T>::KernelFit1D(const std::vector<T> &x, const std::vector<T> &y,
	const T &bandwidth){
		
	// Constructor for the KernelFit1D object. Save the `x`, `y` data
	// and set an initial value for the `bandwidth`.
	
	if ( x.empty() || y.empty() )
		throw KernelFitError("From KernelFit1D::KernelFit1D(), "
			"one or both input vectors are empty!");
	
	if ( x.size() != y.size() )
		throw KernelFitError("From KernelFit1D::KernelFit1D(), input vectors "
			"must be equal in length!");
	
	if ( bandwidth <= 0.0 )
		throw KernelFitError("From KernelFit1D::KernelFit1D(), the bandwidth "
			"must be greater than zero!");
	
	_x = x;
	_y = y;
	
	_bandwidth = bandwidth;
	
}

template<class T>
std::vector<T> KernelFit1D<T>::Solve(const std::vector<T> &x){
	//
	// solve for the smooth profile through the data at all `x`
	//
	
	if ( x.empty() )
		throw KernelFitError("From KernelFit1D::Solve(), `x` was empty!");
	
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
	const std::vector<T> &z, const T &bandwidth){
		
	// Constructor for the KernelFit1D object. Save the `x`, `y` data
	// and set an initial value for the `bandwidth`.
	
	if ( x.empty() || y.empty() || z.empty() )
		throw KernelFitError("From KernelFit2D::KernelFit2D(), one or more "
			"input vectors were empty!");
	
	if ( x.size() != y.size() || x.size() != z.size() )
		throw KernelFitError("From KernelFit2D::KernelFit2D(), input vectors "
			"must be equal in length!");
	
	if ( bandwidth <= 0.0 )
		throw KernelFitError("From KernelFit2D::KernelFit2D(), the bandwidth "
			"must be greater than zero!");
	
	_x = x;
	_y = y;
	_z = z;
	
	_bandwidth = bandwidth;
	
}

template<class T>
std::vector< std::vector<T> > KernelFit2D<T>::Solve(const std::vector<T> &x,
	const std::vector<T> &y){
	
	//
	// solve for the smooth surface through the data at all (x, y)
	//
	
	if ( x.empty() || y.empty() )
		throw KernelFitError("From KernelFit2D::Solve(), one or both of "
			"`x` and `y` were empty!");

	// initialize f[x][y] to zeros with proper dimensions
	std::vector< std::vector<T> > f( x.size(), std::vector<T>( y.size(), 0.0));
	
	// omp_set_num_threads() should be called prior to here!
	#pragma omp parallel for shared(f)
	for (std::size_t i = 0; i < x.size(); i++)
	for (std::size_t j = 0; j < y.size(); j++){
		
		T sum = 0.0;
		
		for (std::size_t k = 0; k < _x.size(); k++){
			
			T sep    = sqrt( pow(x[i] - _x[k], 2.0) + pow(y[j] - _y[k], 2.0) );
			T W      = Kernel(sep);
			f[i][j] += W * _z[k];
			sum     += W;
		}
		
		f[i][j] /= sum;
	}
	
	return f;
}

template class KernelFit1D<double>;
template class KernelFit2D<double>;