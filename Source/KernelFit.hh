// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// KernelFit/Source/KernelFit.hh

// This header file contains the template for the KernelFit objects.


#ifndef _KERNELFIT_HH_
#define _KERNELFIT_HH_

#include <string>
#include <vector>
#include <exception>
#include <cmath>

template<class T>
class KernelFit1D {

public:

	KernelFit1D(){}
	KernelFit1D(const std::vector<T> &x, const std::vector<T> &y,
		const double bandwidth);
	
	T Dispersion();
	
	void SetBandwidth(const double multiple);
	
	// kernel function is `inlined`
	T Kernel(const T x){
		return exp( -x * x / (2 * _bandwidth * _bandwidth));
	}
	
	std::vector<T> Solve(const std::vector<T> &x);

protected:
	
	double _bandwidth;
	std::vector<T> _x, _y;

};

template<class T>
class KernelFit2D {

public:

	KernelFit2D(){}
	KernelFit2D(const std::vector<T> &x, const std::vector<T> &y,
		const std::vector<T> &z, const double bandwidth);
	
	T Dispersion();
	
	void SetBandwidth(const double multiple);
	
	// kernel function is `inlined`
	T Kernel(const T r){
		return exp( -r * r / (2 * _bandwidth * _bandwidth));
	}
	
	std::vector<T> Solve(const std::vector<T> &x, const std::vector<T> &y);

protected:
	
	double _bandwidth;
	std::vector<T> _x, _y, _z;

};

#endif