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
		const T &bandwidth);
	
	// kernel function used by default
	T Kernel(const T &x){return exp( -x * x / (2 * _bandwidth * _bandwidth));}
	
	// solve for smooth curve through data
	std::vector<T> Solve(const std::vector<T> &x);

protected:
	
	T _bandwidth;
	std::vector<T> _x, _y;

};

template<class T>
class KernelFit2D {

public:

	KernelFit2D(){}
	KernelFit2D(const std::vector<T> &x, const std::vector<T> &y,
		const std::vector<T> &z, const T &bandwidth);
	
	// kernel function used by default
	T Kernel(const T &r){return exp( -r * r / (2 * _bandwidth * _bandwidth));}
	
	// solve for the smooth surface through the data
	std::vector< std::vector<T> > Solve(const std::vector<T> &x, 
		const std::vector<T> &y);

protected:
	
	T _bandwidth;
	std::vector<T> _x, _y, _z;

};

// base exception class for KernelFit objects
class KernelException : public std::exception {	

public:
	
	explicit KernelException(const std::string& msg): _msg(msg){ }
	virtual ~KernelException() throw() { }
	virtual const char* what() const throw(){ return _msg.c_str(); }

protected:

	std::string _msg;
};

// exception thrown by KernelFit objects
class KernelFitError : public KernelException {
public:
	
	KernelFitError(const std::string& msg): KernelException(
		" --> KernelFitError: " + msg){ }
};

#endif