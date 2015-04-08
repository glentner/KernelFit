// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// KernelFit/Source/KernelFit.hh

// This header file contains the template for the KernelFit objects.


#ifndef _KERNELFIT_HH_
#define _KERNELFIT_HH_

template<class T>
class KernelFit1D {

public:

	KernelFit1D(){}
	KernelFit1D(const std::vector<T> &x, const std::vector<T> &y,
		const double &bandwidth);
	
	Dispersion();
	SetBandwidth();
	Solve(const std::vector<T> &x);

protected:
	
	double _bandwidth;
	std::vector<T> _x, _y;

};

#endif