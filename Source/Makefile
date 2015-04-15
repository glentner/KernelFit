# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE file (GPLv3)
# KernelSmoother/Source/Makefile

# Makefile for building the KernelSmoother source (otherwise, just add the
# KernelSmoother.hh and KernelSmoother.cc files to your project as you like).

ccflags = -std=c++11 -fopenmp
inc     = .

# compile and run the test programs
Test: Run1D Run2D
	@echo
	@echo " The following output files have been created: "
	@echo
	@echo " Test/raw-1D.dat --> x, y positions of raw data "
	@echo " Test/fit-1D.dat --> x, y positions of fitted curve "
	@echo " Test/std-1D.dat --> x, y positions of standard deviations"
	@echo
	@echo " Test/raw-2D.dat --> x, y, z positions of raw data"
	@echo " Test/fit-2D.dat --> z matrix of values for fitted surface"
	@echo " Test/std-2D.dat --> matrix of values for standard deviations"
	@echo " Test/lattice-2D.dat --> x, y vectors for mesh-grid"
	@echo

Run1D: Test1D
	@echo
	@echo " Running the 1D test ... "
	@echo
	@Test/TestKernelFit1D

Run2D: Test2D
	@echo
	@echo " Running the 2D test ... "
	@echo
	@Test/TestKernelFit2D

Test1D: Test/TestKernelFit1D.cc
	g++ -g -o Test/TestKernelFit1D Test/TestKernelFit1D.cc KernelFit.cc \
	 	-I${inc} ${ccflags}

Test2D: Test/TestKernelFit2D.cc
	g++ -o Test/TestKernelFit2D Test/TestKernelFit2D.cc KernelFit.cc \
		-I${inc} ${ccflags}
