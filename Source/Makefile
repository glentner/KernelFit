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
	@echo
	@echo " Test/raw-2D.dat --> x, y, z positions of raw data"
	@echo " Test/fit-2D.dat --> x, y, z positions of fitted surface"
	@echo
		
Run1D: Test1D
	@echo
	@echo " Running the 1D test ... "
	@echo
	Test/TestKernelFit1D
	
Run2D: Test2D
	@echo
	@echo " Running the 2D test ... "
	@echo
	Test/TestKernelFit2D
	
Test1D: Test/TestKernelFit1D.cc
	@echo
	@echo " Compiling TestKernelFit1D.cc ... "
	@echo
	g++ -g -o Test/TestKernelFit1D Test/TestKernelFit1D.cc KernelFit.cc \
	 	-I${inc} ${ccflags}
	
Test2D: Test/TestKernelFit2D.cc
	@echo
	@echo " Compiling TestKernelFit2D.cc ... "
	@echo
	g++ -o Test/TestKernelFit2D Test/TestKernelFit2D.cc KernelFit.cc \
		-I${inc} ${ccflags}
