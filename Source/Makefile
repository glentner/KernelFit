# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE file (GPLv3)
# KernelSmoother/Source/Makefile

# Makefile for building the KernelSmoother source (otherwise, just add the
# KernelSmoother.hh and KernelSmoother.cc files to your project as you like).

# compile the test program and link to the object
Test: Objects
	${compile} ${debug} -o ${exe} Main.o -I${inc} -L. -L${lib}