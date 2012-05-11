#! /usr/bin/env python

from dolfin import *
from numpy import *
import distribution_diffusion as mod

# some constants
T = 100.0
dt = 0.1
dim = 1
n_ele = 5
N = 2
A_cond = 1e-8
A_pert = 5e-3
G = 1.0
D_x = 1e-7
weighted_abscissa_0 = [Expression('1.0', G=G),Expression('1.0', G=G)]
weight_0 = [Expression('0.5'),Expression('0.5')]

# Create files for storing solution    
files_a = [File("results/weight_" + str(i) + ".pvd") for i in range(N)]   
files_b = [File("results/abscissa_" + str(i) + ".pvd") for i in range(N)]
files_c = [File("results/density.pvd"),
           File("results/mean.pvd"),
           File("results/std_dev.pvd"),
files = files_a + files_b + files_c

mod.main(T, dt, dim, n_ele, N, A_cond, A_pert, D_x, weighted_abscissa_0, weight_0, files)
