#! /usr/bin/env python

from dolfin import *
from numpy import *
import distribution_diffusion as mod

# some constants
T = 5.0
dt = 0.1
dim = 2
n_ele = 1
N = 2
A_cond = 1e-11
A_pert = 1e-4
G = 1.0
D_x = 1e-7
weighted_abscissa_0 = Expression('0.5*G*x[0]', G=G)
weight_0 = Expression('0.5')

# Create files for storing solution
files = [File("results_1/mean.pvd"),
         File("results_1/std_dev.pvd")]

mod.main(T, dt, dim, n_ele, N, A_cond, A_pert, D_x, weighted_abscissa_0, weight_0, files)
