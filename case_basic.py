#! /usr/bin/env python

from dolfin import *
from numpy import *
import distribution_diffusion as mod

def main():

    # some constants
    T = 0.1
    dt = 0.1
    dim = 1
    n_ele = 20
    N = 2
    A_cond = 1e-11
    A_pert = 5e-4
    D_x = 1e-3
    G = 1.0
    minField = 1e-1
    abscissa_0 = [Expression('1.5'), 
                  Expression('0.5')]
    weight_0 = [Expression('0.5*floor(x[0]+0.5)' +
                           '+ minField', minField=minField), 
                Expression('0.5*floor(x[0]+0.5)' +
                           '+ minField', minField=minField)]

    # Create files for storing solution    
    files_a = [File("results_basic/weight_" + str(i) + ".pvd") for i in range(N)]   
    files_b = [File("results_basic/abscissa_" + str(i) + ".pvd") for i in range(N)]
    files_c = [File("results_basic/weighted_abscissa_" + str(i) + ".pvd") for i in range(N)]
    files_d = [File("results_basic/density.pvd"),
               File("results_basic/mean.pvd"),
               File("results_basic/std_dev.pvd")]
    files = files_a + files_b + files_c + files_d

    # import pdb; pdb.set_trace()

    mod.main(T, dt, dim, n_ele, N, A_cond, A_pert, D_x, abscissa_0, weight_0, files)

main()
