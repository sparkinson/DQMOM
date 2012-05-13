#! /usr/bin/env python

from dolfin import *
from numpy import *
import distribution_diffusion as mod

def main():

    # some constants
    T = 5.0
    dt = 0.1
    dim = 2
    n_ele = 40
    N = 2
    A_cond = 1e-11
    A_pert = 5e-4
    D_x = 1e-3
    G = 1.0
    minField = 1e-10
    abscissa_0 = [Expression('1.2*floor(x[0]+0.45)' + 
                             '+ 1.0*(1-ceil(x[0]-0.45))'),
                  Expression('1.0*floor(x[0]+0.45)' + 
                             '+ 0.8*(1-ceil(x[0]-0.45))')]
    weight_0 = [Expression('0.45*floor(x[0]+0.45)' + 
                           '+ 0.45*(1-ceil(x[0]-0.45))' +
                           '+ minField', minField=minField), 
                Expression('0.55*floor(x[0]+0.45)' + 
                           '+ 0.55*(1-ceil(x[0]-0.45))' +
                           '+ minField', minField=minField)]

    # Create files for storing solution    
    files_a = [File("results_4/weight_" + str(i) + ".pvd") for i in range(N)]   
    files_b = [File("results_4/abscissa_" + str(i) + ".pvd") for i in range(N)]
    files_c = [File("results_4/weighted_abscissa_" + str(i) + ".pvd") for i in range(N)]
    files_d = [File("results_4/weight_S_" + str(i) + ".pvd") for i in range(N)]
    files_e = [File("results_4/weighted_abscissa_S_" + str(i) + ".pvd") for i in range(N)]
    files_f = [File("results_4/density.pvd"),
               File("results_4/mean.pvd"),
               File("results_4/std_dev.pvd"),
               File("results_4/skew.pvd")]
    files = files_a + files_b + files_c + files_d + files_e + files_f

    # import pdb; pdb.set_trace()

    mod.main(T, dt, dim, n_ele, N, A_cond, A_pert, D_x, abscissa_0, weight_0, files)

main()
