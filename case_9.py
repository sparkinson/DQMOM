#! /usr/bin/env python

from dolfin import *
from numpy import *
import distribution_diffusion as mod

def main():

    # some constants
    T = 5.0
    dt = 0.1
    dim = 2
    n_ele = 21
    N = 2
    A_cond = 1e-11
    A_pert = 5e-3
    D_x = 1e-3
    minField = 0.5
    abscissa_0 = [Expression('1.0+(floor(x[0]+0.5)+floor(x[1]+0.5))*1.0'),
                  Expression('0.5+(floor(x[0]+0.5)+floor(x[1]+0.5))*0.5')]
    weight_0 = [Expression('minField+1.0*(floor(x[0]+0.5)+floor(x[1]+0.5))', minField=minField), 
                Expression('minField+1.0*(floor(x[0]+0.5)+floor(x[1]+0.5))', minField=minField)]

    # import pdb; pdb.set_trace()

    mod.main(T, dt, dim, n_ele, N, A_cond, A_pert, D_x, abscissa_0, weight_0, 'results_9')

main()
