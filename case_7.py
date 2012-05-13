#! /usr/bin/env python

from dolfin import *
from numpy import *
import distribution_diffusion as mod

def main():

    # some constants
    T = 5.0
    dt = 0.1
    dim = 1
    n_ele = 41
    N = 2
    A_cond = 1e-11
    A_pert = 5e-4
    D_x = 1e-3
    G = 1.0
    minField = 1e-10
    abscissa_0 = [Expression('0.6*floor(x[0]+0.5) + 0.4*(1.0-floor(x[0]+0.5))'),
                  Expression('0.4*floor(x[0]+0.5) + 0.2*(1.0-floor(x[0]+0.5))')]
    weight_0 = [Expression('0.25*floor(x[0]+0.5) + 0.25*(1.0-floor(x[0]+0.5))'), 
                Expression('0.25*floor(x[0]+0.5) + 0.25*(1.0-floor(x[0]+0.5))')]

    # import pdb; pdb.set_trace()

    mod.main(T, dt, dim, n_ele, N, A_cond, A_pert, D_x, abscissa_0, weight_0)

main()
