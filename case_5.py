#! /usr/bin/env python

from dolfin import *
from numpy import *
import distribution_diffusion as mod

def main():

    # some constants
    T = 1.0
    dt = 0.1
    dim = 1
    n_ele = 10
    N = 2
    A_cond = 1e-11
    A_pert = 5e-4
    D_x = 1e-4
    G = 10.0
    minField = 1e-10
    abscissa_0 = [Expression('G*x[0]', G=G), 
                  Expression('G*x[0]-1.0', G=G)]
    weight_0 = [Expression('0.5'), 
                Expression('0.5')]

    # Create files for storing solution    
    files_a = [File("results_5/weight_" + str(i) + ".pvd") for i in range(N)]   
    files_b = [File("results_5/abscissa_" + str(i) + ".pvd") for i in range(N)]
    files_c = [File("results_5/weighted_abscissa_" + str(i) + ".pvd") for i in range(N)]
    files_d = [File("results_5/density.pvd"),File("results_5/mean.pvd"),File("results_5/std_dev.pvd")]
    files_e = [File("results_5/t_density.pvd"),File("results_5/t_mean.pvd"),File("results_5/t_std_dev.pvd")]
    files = files_a + files_b + files_c + files_d + files_e

    # import pdb; pdb.set_trace()

    mod.main(T, dt, dim, n_ele, N, A_cond, A_pert, D_x, abscissa_0, weight_0, files[:-3])

    def calculateTheoryStats():
        theoretical_moments = zeros([mesh.num_vertices(),2*N])
        for i in range(2*N):
            for j in range(N):
                theoretical_moments[:,i] += 0.5 * (G*mesh.coordinates()[:,0] - 1.0*j +
                                                   sin(3**j*pi/2) * 
                                                   sqrt(2*D_x.vector().array()[:]*G**2*t))**i
        TheoryDensity = Function(FS)
        TheoryMean = Function(FS)
        TheoryStd = Function(FS)
        TheoryDensity.vector()[:] = array(theoretical_moments[:,0])
        TheoryMean.vector()[:] = array(theoretical_moments[:,1]/theoretical_moments[:,0])
        var = (theoretical_moments[:,2]/theoretical_moments[:,0]
               - (theoretical_moments[:,1]/theoretical_moments[:,0])**2)
        var[var<0.0] = 0.0
        TheoryStd.vector()[:] = array(sqrt(var))
        
        files[-3] << TheoryDensity
        files[-2] << TheoryMean
        files[-1] << TheoryStd

    # Create mesh and define function space
    if dim == 1:
        mesh = UnitInterval(n_ele)
    if dim == 2:
        mesh = UnitSquare(n_ele, n_ele)
    if dim ==3:
        mesh = UnitCube(n_ele, n_ele, n_ele)
    
    FS = FunctionSpace(mesh, 'Lagrange', 1)
    D_x = project(Expression('D_x', D_x=D_x), FS)
    
    t = 0
    # Time step loop for theory
    while t <= T:
        calculateTheoryStats()
        t += dt

main()
