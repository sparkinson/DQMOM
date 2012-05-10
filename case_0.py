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
A_cond = 1e-11
A_pert = 1e-4
G = 1.0
D_x = 1e-7
weighted_abscissa_0 = Expression('0.5*G*x[0]', G=G)
weight_0 = Expression('0.5')

# Create files for storing solution
files = [File("results/mean.pvd"),
         File("results/std_dev.pvd"),
         File("results/t_mean.pvd"),
         File("results/t_std_dev.pvd")]

def calculateError():
      for i, q in enumerate(quads):
            a = []; b = []
            for v in range(mesh.num_vertices()):
                  a.append(0.5)
                  b.append(G*mesh.coordinates()[v][0] + sin(3**i*pi/2) * sqrt(2*D_x.vector().array()[v]*G**2*t))
            # # print list(q.weight_1.vector().array())
            # # print a
            print (abs(q.weight_1.vector().array() - a)/q.weight_1.vector().array() * 100.0)
            # # print list(q.abscissa.vector().array())
            # # print b
            print (abs(q.abscissa.vector().array() - b)/q.abscissa.vector().array() * 100.0)
            # q.Print() 

def calculateTheoryStats():
      theoretical_moments = zeros([mesh.num_vertices(),2*N])
      for i in range(2*N):
            for j in range(N):
                  theoretical_moments[:,i] += 0.5 * (G*mesh.coordinates()[:,0] + 
                                                     sin(3**j*pi/2) * 
                                                     sqrt(2*D_x.vector().array()[:]*G**2*t))**i
      TheoryMean = Function(FS)
      TheoryStd = Function(FS)
      TheoryMean.vector()[:] = array(theoretical_moments[:,1])
      TheoryStd.vector()[:] = array(sqrt(theoretical_moments[:,2] - theoretical_moments[:,1]**2))
      
      files[2] << TheoryMean
      files[3] << TheoryStd

# mod.main(T, dt, dim, n_ele, N, A_cond, A_pert, D_x, weighted_abscissa_0, weight_0, files)

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

