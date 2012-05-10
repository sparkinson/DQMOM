#! /usr/bin/env python

from dolfin import *
from numpy import *
      
class quad:
      def __init__(self, FS):
            self.abscissa = Function(FS)
            self.weight = TrialFunction(FS)
            self.weighted_abscissa = TrialFunction(FS)
            self.weight_1 = None
            self.weighted_abscissa_1 = None
            self.weight_tf = TestFunction(FS)
            self.weighted_abscissa_tf = TestFunction(FS)
            self.weight_a = None
            self.weighted_abscissa_a = None
            self.weight_L = None
            self.weighted_abscissa_L = None
            self.weight_A = None
            self.weighted_abscissa_A = None
            self.weight_b = None
            self.weighted_abscissa_b = None
            self.weight_S = Function(FS)
            self.weighted_abscissa_S = Function(FS)

      def Print(self):
            print self.weight_1.vector().array(), self.abscissa.vector().array()

def calculateDiagnostics():
      g_abscissa = []
      for i, q in enumerate(quads):
            q.abscissa.vector()[:] = (q.weighted_abscissa_1.vector().array()/q.weight_1.vector().array())
            g_abscissa.append(project(grad(q.abscissa), FunctionSpace(mesh, 'Lagrange', 1)))

      for v in range(mesh.num_vertices()):
            A_abscissa = array([quads[0].abscissa.vector().array()[v],quads[1].abscissa.vector().array()[v]])
            if abs(A_abscissa[1]-A_abscissa[0])<1e-4:
                  if sign(A_abscissa[1]-A_abscissa[0])==0.0:
                        A_abscissa[1] = A_abscissa[1] + 1e-4
                  else:
                        A_abscissa[1] = (A_abscissa[1] + sign(A_abscissa[1]-A_abscissa[0])*
                                         (1e-4 - abs(A_abscissa[1]-A_abscissa[0])))

            A = zeros([2*N,2*N])
            A_3 = zeros([2*N,N])
            C = zeros([N])
            for i in range(2*N):
                  for j in range(N):
                        A[i,j] = (1-i)*A_abscissa[j]**(i)
                        A[i,j+N] = i*A_abscissa[j]**(i-1)
                        A_3[i,j] = i*(i-1)*A_abscissa[j]**(i-2)
                        C[j] = (quads[j].weight_1.vector().array()[v]*D_x.vector().array()[v]*
                                g_abscissa[j].vector().array()[v]**2)

            sources = linalg.solve(A, dot(A_3, C))
            
            for i, q in enumerate(quads):
                  q.weight_S = sources[i]
                  q.weighted_abscissa_S = sources[i+N]         

      # print 'mean is: ', moments[:,1]
      # # print 'var is:  ', moments[:,2] - moments[:,1]**2
      # print 'std is:  ', sqrt(moments[:,2] - moments[:,1]**2)

def calculateError():
      for i, q in enumerate(quads):
            a = []; b = []
            for v in range(mesh.num_vertices()):
                  a.append(0.5)
                  b.append(G*mesh.coordinates()[v][0] + sin(3**i*pi/2) * sqrt(2*D_x.vector().array()[v]*G**2*t))
            # # print list(q.weight_1.vector().array())
            # # print a
            print max((abs(q.weight_1.vector().array() - a)/q.weight_1.vector().array() * 100.0)[1:])
            # # print list(q.abscissa.vector().array())
            # # print b
            print max((abs(q.abscissa.vector().array() - b)/q.abscissa.vector().array() * 100.0)[1:])
            # q.Print()    

def calculateStats(FS):     
      moments = zeros([mesh.num_vertices(),2*N])
      for i in range(2*N):
            for q in quads:
                  moments[:,i] += q.weight_1.vector().array()*q.abscissa.vector().array()**i
      theoretical_moments = zeros([mesh.num_vertices(),2*N])
      for i in range(2*N):
            for j, q in enumerate(quads):
                  theoretical_moments[:,i] += 0.5 * (G*mesh.coordinates()[:,0] + 
                                                     sin(3**j*pi/2) * 
                                                     sqrt(2*D_x.vector().array()[:]*G**2*t))**i
      Mean = Function(FS)
      Std = Function(FS)
      TheoryMean = Function(FS)
      TheoryStd = Function(FS)
      Mean.vector()[:] = array(moments[:,1])
      TheoryMean.vector()[:] = array(theoretical_moments[:,1])
      Std.vector()[:] = array(sqrt(moments[:,2] - moments[:,1]**2))
      TheoryStd.vector()[:] = array(sqrt(theoretical_moments[:,2] - theoretical_moments[:,1]**2))

      # Save to file
      file1 << Mean
      file2 << Std
      file3 << TheoryMean
      file4 << TheoryStd

# some constants
G = 1.0
T = 1.0
dt = 0.001
N = 2
t = 0.0

# Create mesh and define function space
mesh = UnitInterval(30)
#mesh = UnitSquare(6, 4)
#mesh = UnitCube(6, 4, 5)
P1 = FunctionSpace(mesh, 'Lagrange', 1)

# Define initial values function
weighted_abscissa_0 = Expression('0.5*G*x[0]', G=G)
weight_0 = Expression('0.5')
D_x = project(Expression('D_x', D_x=1e-7), P1)

quads = []
for i in range(N):
      quads.append(quad(P1))
      quads[i].weighted_abscissa_1 = project(weighted_abscissa_0, P1) 
      quads[i].weight_1 = project(weight_0, P1)

calculateDiagnostics()

# Create files for storing solution
file1 = File("results/mean.pvd")
file2 = File("results/std_dev.pvd")
file3 = File("results/t_mean.pvd")
file4 = File("results/t_std_dev.pvd")

# Define variational problem
for q in quads:
      q.weighted_abscissa_a = (q.weighted_abscissa*q.weighted_abscissa_tf*dx + 
                               dt*inner(nabla_grad(q.weighted_abscissa), nabla_grad(D_x*q.weighted_abscissa_tf))*dx)
      q.weighted_abscissa_L = (q.weighted_abscissa_1 + dt*q.weighted_abscissa_S)*q.weighted_abscissa_tf*dx
      q.weight_a = (q.weight*q.weight_tf*dx + 
                    dt*inner(nabla_grad(q.weight), nabla_grad(D_x*q.weight_tf))*dx)
      q.weight_L = (q.weight_1 + dt*q.weight_S)*q.weight_tf*dx
      q.weighted_abscissa_A = assemble(q.weighted_abscissa_a) # assemble only once, before the time stepping
      q.weight_A = assemble(q.weight_a)
      q.weighted_abscissa = Function(P1)
      q.weight = Function(P1)

t = dt
# Time step loop
while t <= T:
      for q in quads:
            q.weighted_abscissa_b = assemble(q.weighted_abscissa_L, tensor=q.weighted_abscissa_b) 
            q.weight_b = assemble(q.weight_L, tensor=q.weight_b) 
            solve(q.weighted_abscissa_A, q.weighted_abscissa.vector(), q.weighted_abscissa_b)
            solve(q.weight_A, q.weight.vector(), q.weight_b)
            q.weighted_abscissa_1.assign(q.weighted_abscissa)
            q.weight_1.assign(q.weight)
            
      calculateDiagnostics()
      calculateStats(P1)
      t += dt

calculateError()
