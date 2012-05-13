from dolfin import *
from numpy import *

degree = 1
      
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

def calculateDiagnostics(quads, N, A_cond, A_pert, files, FS, mesh, D_x, dim, n_pert):

      # calculate abscissa's
      g_abscissa = []
      for i, q in enumerate(quads):
            q.abscissa.vector()[:] = (q.weighted_abscissa_1.vector().array()/q.weight_1.vector().array())

      # calculate sum of gradients squared using FE approximation of gradient at vertices
      g_sum_square = []
      for i, q in enumerate(quads):
            if dim == 1:
                  g_abscissa_d = project(grad(q.abscissa), FunctionSpace(mesh, 'Lagrange', degree))
            if dim == 2:
                  g_abscissa_d = project(grad(q.abscissa), VectorFunctionSpace(mesh, 'Lagrange', degree))
            if dim == 3:
                  g_abscissa_d = project(grad(q.abscissa), TensorFunctionSpace(mesh, 'Lagrange', degree))
            g_abscissa_d_split = g_abscissa_d.split(deepcopy=True) 
            try:
                  g_sum_square.append(g_abscissa_d_split[0].vector().array()**2)
                  for i in range(1, dim):
                        g_sum_square[-1] += g_abscissa_d_split[i].vector().array()**2
            except:
                  g_sum_square.append(g_abscissa_d.vector().array()**2)   
            
      print 'solving linear system'   

      Sources = [Function(FS) for i in range(2*N)]
      # calculate source terms
      for v in range(mesh.num_vertices()):
            abscissa = []
            weights = []
            for j in range(N):
                  abscissa.append(quads[j].abscissa.vector().array()[v])
                  weights.append(quads[j].weight_1.vector().array()[v])
                        
            A, A_3, C = constructMatrices(weights, abscissa, g_sum_square, v, N, D_x)
            counter = 0
            while (1/linalg.cond(A)) < A_cond:
                  n_pert += 1
                  abscissa = abscissa + (random.random_sample((N,))*2. -1.)*A_pert
                  A, A_3, C = constructMatrices(weights, abscissa, g_sum_square, v, N, D_x)
                  counter += 1
                  if counter > 10:
                        print 'Struggling to find well conditioned matrix' 
                        counter = 0
            
            sources = linalg.solve(A, dot(A_3, C))
            
            for i, q in enumerate(quads):
                  q.weight_S.vector()[v] = sources[i]
                  q.weighted_abscissa_S.vector()[v] = sources[i+N]
            
      print 'finished solving linear system'   

      # calculate mean and standard deviation 
      moments = zeros([mesh.num_vertices(),2*N])
      for i in range(2*N):
            for q in quads:
                  moments[:,i] += q.weight_1.vector().array()*q.abscissa.vector().array()**i
      Density = Function(FS)
      Mean = Function(FS)
      Std = Function(FS)
      Skew = Function(FS)
      d = Density.vector()[:] = array(moments[:,0])
      mu = Mean.vector()[:] = array(moments[:,1]/moments[:,0])
      var = moments[:,2]/moments[:,0] - (moments[:,1]/moments[:,0])**2
      var[var<0.0] = 0.0
      std = Std.vector()[:] = array(sqrt(var))
      Skew.vector()[:] = (moments[:,3]/moments[:,0] - 
              3.*mu*std**2 - 
              mu**3.) / std**3

      # Save to file
      for i in range(N):
            files[i] << quads[i].weight_1
            files[i+N] << quads[i].abscissa
            files[i+2*N] << quads[i].weighted_abscissa_1
            files[i+3*N] << quads[i].weight_S
            files[i+4*N] << quads[i].weighted_abscissa_S
      files[-4] << Density
      files[-3] << Mean
      files[-2] << Std  
      files[-1] << Skew 

      return quads, n_pert

def constructMatrices(weights, abscissa, g_sum_square, v, N, D_x):
      
      A = zeros([2*N,2*N])
      A_3 = zeros([2*N,N])
      C = zeros([N])
      for i in range(2*N):
            for j in range(N):
                  A[i,j] = (1-i)*abscissa[j]**(i)
                  A[i,j+N] = i*abscissa[j]**(i-1)
                  A_3[i,j] = i*(i-1)*abscissa[j]**(i-2)
                  C[j] = (weights[j]*D_x.vector().array()[v]*
                          g_sum_square[j][v])
      
      return A, A_3, C    

def main(T, dt, dim, n_ele, N, A_cond, A_pert, D_x, abscissa_0, weight_0):

      # Create files for storing solution    
      files_a = [File("results_6/weight_" + str(i) + ".pvd") for i in range(N)]   
      files_b = [File("results_6/abscissa_" + str(i) + ".pvd") for i in range(N)]
      files_c = [File("results_6/weighted_abscissa_" + str(i) + ".pvd") for i in range(N)]
      files_d = [File("results_6/weight_S_" + str(i) + ".pvd") for i in range(N)]
      files_e = [File("results_6/weighted_abscissa_S_" + str(i) + ".pvd") for i in range(N)]
      files_f = [File("results_6/density.pvd"),
                 File("results_6/mean.pvd"),
                 File("results_6/std_dev.pvd"),
                 File("results_6/skew.pvd")]
      files = files_a + files_b + files_c + files_d + files_e + files_f
      
      # Counter for number of perturbations
      n_pert = 0

      # Create mesh and define function space
      if dim == 1:
            mesh = UnitInterval(n_ele)
      if dim == 2:
            mesh = UnitSquare(n_ele, n_ele)
      if dim == 3:
            mesh = UnitCube(n_ele, n_ele, n_ele)

      FS = FunctionSpace(mesh, 'Lagrange', degree)
      D_x = project(Expression('D_x', D_x=D_x), FS)

      # Define quad classes
      quads = []
      for i in range(N):
            quads.append(quad(FS))
            quads[i].weighted_abscissa_1 = project(abscissa_0[i], FS) 
            quads[i].weight_1 = project(weight_0[i], FS)
            quads[i].weighted_abscissa_1.vector()[:] = (quads[i].weight_1.vector().array()*
                                                        quads[i].weighted_abscissa_1.vector().array())

      # Calculate diagnostics
      t = 0.0
      quads, n_pert = calculateDiagnostics(quads, N, A_cond, A_pert, files, FS, mesh, D_x, dim, n_pert)

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
            q.weighted_abscissa = Function(FS)
            q.weight = Function(FS)

      t = dt
      # Time step loop
      while t <= T:
            print t
            print 'start solve'

            for q in quads:
                  q.weighted_abscissa_b = assemble(q.weighted_abscissa_L, tensor=q.weighted_abscissa_b) 
                  q.weight_b = assemble(q.weight_L, tensor=q.weight_b) 
                  solve(q.weighted_abscissa_A, q.weighted_abscissa.vector(), q.weighted_abscissa_b)
                  solve(q.weight_A, q.weight.vector(), q.weight_b)
                  q.weighted_abscissa_1.assign(q.weighted_abscissa)
                  q.weight_1.assign(q.weight)

            print 'end solve'
            
            print 'start diagnostics'

            # Calculate diagnostics and increase time step
            quads, n_pert = calculateDiagnostics(quads, N, A_cond, A_pert, files, FS, mesh, D_x, dim, n_pert)

            print 'end diagnostics'

            t += dt

      print 'number of ill-conditioned matrices = ', n_pert
