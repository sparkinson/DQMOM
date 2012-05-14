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

def calculateDiagnostics(quads, N, A_cond, A_pert, files, FS, mesh, D_x, n_pert):

      # calculate abscissa's
      for i, q in enumerate(quads):
            q.abscissa.vector()[:] = (q.weighted_abscissa_1.vector().array()/
                                      q.weight_1.vector().array())  

      print 'constructing linear system'
            
      # get gradient squared sum
      g_sum_square = getGradientSumSquare(quads, mesh)

      # construct linear system matrices
      A, A_3, C = constructMatrices(quads, g_sum_square, D_x)

      print 'solving linear system'   

      # calculate source terms
      for v in range(mesh.num_vertices()):

            # check for ill-conditioned matrix
            if (1/linalg.cond(A[:,:,v])) < A_cond:
                  print 'ill conditioned matrix found' 

                  abscissa = []
                  for j in range(N):
                        abscissa.append(quads[j].abscissa.vector().array()[v])
                       
                  # perturbate abscissa to generate well-conditioned matrix
                  counter = 0
                  while (1/linalg.cond(A)) < A_cond:

                        counter += 1
                        if counter > 10:
                              print '...struggling to find well conditioned matrix' 
                              counter = 0

                        abscissa = abscissa + (random.random_sample((N,))*2. -1.)*A_pert
                        A[:,:,v] = reconstructMatrix(abscissa)
                        n_pert += 1
            
            # solve linear system
            sources = linalg.solve(A[:,:,v], dot(A_3[:,:,v], C[:,v]))
            
            for i, q in enumerate(quads):
                  q.weight_S.vector()[v] = sources[i]
                  q.weighted_abscissa_S.vector()[v] = sources[i+N]
            
      print 'calculating statistics'   

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
            
      print 'writing to file'   

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

def getGradientSumSquare(quads, mesh):
      dim = mesh.geometry().dim()

      # calculate sum of gradients squared using FE approximation of gradient at vertices
      g_sum_square = []
      for q in quads:
            if dim == 1:
                  g_abscissa_d = project(grad(q.abscissa), 
                                         FunctionSpace(mesh, 'Lagrange', degree))
            if dim == 2:
                  g_abscissa_d = project(grad(q.abscissa), 
                                         VectorFunctionSpace(mesh, 'Lagrange', degree))
            if dim == 3:
                  g_abscissa_d = project(grad(q.abscissa), 
                                         TensorFunctionSpace(mesh, 'Lagrange', degree))
            g_abscissa_d_split = g_abscissa_d.split(deepcopy=True) 
            try:
                  g_sum_square.append(g_abscissa_d_split[0].vector().array()**2)
                  for i in range(1, dim):
                        g_sum_square[-1] += g_abscissa_d_split[i].vector().array()**2
            except:
                  g_sum_square.append(g_abscissa_d.vector().array()**2) 

      return g_sum_square

def constructMatrices(quads, g_sum_square, D_x):
      N = len(quads)
      n = len(quads[0].abscissa.vector().array())
      
      # construct A matrices
      A_1 = zeros([2*N,N,n]) 
      A_2 = zeros([2*N,N,n]) 
      A_3 = zeros([2*N,N,n])
      for i in range(2*N):
            for j in range(N):
                  A_1[i,j] = (1-i)*quads[j].abscissa.vector().array()**(i)
                  A_2[i,j] = i*quads[j].abscissa.vector().array()**(i-1)
                  A_3[i,j] = i*(i-1)*quads[j].abscissa.vector().array()**(i-2)

      # construct matrix C
      C = zeros([N,n])
      for j in range(N):
            C[j] = quads[j].weight_1.vector().array()*D_x.vector().array()*g_sum_square[j]
      
      return append(A_1, A_2, 1), A_3, C

def reconstructMatrix(abscissa):
      N = len(abscissa)

      for i in range(2*N):
            for j in range(N):
                  A[i,j] = (1-i)*abscissa[j]**(i)
                  A[i,j+N] = i*abscissa[j]**(i-1)
      
      return A 

def main(T, dt, dim, n_ele, N, A_cond, A_pert, D_x, abscissa_0, weight_0, save_folder):

      # Create files for storing solution    
      files_a = [File(save_folder + "/weight_" + str(i) + ".pvd") for i in range(N)]   
      files_b = [File(save_folder + "/abscissa_" + str(i) + ".pvd") for i in range(N)]
      files_c = [File(save_folder + "/weighted_abscissa_" + str(i) + ".pvd") for i in range(N)]
      files_d = [File(save_folder + "/weight_S_" + str(i) + ".pvd") for i in range(N)]
      files_e = [File(save_folder + "/weighted_abscissa_S_" + str(i) + ".pvd") for i in range(N)]
      files_f = [File(save_folder + "/density.pvd"),
                 File(save_folder + "/mean.pvd"),
                 File(save_folder + "/std_dev.pvd"),
                 File(save_folder + "/skew.pvd")]
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
      quads, n_pert = calculateDiagnostics(quads, N, A_cond, A_pert, files, FS, mesh, D_x, n_pert)

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
            print 't = ', t
            print 'start solve'

            for q in quads:
                  q.weighted_abscissa_b = assemble(q.weighted_abscissa_L, tensor=q.weighted_abscissa_b) 
                  q.weight_b = assemble(q.weight_L, tensor=q.weight_b) 
                  solve(q.weighted_abscissa_A, q.weighted_abscissa.vector(), q.weighted_abscissa_b)
                  solve(q.weight_A, q.weight.vector(), q.weight_b)
                  q.weighted_abscissa_1.assign(q.weighted_abscissa)
                  q.weight_1.assign(q.weight)

            # Calculate diagnostics and increase time step
            quads, n_pert = calculateDiagnostics(quads, N, A_cond, A_pert, files, FS, mesh, D_x, n_pert)

            t += dt

      print 'number of ill-conditioned matrices = ', n_pert
