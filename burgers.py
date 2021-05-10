def solve_burgers(eta,profile): # Solve eq. (10) Bollati et al. 2021

   profile = profile * (s.indad+1) * s.betap / 2**(3/4) # Eq. (15) Bollati et al. 2021

   deta = eta[1]-eta[0]

   tf_th = 300 # time required to develop N-wave for betap = 1
   tf = tf_th/s.betap # time required to display N-wave for generic betap, Eq. (39) Rafikov 2002

   eta_min = -s.eta_tilde-np.sqrt(2*s.C*tf) - 3  # Eq. (18) Bollati et al. 2021
   eta_max = -s.eta_tilde+np.sqrt(2*s.C*tf) + 3

   # extend the profile domain due to profile broadening during evolution ....
   extr = eta[-1]
   while extr<eta_max:
      eta = np.append(eta,eta[-1]+deta)
      extr = eta[-1]
      profile = np.append(profile,0)
   extr = eta[0]
   while extr>eta_min:
      eta = np.insert(eta,0,eta[0]-deta)
      extr = eta[0]
      profile = np.insert(profile,0,0)

   Neta = len(eta) # number of centers
   a = eta[0] - deta/2
   b = eta[-1] + deta/2

   L  = b-a   # spatial grid length
   T  = tf    # final time

   # set time array

   x = np.zeros(Neta+1) # cells edges
   for i in range(0,Neta+1):
       x[i] = a+i*deta

   solution = [np.zeros((Neta),  dtype=float)]
   time = [0]

   # linear solution as initial condition
   solution[0] = profile

   # define flux vector
   F = np.zeros(Neta+1) #there is a flux at each cell edge!!!!!

   # define the flux function of Burgers equation
   def flux(u):
       return 0.5*u**2

   # define the Central difference numerical flux---> ritorna la media aritmetica dei flux sx e dx passati
   def CentralDifferenceFlux(uL,uR):
       # compute physical fluxes at left and right state
       FL = flux(uL)
       FR = flux(uR)
       return 0.5*(FL+FR)

   def GodunovNumericalFlux(uL,uR):
     # compute physical fluxes at left and right state
     FL = flux(uL)
     FR = flux(uR)
     # compute the shock speed
     s = 0.5*(uL + uR)
     # from Toro's book
     if (uL >= uR):
         if (s > 0.0):
             return FL
         else:
             return FR
     else:
         if (uL > 0.0):
             return FL
         elif (uR < 0.0):
             return FR
         else:
             return 0.0

   def NumericalFlux(uL,uR):
       # return CentralDifferenceFlux(uL,uR)
       return GodunovNumericalFlux(uL,uR)

   # time integrate
   lapsed_time = 0
   counter = 0
   while lapsed_time < tf:
   #while counter < 10:
      counter += 1
      dt = min(deta * s.CFL / (max(abs(solution[-1])) + 1e-8), 0.02)
      #print('max of sol = ', max(abs(solution[-1])), ', dt = ', dt)
      time.append(time[-1]+dt)
      lapsed_time += dt

      # compute the interior fluxes
      for i in range(1,Neta):
          uL = solution[-1][i-1]
          uR = solution[-1][i]
          F[i] = NumericalFlux(uL,uR)

      # compute the left boundary flux
      if solution[-1][0] < 0.0:
          uL = 2.0*solution[-1][0] - solution[-1][1]
      else:
          uL = solution[0][0]
      uR = solution[-1][0]
      F[0] = NumericalFlux(uL,uR)

      # compute the right boundary flux
      if solution[-1][Neta-1] > 0.0:
          uR = 2.0 * solution[-1][Neta-1] - solution[-1][Neta-2]
      else:
          uR = solution[0][Neta-1]
      uL = solution[-1][Neta-1]
      F[Neta] = NumericalFlux(uL,uR)

      solution.append(solution[-1][0:Neta] - dt / deta * (F[1:Neta+1] - F[0:Neta]))

   solution = np.array(solution).transpose()
   time = np.array(time)
   # plot Fig. 3 Bollati et al. 2021
   """
   plt.plot(eta, solution[:,0], label = "$t=t_0+$ ")#+str(0*dt))
   plt.plot(eta, solution[:,100], label = "$t=t_0+$ ")#+str(round(dt*Nt/10,2)))
   plt.plot(eta, solution[:,200], label = "$t=t_0+$ ")#+str(round(dt*Nt/6,2)))
   plt.plot(eta, solution[:,300], label = "$t=t_0+$ ")#+str(round(dt*Nt/3,2)))
   plt.plot(eta, solution[:,400], label = "$t=t_0+$ ")#+str(round(T,2)))
   plt.legend()
   plt.xlabel("$\eta$")
   plt.ylabel("$\chi(t,\eta)$")
   plt.grid(True)
   plt.title('$\chi$ "evolution" $r > r_p$')
   plt.show()
   """

   """
   print(np.shape(solution))
   plt.contourf(time[:3000], eta, solution[:,:3000], levels=np.arange(-30,30,0.05), cmap='RdBu')
   plt.colorbar()
   plt.show()
   """

   Nt = np.shape(solution)[1]

   solution_inner = np.zeros(solution.shape) # solution for r < Rp (and disc rotating counterclockwise)
   for i in range(Neta):
      for j in range(Nt):
         solution_inner[i,j] = solution[int(Neta-1-i),j]

   eta_inner = - eta[::-1]

   return time, eta, solution, eta_inner, solution_inner