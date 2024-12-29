import numpy as np
import matplotlib.pyplot as plt

#Variables
D = 0.540*24 # Glutamte diffusion within the biofilm 1/h
DE = 0.9*24 # Glutamate diffusion between biofilm and exterior 1/h
GE = 30 #Glutamate concentration in the external medium mMole
Ginmax = 20 #mMole Concentration in cell

DK = 0.497*24# mm^2/hour
DKE = 4.97*24# mm^2/hour
KE = 8
 
deltaK = 1#mV/mM
deltaL = 60#mV/mM

gK = 180 # 1/hour
gL = 1.2# 1/hour
VK0 = -380#-380 #mV
VL0 = -156#-156 #mV
Vth = -150#mMole
F = 5.6 #mMole/mV

alpha = 5 #OPENING RATE
beta = 2.5
gamma_G = 10#Ginmax*1.125 #GLUTAMINE PRODUCTION RATE
gamma_K = 1#.025
delta_G = 10
G_u = 18
V_low = -180
gamma_V = 20
eta_V = 20
eta_K = 30
PotassiumSetAmount = 300
delta_grow = 0.005
alphaT = 20
gammaT = 10
gT = .3
V0T = -170
r_b = 0.1

L = .8#.08#.06#.04 #/Lbar #Depth of Biofilm in mm #
Lout = .5#/Lbar#Diffusion layer

# Simulation Parameters
N = 100
dx = 1/N
x = np.linspace(0,1,N+1)
T = 20#/tbar#Terminal Time
dt = .001#/tbar#Time Step
maxiters = round(T/dt)
 
G = np.ones((N+1,1))*GE
GrowMode = np.ones((N+1,1))
K = np.ones((N+1,1))*KE
K_acc = np.ones((N+1,1))*KE
Kin = np.ones((N+1,1))*PotassiumSetAmount
V = np.ones((N+1,1))*Vth
nK = np.ones((N+1,1))*.1
Gin = np.ones((N+1,1))*Ginmax
ChangeGin = np.ones((N+1,1))*0
S = np.ones((N+1,1))*0
VK = np.ones((N+1,1))*VK0
VL = np.ones((N+1,1))*VL0
ThT = np.ones((N+1,1))*0
APG = np.ones((N+1,1))*0
Uin = np.ones((N+1,1))*0
 
U = np.zeros((N+1,1))
dL = np.zeros((maxiters,1))
dL2 = np.zeros((maxiters,1))
MLGdiff = np.zeros((N+1,N+1)) #Matrices for Glutamate
MRGdiff = np.zeros((N+1,N+1))
bGdiff = np.zeros((N+1,1))
MLKdiff = np.zeros((N+1,N+1)) #Matrices for Potassium
MRKdiff = np.zeros((N+1,N+1))
bKdiff = np.zeros((N+1,1))
 
 
#Set up Matrices
MLGdiff[0,0] = 1/dx
MLGdiff[0,1] = -1/dx
MLKdiff[0,0] = 1/dx
MLKdiff[0,1] = -1/dx

#Plotting utils
fig = plt.figure()
plt.ion() 
ax = fig.add_subplot(111) 
line1, = ax.plot(x*L,V, 'b-')
plt.ylim([-300,-100])
plt.show()
plot_freq = 100 # steps

for iters in range(maxiters):
    #Set Up Matrices, which depend on L, which changes!!!
    for matrix_iter in range(N-1):
        i = matrix_iter + 1
        MLGdiff[i,i] = 1 + 2*D*dt/np.pow(dx,2)/2/(np.pow(L,2))
        MLGdiff[i,i-1] = -D*dt/np.pow(dx,2)/2/(np.pow(L,2))
        MLGdiff[i,i+1] = -D*dt/np.pow(dx,2)/2/(np.pow(L,2))
   
        MRGdiff[i,i] = 1 + -2*D*dt/np.pow(dx,2)/2/(np.pow(L,2))
        MRGdiff[i,i-1] = D*dt/np.pow(dx,2)/2/(np.pow(L,2))
        MRGdiff[i,i+1] = D*dt/np.pow(dx,2)/2/(np.pow(L,2))
       
        MLKdiff[i,i] = 1 + 2*DK*dt/np.pow(dx,2)/2/(np.pow(L,2))
        MLKdiff[i,i-1] = -DK*dt/np.pow(dx,2)/2/(np.pow(L,2))
        MLKdiff[i,i+1] = -DK*dt/np.pow(dx,2)/2/(np.pow(L,2))
   
        MRKdiff[i,i] = 1 + -2*DK*dt/np.pow(dx,2)/2/(np.pow(L,2))
        MRKdiff[i,i-1] = DK*dt/np.pow(dx,2)/2/(np.pow(L,2))
        MRKdiff[i,i+1] = DK*dt/np.pow(dx,2)/2/(np.pow(L,2))
  
    MLGdiff[N,N] = D/L/dx+DE/Lout
    MLGdiff[N,N-1] = -D/L/dx
    bGdiff[N] = DE/Lout*GE
 
    MLKdiff[N,N] = (DK/L/dx+DKE/Lout)
    MLKdiff[N,N-1] = -DK/L/dx
    bKdiff[N] = DKE/Lout*KE

    # Glutamate in Surrounding Area
    CellIntake =  delta_G*(1/(1+np.exp((V-Vth)))*G*(Ginmax-Gin))
 
    b1 = MRGdiff@G -CellIntake*dt
    bGdiff[1:N-1] = b1[1:N-1]
    G = np.linalg.solve(MLGdiff,bGdiff)
    G = np.clip(G, 0, 10000)
   
   # Glutamate Consuption
    GrowMode = (Gin/(Gin+G_u))/(Gin/(Gin+G_u)+eta_V*((np.tanh((V/V_low-1)*gamma_V))+1))
    Growth = Gin*(GrowMode+r_b)
 
    # Glutamate inside cells
    ChangeGin = (CellIntake-gamma_G*Growth)
    ChangeGin[1:] = ChangeGin[1:] - (Uin[1:]*Gin[1:]-Uin[0:-1]*Gin[0:-1])/dx/L
    Gin  = Gin + ChangeGin*dt
    Gin = np.clip(Gin,0,Ginmax)
 
    # Potassium Gating
    changenK = (np.pow(1-Gin/Ginmax,2)/(np.pow((1-.5),2)+np.pow(1-Gin/Ginmax,2))*alpha)*(1-nK) - nK*beta
    changenK[1:] -= (Uin[1:]*nK[1:]-Uin[:-1]*nK[:-1])/dx/L

    nK = nK + changenK*dt
    nK = np.clip(nK,0,1)
 
    # Potassium inside cells
    changeK = F*gK*np.pow(nK,4)*(V-VK) + F*gL*(V-VL) - np.clip(gamma_K*K*(300-Kin),-10000,0)
    ChangeKin = -changeK
    ChangeKin[1:] = ChangeKin[1:] - (Uin[1:]*Kin[1:]-Uin[:-1]*Kin[:-1])/dx/L
    Kin = Kin + ChangeKin*dt
    Kin = np.clip(Kin,0,10000)
   
    # Calculate Ambient Potassium
    b1 = MRKdiff@K + changeK*dt
    bKdiff[1:-1] = b1[1:-1]
    K = np.linalg.solve(MLKdiff,bKdiff)
 
    # Compute K set point for cells
    ChangemeanK = eta_K*(K-K_acc)
    K_acc = K_acc + ChangemeanK*dt
 
    #Calculate Potential
    ChangeV =  -gK*np.pow(nK,4)*(V-VK) - gL*(V-VL)
    ChangeV[1:] -= (Uin[1:]*V[1:]-Uin[:-1]*V[:-1])/dx/L
    V = V + ChangeV*dt
    V = np.clip(V,-10000,0)       
    VK = VK0 + deltaK*K
    VL = VL0 + deltaL*(K-K_acc)
 
    #Calculate luminescence
    ChangeAPG = (.5*(K) - 1*APG)
    ChangeAPG[1:] = ChangeAPG[1:] - (Uin[1:]*APG[1:]-Uin[-1]*APG[:-1])/dx/L
    APG = APG + ChangeAPG*dt
   
    ChangeThT = (alphaT/(1+np.exp(gT*(V-V0T))) -gammaT*ThT)
    ChangeThT[1:] = ChangeThT[1:] - (Uin[1:]*ThT[1:]-Uin[:-1]*ThT[:-1])/dx/L
    ThT = ThT + ChangeThT*dt

    # Compute Change in biofilm length from growth
    Uin = np.cumsum(Gin*GrowMode/N*L)*delta_grow
    Uin = np.expand_dims(Uin,1)
    U = Uin[-1,0]
    L = L+U*dt

    if iters%plot_freq ==0:
        x_plot = x*L
        line1.set_xdata(x_plot)
        line1.set_ydata(V)
        plt.xlim(x_plot[0], x_plot[-1])
        plt.ylim([-300, -100])
        fig.canvas.draw() 
        fig.canvas.flush_events() 
 
    dL[iters] = np.mean(V)
    dL2[iters] = np.mean(K)

iter_num = np.linspace(1,maxiters,maxiters)*dt
 
plt.figure()
plt.plot(iter_num,dL, color='steelblue')
plt.title('Mean Voltage')
plt.xlabel('Time')
plt.ylabel('Voltage')
plt.savefig('Voltage')
 
plt.figure()
plt.plot(iter_num,dL2, color='darkorange')
plt.title('Potassium')
plt.xlabel('Time')
plt.ylabel('Potassium')
plt.savefig('Potassium')

 