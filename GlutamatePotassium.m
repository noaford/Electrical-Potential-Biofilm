%clear all
%Variables
alpha = .1; %Glutamate dehydrogenation coefficient 1/mMole/h
deltaGrowth = 100; % Glutamate consumption rate 1/mMole/h
deltaGDecay = .1;
D = .4; % Glutamte diffusion within the biofilm 1/h
DE = .6; % Glutamate diffusion between biofilm and exterior 1/h
DA = 40; %Ammonium Diffusion in interior 1/h
DAE = 60; %Ammonium Diffusion in exterior 1/h
GE = 30; %Glutamate concentration in the external medium mMole
AE = 0; %Ammonium concentration in the external medium mMole
betaH = 50/1000; % GMaximal activation rate of GDH mMole/h
gammaH = 7.5; % deactivation rate of GDH 1/h
KH = 1; %GDP activation threshold mMole
n = 10; %Hill coefficent for GDH activation
betar = .14; % Expression coefficent of ribosomal/housekeeping proteins 1/mMole/h
gammar = 2; %Degradation rate of ribosomal/housekeeping proteins 1/h
eta = 100; %Population growth rate coefficient 1/mMole
Kk = .85; % Glutamate threshold for carrying capacity mMole
%m = 12; %Hill coefficient for carrying capacity
lambdai = 0; %Expansion rate of interior cells 1/h
lambdap = .032; %expansion rate of peripheral cells
Gln = 1; %Concentration of glutamine in the medium mMole
Kalpha = 5E-8; %Glutamine inhibition threshold on GDH activity mMole
Kgamma = 5E-2; %Glutamine inhibition threshold on GS activity mMole
%alpha0 = .03; %Rate of ammonium entering the biofilm from the external medium mMole/h
beta0 = 1.5E-6; %Expression rate of gDH from the additional copy of the gene mMole/h
lambdaH2O2 = 5; %Death rate due to hydrogen peroxide 1/h
KA = 2;%Monod Term for producing ammonia

gK = 30*60; % 1/hour
gL = .2*60;% 1/hour
VK0 = -380; %mV
VL0 = -156; %mV
Sth = 40/1000; %mMole
Vth = -150;%mMole
alpha0 = 2*60; %1/hour
beta = 1.3*60; %1/hour
m = 1;
F = 5.6; %mMole/mV
sigma = .2;%mV
deltaK = 1;%mV/mM
deltaL = 8;%mV/mM
gammaS = .1*60;%1/hour
gammaE = 10*60;%1/hour
gammaT = 4*60;%1/hour
alphas = 1*60/1000;%mmole/(hour mV)
alphat = 1*60/1000;%mmole/(hour mV)
DK = 13.8E-6*10^2*60*60;% mm^2/hour
DKE = 13.8E-6*10^2*60*60;% mm^2/hour
KE = 0;

Lbar = .1;
tbar = GE/(alpha+deltaGrowth+deltaGDecay);
Gbar = GE;
Abar = 1E-4;
Kbar = 1E-4;

%Scale terms
D = D/Lbar^2;
DE = DE/Lbar^2;
DA = DA/Lbar^2;
DAE = DAE/Lbar^2; 
DK = DK/Lbar^2;
DKE = DKE/Lbar^2;
GE = 1;

L = 10/Lbar; %Depth of Biofilm in mm
Lout = .01/Lbar;%Diffusion layer


N = 100;
dx = 1/N;
x = linspace(0,L,N+1);
T = .4/tbar;%Terminal Time
dt = .00001/tbar;%Time Step
maxiters = round(T/dt);

G = ones(N+1,1)*1/Gbar;
A = ones(N+1,1)*0.004/Abar;
K = ones(N+1,1).*heaviside(-x'+L/40)/Kbar/10;
V = ones(N+1,1)*(Vth*0-155);
nK = ones(N+1,1)*0;
S = ones(N+1,1)*0;
VK = ones(N+1,1)*VK0;
VL = ones(N+1,1)*VL0;
Gm = ones(N+1,1)*1/Gbar;
Am = ones(N+1,1)*0.004/Abar;
Km = ones(N+1,1).*heaviside(x'-L/2)/Kbar;
growth = zeros(N+1,1);


U = zeros(N+1,1);
dL = zeros(maxiters,1);
MLGdiff = zeros(N+1); %Matrices for Glutamate
MRGdiff = zeros(N+1);
bGdiff = zeros(N+1,1);
MLAdiff = zeros(N+1); %Matrices for Ammonium
MRAdiff = zeros(N+1);
bAdiff = zeros(N+1,1);
MLKdiff = zeros(N+1); %Matrices for Potassium
MRKdiff = zeros(N+1);
bKdiff = zeros(N+1,1);

%Set up Matrices
MLGdiff(1,1) = 1/dx;
MLGdiff(1,2) = -1/dx;
MLAdiff(1,1) = 1/dx;
MLAdiff(1,2) = -1/dx;
MLKdiff(1,1) = 1/dx;
MLKdiff(1,2) = -1/dx;


for iters = 1:maxiters
    %Set Up Matrices, which depend on L, which changes!!!
    for i = 2:N
        MLGdiff(i,i) = 1 + 2*D*dt/dx^2/2/(L^2);
        MLGdiff(i,i-1) = -D*dt/dx^2/2/(L^2);
        MLGdiff(i,i+1) = -D*dt/dx^2/2/(L^2);
    
        MRGdiff(i,i) = 1 + -2*D*dt/dx^2/2/(L^2);
        MRGdiff(i,i-1) = D*dt/dx^2/2/(L^2);
        MRGdiff(i,i+1) = D*dt/dx^2/2/(L^2);
    
        MLAdiff(i,i) = 1 + 2*DA*dt/dx^2/2/(L^2);
        MLAdiff(i,i-1) = -DA*dt/dx^2/2/(L^2);
        MLAdiff(i,i+1) = -DA*dt/dx^2/2/(L^2);
    
        MRAdiff(i,i) = 1 + -2*DK*dt/dx^2/2/(L^2);
        MRAdiff(i,i-1) = DK*dt/dx^2/2/(L^2);
        MRAdiff(i,i+1) = DK*dt/dx^2/2/(L^2);
        
        MLKdiff(i,i) = 1 + 2*DK*dt/dx^2/2/(L^2);
        MLKdiff(i,i-1) = -DK*dt/dx^2/2/(L^2);
        MLKdiff(i,i+1) = -DK*dt/dx^2/2/(L^2);
    
        MRKdiff(i,i) = 1 + -2*DK*dt/dx^2/2/(L^2);
        MRKdiff(i,i-1) = DK*dt/dx^2/2/(L^2);
        MRKdiff(i,i+1) = DK*dt/dx^2/2/(L^2);
    end

    MLGdiff(N+1,N+1) = D/L/dx+DE/Lout;
    MLGdiff(N+1,N) = -D/L/dx;
    bGdiff(N+1) = DE/Lout*GE;

    MLAdiff(N+1,N+1) = DA/L/dx+DAE/Lout;
    MLAdiff(N+1,N) = -DA/L/dx;
    bAdiff(N+1) = DAE/Lout*AE;
    
    MLKdiff(N+1,N+1) = DK/L/dx+DKE/Lout;
    MLKdiff(N+1,N) = -DK/L/dx;
    bKdiff(N+1) = DKE/Lout*KE;
    

    b1 = MRGdiff*G - alpha*(3/2*G.*KH^n./(KH^n+G.^n)-1/2*Gm.*KH^n./(KH^n+Gm.^n))*dt-deltaGrowth*(3/2*G.*A-1/2*Gm.*Am)*dt-deltaGDecay*(3/2*G-1/2*Gm)*dt;
    bGdiff(2:N) = b1(2:N);
    Gm = G;
    G = MLGdiff\bGdiff;
    G = G.*(G>0);
    
    b1 = MRAdiff*A + alpha*(3/2*G.*KH^n./(KH^n+G.^n)-1/2*Gm.*KH^n./(KH^n+Gm.^n))*dt -deltaGrowth*(3/2*G.*A-1/2*Gm.*Am)*dt;
    bAdiff(2:N) = b1(2:N);
    Am = A;
    A = MLAdiff\bAdiff;
    A = A.*(A>0);
    
    VK = VK0 + deltaK*K;
    VL = VL0 + deltaL*K;
    V = V + (-gK*nK.^4.*(V-VK)-gL*(V-VL))*dt;
    nK = nK + (alpha0*S.^m./(Sth^m+S.^m).*(1-nK)-beta*nK)*dt;
    S = S +(alphas*(Vth-V)./(exp((Vth-V)/sigma)-1)-gammaS*S)*dt;
    
    
    b1 = MRKdiff*K + (F*gK*nK.^4.*(V-VK)-gammaE*K)*dt;
    bKdiff(2:N) = b1(2:N);
    Km = K;
    K = MLKdiff\bKdiff;
    K = K.*(K>0);
    
    
    U = sum((betar*(3/2*A.*G-1/2*Am.*Gm))*dt)/(N+1);
    L = L+U;
    dL(iters) = sum(V);
    
    %dlong(iters,:) = G.*A;
  %  plot(x,K);
  %  title('Potassium');
end

%Rescale Up
A = A*Abar;
G = G*Gbar;
K = K*Lbar;
x = x*Lbar;
L = L*Lbar;
dL = dL*Lbar;

 
 iter = linspace(1,maxiters,maxiters)*dt;
% figure
% plot(x,G,x,A,'LineWidth',2);
% legend('Glutamate','Amonium')

% figure
% plot(iter(maxiters/2:maxiters),dL(maxiters/2:maxiters));
% title('Voltage');

figure
plot(iter,dL);
title('Mean Voltage');

% figure
% plot(x,A);
% title('Ammonium');

figure
plot(x,K);
title('Potassium');

figure
plot(x,V);
title('Voltage');

figure
plot(x,S);
title('Stress');


figure
plot(x,(alphas*(Vth-V)./(exp((Vth-V)/sigma)-1)-gammaS*S));
title('Change in Stress');







%mesh(dlong)