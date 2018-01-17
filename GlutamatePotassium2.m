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

Lbar = .1;
tbar = GE/(alpha+deltaGrowth+deltaGDecay);
Gbar = GE;
Abar = 1E-4;
Kbar = 1E-4*0+1;

gK = 30*60; % 1/hour
gL = .2*60;% 1/hour
VK0 = -380;%-380; %mV
VL0 = -156; %mV
Sth = .8;%40/1000; %mMole
Vth = -150;%mMole
Gth = 5;
alpha0 = 2*60; %1/hour
beta = 1.3*60; %1/hour
m = 10;
F = 5.6; %mMole/mV
sigma = .2;%mV
deltaK = 1;%mV/mM
deltaL = 8;%mV/mM
gammaS = .1*60;%.1*60;%1/hour
gammaE = 10*60;%1/hour
gammaT = 4*60;%1/hour
alphas = 1*60/1000;%mmole/(hour mV)
alphat = 1*60/1000;%mmole/(hour mV)
DK = 13.8E-6*10^2*60*60/100;% mm^2/hour
DKE = 13.8E-6*10^2*60*60/100;% mm^2/hour
KE = 1/Kbar/1;


%Scale terms
D = D/Lbar^2;
DE = DE/Lbar^2;
DA = DA/Lbar^2;
DAE = DAE/Lbar^2; 
DK = DK/Lbar^2;
DKE = DKE/Lbar^2;
GE = 1;
Gth = Gth/Gbar;

L = .1/Lbar; %Depth of Biofilm in mm %
Lout = .01/Lbar;%Diffusion layer


N = 100;
dx = 1/N;
x = linspace(0,L,N+1);
T = .5/tbar;%Terminal Time
dt = .00001/tbar;%Time Step
maxiters = round(T/dt);

G = ones(N+1,1)*1/Gbar;
A = ones(N+1,1)*0.004/Abar;
K = ones(N+1,1).*heaviside(-x'+L/20)/Kbar/1*100*0+ones(N+1,1)*KE/Kbar;
Kin = ones(N+1,1)/Kbar*100;
V = ones(N+1,1)*(Vth*0-156);
nK = ones(N+1,1)*.5;
S = ones(N+1,1)*0;
VK = ones(N+1,1)*VK0;
VL = ones(N+1,1)*VL0;
Gm = ones(N+1,1)*1/Gbar;
Am = ones(N+1,1)*0.004/Abar;
Km = ones(N+1,1).*heaviside(-x'+L/20)/Kbar/1*100*0+ones(N+1,1)*KE/Kbar;
growth = zeros(N+1,1);


U = zeros(N+1,1);
dL = zeros(maxiters,1);
MLGdiff = zeros(N+1); %Matrices for Glutamate
MRGdiff = zeros(N+1);
bGdiff = zeros(N+1,1);
% MLAdiff = zeros(N+1); %Matrices for Ammonium
% MRAdiff = zeros(N+1);
% bAdiff = zeros(N+1,1);
MLKdiff = zeros(N+1); %Matrices for Potassium
MRKdiff = zeros(N+1);
bKdiff = zeros(N+1,1);
Sswitch = zeros(N+1,1);

%Set up Matrices
MLGdiff(1,1) = 1/dx;
MLGdiff(1,2) = -1/dx;
% MLAdiff(1,1) = 1/dx;
% MLAdiff(1,2) = -1/dx;
MLKdiff(1,1) = 1/dx;
MLKdiff(1,2) = -1/dx;

Smean = zeros(N+1,1);
for iters = 1:maxiters
    %Set Up Matrices, which depend on L, which changes!!!
    for i = 2:N
        MLGdiff(i,i) = 1 + 2*D*dt/dx^2/2/(L^2);
        MLGdiff(i,i-1) = -D*dt/dx^2/2/(L^2);
        MLGdiff(i,i+1) = -D*dt/dx^2/2/(L^2);
    
        MRGdiff(i,i) = 1 + -2*D*dt/dx^2/2/(L^2);
        MRGdiff(i,i-1) = D*dt/dx^2/2/(L^2);
        MRGdiff(i,i+1) = D*dt/dx^2/2/(L^2);
    
%         MLAdiff(i,i) = 1 + 2*DA*dt/dx^2/2/(L^2);
%         MLAdiff(i,i-1) = -DA*dt/dx^2/2/(L^2);
%         MLAdiff(i,i+1) = -DA*dt/dx^2/2/(L^2);
%     
%         MRAdiff(i,i) = 1 + -2*DK*dt/dx^2/2/(L^2);
%         MRAdiff(i,i-1) = DK*dt/dx^2/2/(L^2);
%         MRAdiff(i,i+1) = DK*dt/dx^2/2/(L^2);
        
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

%     MLAdiff(N+1,N+1) = DA/L/dx+DAE/Lout;
%     MLAdiff(N+1,N) = -DA/L/dx;
%     bAdiff(N+1) = DAE/Lout*AE;
%     
    MLKdiff(N+1,N+1) = DK/L/dx+DKE/Lout;
    MLKdiff(N+1,N) = -DK/L/dx;
    bKdiff(N+1) = DKE/Lout*KE;
    

    b1 = MRGdiff*G -deltaGrowth*(max(V/Vth,0)).*(3/2*(G)./(G+Gth)*2/3-1/2*Gm./(Gm+Gth))*dt-deltaGDecay*(max(0,V/Vth)).*(3/2*(G)*2/3-1/2*Gm)*dt;
    bGdiff(2:N) = b1(2:N);
    Gm = G;
    G = MLGdiff\bGdiff;
    G = G.*(G>0);
    
%     b1 = MRAdiff*A + alpha*(3/2*G.*KH^n./(KH^n+G.^n)-1/2*Gm.*KH^n./(KH^n+Gm.^n))*dt -2*deltaGrowth*(V/Vth).*(3/2*G.*A-1/2*Gm.*Am)*dt;
%     bAdiff(2:N) = b1(2:N);
%     Am = A;
%     A = MLAdiff\bAdiff;
%     A = A.*(A>0);
    
    VK = VK0 + deltaK*K;
    VL = VL0 + deltaL*K;
    V = V + (-gK.*nK.^4.*(V-VK)-gL*(V-VL))*dt;
    %(-gK*nK.^4.*heaviside(V-VK).*(V-VK)-gL*heaviside(V-VL).*(V-VL))*dt+ (-gK*nK.^4.*heaviside(VK-V).*(V-VK)-gL*heaviside(VL-V).*(V-VL))*dt;
    S = (S + (log(exp((1-max(V/Vth,0).*G)+1)).*((V+380).^10./((V+380).^10+200^10))*1000 - gammaS.*(S)*10)*dt);
%     S = (S + (log(exp((1-max(V/Vth,0).*G)+1)).*(1-Smean) - gammaS.*(S)/10)*dt);
%     Smean = (Smean*1000+Sswitch)/(1000+1);
%     Sswitch = Sswitch + ((1-Sswitch).*heaviside(S-.5) - Sswitch.*(heaviside(.01-S)))*100*dt;
    S = S.*(S>0);
    nK = nK + (alpha0*(S).^2./(1^2+(S).^2).*(1-nK)+.01*(10-V)./(exp((10-V)/100)-1).*(1-nK)-.125*(exp(-V/80)).*nK)*dt; 
    %nK = nK + (alpha0*(S.*Sswitch).^2./(.1^2+(S.*Sswitch).^2).*(1-nK)+.01*(10-V)./(exp((10-V)/100)-1).*(1-nK)-.125*(exp(-V/80)).*nK)*dt; 

    
    b1 = MRKdiff*K + (F.*gK.*nK.^4.*(V-VK)+F*gL*(V-VL).*K-gammaE*K*0)*dt;
    %+ (-gK*nK.^4.*heaviside(V-VK).*(V-VK)-gL*heaviside(V-VL).*(V-VL))*dt+(-gK*nK.^4.*heaviside(VK-V).*(V-VK)-gL*heaviside(VL-V).*(V-VL))*dt;
    bKdiff(2:N) = b1(2:N);
    Km = K;
    K = MLKdiff\bKdiff;
    K = K.*(K>0);
   % Kin = Kin - (F.*gK.*nK.^4.*(V-VK)+gL*(V-VL)*F.*K-gammaE*K)*dt;
    %- (-gK*nK.^4.*heaviside(V-VK).*(V-VK)-gL*heaviside(V-VL).*(V-VL))*dt - (-gK*nK.^4.*heaviside(VK-V).*(V-VK)-gL*heaviside(VL-V).*(V-VL))*dt;
   % Kin = Kin.*(Kin>0);
    
    
  %   U = sum((betar*(3/2*A.*G-1/2*Am.*Gm))*dt)/(N+1)/1000;
   % L = L+U;
[value, place] = min(V);
   dL(iters) = mean(place/N*L*Lbar);
%   dL(iters) = mean(V);
    
    %dlong(iters,:) = G.*A;
  %  plot(x,K);
  %  title('Potassium');
end

%Rescale Up
% A = A*Abar;
% G = G*Gbar;
% K = K*Kbar;
% Kin = Kin*Kbar;
% x = x*Lbar;
% L = L*Lbar;
%dL = dL*Lbar;

 
 iter = linspace(1,maxiters,maxiters)*dt;
% figure
% plot(x,G,x,A,'LineWidth',2);
% legend('Glutamate','Amonium')

% figure
% plot(iter(maxiters/2:maxiters),dL(maxiters/2:maxiters));
% title('Voltage');

figure
plot(iter(maxiters/2*0+1:maxiters),(dL(maxiters/2*0+1:maxiters)));
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
plot(x,G);
title('Glutamate');


figure
plot(x,(log(exp((1-max(V/Vth,0).*G)+1)).*(1-Sswitch) - gammaS.*(S)));
%plot(x,(1.-(max(V/Vth,0).*G)./((max(V/Vth,0).*G)+.00001)-0*gammaS*S));
title('Change in Stress');

figure
plot(x,((1-Sswitch).*(heaviside(S-1)+heaviside(S-Smean-.5)) - Sswitch.*(1-heaviside(S-.001))));






%mesh(dlong)