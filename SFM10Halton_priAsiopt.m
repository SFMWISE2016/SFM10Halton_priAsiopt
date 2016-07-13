% clearing work&preparing
clear
close all
clc
% parameter setting
S0 =50;       % Price of underlying today
X = 50;       % Strike at expiry
mu = 0.04;    % expected return
sig = 0.1;    % expected vol.
r = 0.03;     % Risk free rate
dt = 1/365;   % time steps
steps = 50;   % days to expiry
T = dt*steps; % years to expiry
nsims = 1000; % Number of simulated paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function hs = HaltonSequence(n,b)    
hs = zeros(n,1);
for idx = 1:n
    hs(idx) = localHaltonSingleNumber(idx,b);
end
function hn = localHaltonSingleNumber(n,b)
n0 = n;
hn = 0;
f = 1/b;
while (n0>0)
    n1 = floor(n0/b);
    r = n0-n1*b;
    hn = hn + f*r;
    f = f/b;
    n0 = n1;
end
%%%%%%%%%%% Function to generates the first n numbers in Halton's low discrepancy sequence with base b

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sequence = BoxMueller(hs1,hs2)
R = sqrt(-2*log(hs1));
Theta = 2*pi*hs2;
P = R.*sin(Theta);
Q = R.*cos(Theta);
sequence = [P;Q];
%%%%%%%%%%Function to convert uniform variates into normal variates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = AssetPathsHalton(S0,mu,sig,dt,steps,nsims) 
nu = mu - sig*sig/2;
pVec = primes(1e5); % This will give 9592 numbers
bases = pVec(1:2*steps);
epsilon = nan(steps,nsims);
for idx = 1:steps
    epsH1 = HaltonSequence(nsims/2,bases(idx));
    epsH2 = HaltonSequence(nsims/2,bases(steps+idx));
    epsilon(idx,:) = BoxMueller(epsH1,epsH2)';
end
% Generate potential paths
S = S0*[ones(1,nsims); ...
   cumprod(exp(nu*dt+sig*sqrt(dt)*epsilon))];
%%%%%%%%%%%% Function to generate sample paths for assets assuming geometric 
%%%%%%%%%%%% Brownian motion and using Halton's quasi-random sequence

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = AssetPathsHalton(S0,mu,sig,dt,steps,nsims);
time = steps:-1:0;
plot(time,S,'Linewidth',2);
set(gca,'XDir','Reverse','FontWeight','bold','Fontsize',24);
xlabel('Days to Expiry','FontWeight','bold','Fontsize',24);
ylabel('Asset Price','FontWeight','bold','Fontsize',24);
title('Simulated Asset Paths','FontWeight','bold','Fontsize',24);
grid on
set(gcf,'Color','w');
%%%%%%%%%%%% Plot the asset paths
PutPayoffT = max(X-mean(S),0);
CallPayoffT = max(mean(S)-X,0);
putPrice = mean(PutPayoffT)*exp(-r*T)
callPrice = mean(CallPayoffT)*exp(-r*T)
%%%%%%%%%%%%Price a standard Asian Put and Call option





