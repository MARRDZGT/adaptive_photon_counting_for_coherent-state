close all;clear all;clc
Nphi = 1024;
PhiStep = 2*pi/Nphi;
phi = (0: PhiStep : 2*pi - PhiStep);
PNR = 3;
M = 10000; %number of repetitions
L = 200; %number of daptive steps
Atot = (1); %total MPNs to run
DE = 1;
Vis = 1;
DC = 0;
LAST = M;
LAST1 = length(Atot);

%For computing S and MI
K = PNR + 1;
Kv = (0:1:K-1);
Kv2 = (0:1:K-2);
pr = repmat(phi, PNR+1, 1);
xx = repmat((0:1:PNR-1), Nphi, 1);

%Initializing empty estimate matricies
PhiEstS = zeros([length(Atot),M,L]);
PhiEstMI = zeros([length(Atot),M,L]);

%options for optimization
options = optimoptions('fmincon', 'algorithm', 'interior-point', 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxFunEvals', 10000, 'Display', 'Off');

Pn = @(n, mpn) (exp(-mpn).*(mpn.^n))./factorial(n);
mpnMI = @(a, b, phi, theta) DE*(a.^2 + b.^2 - 2*abs(a.*b)*Vis.*cos(phi - theta)) + DC;
PnMI = @(n, a, b, phi, theta) exp(-mpnMI(a, b, phi, theta)) .* (mpnMI(a, b, phi, theta).^n)./factorial(n);

ap = pi; %Actual phase
VHS = zeros(L,PNR);
tic

global VarHS;
VarHS = zeros(length(Atot), L);
global VarHMI;
VarHMI = zeros(length(Atot), L);

%Initializing empty estimate matricies
%global Est_MI_A;
%Est_MI_A = zeros(M,1);
%PhiEstS = zeros([length(Atot),M,L]);
%global Est_SH_A;
%Est_SH_A = zeros(M,1);
%PhiEstMI = zeros([length(Atot),M,L]);

for kk = 1:1:length(Atot)

    a2 = Atot(kk)./L; %MPN per adaptive step
    parfor m = 1:M
       HEL(kk,m,Nphi,PhiStep,PNR,L,Atot, DE, Vis, DC, K, Kv,Kv2,pr,xx,PhiEstS, PhiEstMI,options,Pn,mpnMI,PnMI,ap,a2,phi,LAST,LAST1);
    end
       %VHS(kk) = VarHS;
       %VHMI(kk) = VarHMI;

    clear a2
end
toc


