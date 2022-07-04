function [] = HEL(kk,m,Nphi,PhiStep,PNR,L,Atot, DE, Vis, DC, K, Kv,Kv2,pr,xx,PhiEstS, PhiEstMI,options,Pn,mpnMI,PnMI,ap,a2,phi,LAST,LAST1)
global VarHS;
global VarHMI,
global Est_MI_A,
global Est_SH_A,
global PhiEstMI_A,
global PhiEstS_A,

%Set prior to uniform ditribution
PpriorS = ones([1, Nphi])./Nphi;
PpriorMI = ones([1, Nphi])./Nphi;

for k = 1:1:L
    %     PhAdjS = zeros(kk, m, k);
    %     PhAdjMI = zeros(kk, m, k);
    %     PhiEstS_A = zeros(k);
    %     PhiEstMI_A = zeros(k);
    % %     VarHS = zeros(kk, k);
    % %     VarHMI = zeros(kk, k);
    %     betaMI = zeros(kk, m, k);
    %     betaS = zeros(kk, m, k);
    %     DataMS =  zeros(m, k);
    %     DataMMI=  zeros(m, k);
    %     PhiEstSmax = zeros(kk, m, k);
    %     PhiEstS = zeros(kk, m, k);
    %     PhiEstMI = zeros(kk, m, k);
    %     PhiEstMImax = zeros(kk, m, k);
    %     PpriorS = ones([1, Nphi])./Nphi;
    %     PpriorMI = ones([1, Nphi])./Nphi;
    %     PmaxIndexS = 0;
    %     PmaxIndexMI = 0;


        %[kk m k]
    %==============================================================================
    %Calculate the optimal phase adjustment and local oscillator amplitude

    %To optimize, we need to shift the distribution so that it is centered on
    %the estimate. But not for the first measurement since the prior is uniform
    if k == 1
        PpriorMS = repmat(PpriorS, PNR+1, 1);
        PpriorMMI = repmat(PpriorMI, PNR+1, 1);
    end

    %Shift prior by estimate, so it is centered at phi=0
    if k ~= 1
        PpriorShiftedS = circshift(PpriorS',-PmaxIndexS + 1)';
        PpriorMS = repmat(PpriorShiftedS, PNR+1, 1);
        PpriorShiftedMI = circshift(PpriorMI',-PmaxIndexMI + 1)';
