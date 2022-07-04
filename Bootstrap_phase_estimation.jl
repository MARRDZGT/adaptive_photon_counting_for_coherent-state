using MATLAB
#using Statistics
using DataFrames
using Plots
#using Distributions
using CSV

#---------------------
#Parameters
#---------------------

Adapt_Steps = 80.0
Boot_Reps = 10000.0

a=[4.0;200.0]
#a = [15.0:5.0:50.0;]
#b = [20.0:20.0:100.0;]
#c = [200.0:100.0:500.0;]
#c = []
#d = [900]
#a=[1.0]

#MPNS = [a;c]
MPNS = a
#print(a)
print(MPNS)

#---------------------
#Algorithm
#---------------------

function Adap_Algot(Ad_Steps, N_Rep, MPN)

    mat"""
        close all;clear all;clc
    """

    @mput Ad_Steps
    @mput N_Rep
    @mput MPN

    mat"""
        %addpath('/home/marco/Documentos/Phd/Paper/Quantum_Optimal_Design/Julia-Script')
        Nphi = 1024;
        PhiStep = 2*pi/Nphi;
        phi = (0: PhiStep : 2*pi - PhiStep);
        PNR = 3;

        %display([MPN,Ad_Steps,MPN])
        
        
                %M = 20; %number of repetitions
        M = N_Rep; %number of repetitions
        %L = 10; %number of adaptive steps
        L = Ad_Steps; %number of adaptive steps
        %Atot = (100); %total MPNs to run
        Atot = (MPN); %total MPNs to run

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
        
        
                Res = zeros(M,2);

        for kk = 1:1:length(Atot)

            a2 = Atot(kk)./L; %MPN per adaptive step
            parfor m = 1:M
              Res(m,:) = HEL(kk,m,Nphi,PhiStep,PNR,L,Atot, DE, Vis, DC, K, Kv,Kv2,pr,xx,PhiEstS, PhiEstMI,options,Pn,mpnMI,PnMI,ap,a2,phi,LAST,LAST1);
            end
        %VHS(kk) = VarHS;
        %VHMI(kk) = VarHMI;

        clear a2
        end
        toc

    """
    @mget Res
    return(Res)
end

#-----------------
#Obtain Data
#-----------------

All = [Adap_Algot(Adapt_Steps, Boot_Reps, MPNS[k]) for k in 1:length(MPNS)]

#------------------
#Write Data
#------------------

Data_MI = Array{Float64}(undef, floor(Int64, Boot_Reps),length(MPNS))
for j = 1:length(MPNS)
    Data_MI[:,j] = All[j][:,1]
end

Data_SH = Array{Float64}(undef, floor(Int64, Boot_Reps),length(MPNS))
for j = 1:length(MPNS)
    Data_SH[:,j] = All[j][:,2]
end

D_MI= DataFrame(Data_MI)
D_SH= DataFrame(Data_SH)

col_names = MPNS

rename!(D_MI, Symbol.(col_names))
rename!(D_SH, Symbol.(col_names))

CSV.write("D_MI_PNR3_MPN-L80-3.csv",  D_MI, writeheader=true)
CSV.write("D_SH_PNR3_MPN-L80-3.csv",  D_SH, writeheader=true)

#---------------------
#Histograms
#---------------------

#MI_1 = histogram(Data_MI[:,1], bins=10)
#savefig("h_MI_1.png")

#MI_F = histogram(Data_MI[:,length(MPNS)], bins=10)
#savefig("h_MI_F.png")

#SH_1 = histogram(Data_SH[:,1], bins=10)
#savefig("h_SH_1.png")

#SH_F = histogram(Data_SH[:,length(MPNS)], bins=10)
#savefig("h_SH_F.png")

#---------------------
#Holevo Variance
#---------------------

#phi = pi

#function Holevo_Variance(Dat)
#    sh = [cos(Dat[i]-phi) for i in 1:length(Dat)]
#    HV =mean(sh)^(-2) - 1
#    return(HV)
#end


