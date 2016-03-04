%% Matlab code to evaluate a single cell time evolutions

%% Default command to clear the Workspace and close all the figures
close all
clear all
clc

%% GENES PARAMETERS
ns  = 4;                      % Number of phosphorylation sites for Swi6
nw  = 12;                     % Total number of phosphorylation sites for Whi5
nw1 = 4;                      % Number of Whi5 phosphorylation sites, useful for G1/S regulon activation
nw2 = nw - nw1;               % Number of Whi5 phosphorylation sites, useless for G1/S regulon activation
nd  = ns + 2;                 % Length of dimer-vector Pd
nt  = nd*(nw1 + 1)*(nw2 + 2); % Length of trimer-vector Pt
NA  = 136;                    % Number of A genes
NB  = 63;                     % Number of B genes
NC  = 36;                     % Number of C genes

[Dd,Dk,Fd,Ud,Rt,Rb,Rk,Rw,Fr,Ur,Hr,nr]  = generate_matrices(ns,nw1,nw2,'wt');
n   = [nd, nr, nw1, nt];
N   = [NA,NB,NC];

[T,X,T1,T2,Ps,AG,MA,Vs,alpha] = d1_one_cell(n, N, Dd, Dk, Rb, Rt, Rk, Rw, Fd, Fr, Ud, Ur, Hr);

TG1 = T1 + T2;
[K_hill,n_hill] = Best_Hill(T,AG);

disp(['T1 duration value = ',num2str(T1),' min'])
disp(['T2 duration value = ',num2str(T2),' min'])
disp(['TG1 duration value = ',num2str(TG1),' min'])
disp(['Hill coefficient value = ',num2str(n_hill)])
disp(['Median point value = ',num2str(K_hill)])
disp(['Ps value = ',num2str(Ps/1e10),'1e10'])
disp(['Vs value = ',num2str(Vs)])
disp(['Linear growth rate (alpha) value = ',num2str(alpha)])
%% Relevant plots
Figures(T, X, AG, MA, Fd, Fr, N, n)