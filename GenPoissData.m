clear
clc
close all

%Generate Poisson data as csv file

%Number of nodes
n = 50; 
Zp = n*(n-1)/2;

lambdas = ones(n+Zp,1);
%Generate noise
noiseLam = ones(n,1)*0.5;

nsamples = 1200;
%Probability for ER network;
ERparam = 0.05; 

Aa = ER_Adj(n,ERparam);

A = PoissonDataMaker(Aa,nsamples,lambdas,noiseLam);

csvwrite(['Data/Data_ER_' num2str(ERparam) '_nsamp_' num2str(nsamples) '_' date '.csv'],A)
save(['Data/Data_ER_' num2str(ERparam) '_nsamp_' num2str(nsamples) '_' date '.mat'])