# PoissonCausationEntropy
Code for paper related to Poisson Causation Entropy, please cite "Inference of Polygenic Factors Associated with Breast Cancer Gene Interaction Networks from Discrete Data Utilizing Poisson Multivariate Mutual Information" (Working title)

Functions are provided to perform the Poisson version of optimal causation entropy (oCSE). For more on oCSE and its Gaussian version see: 

"Causal Network Inference by Optimal Causation Entropy" by Jie Sun, Dane Taylor and Erik Bollt in SIAM Journal of Dynamical Systems.

Also a useful reference for this work:

"Causation Entropy Identifies Indirect Influences, Dominance of Neighbors and Anticipatory Couplings" by Jie Sun and Erik Bollt in Physica D.

How to use this code:

First we need to create some data, this is done in the GenPoissData.m file, the code in this file is given below for a single network run (in the paper there were 50 runs):

```
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
```


The data is now saved in the Data folder. This example file has already been run and the example data can be found in the data folder already.
Note that the data is saved both as a csv file and as a mat file. The reason for saving as csv is that we eventually need to read the data into R, as Glasso was run in R (specifically because the R version of Glasso had a prebuilt option for outputting the likelihood which is needed for the BIC calculation so rather than reinventing the wheel I used the preexisting codebase).


Now we want Poisson oCSE, this can be achieved by running test_Poisson_oCSE.m
