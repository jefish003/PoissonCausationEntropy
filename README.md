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


Now we want Poisson oCSE, this can be achieved by running test_Poisson_oCSE.m, the code of which is provided below:

```
clear
clc
close all



Date = '08-Apr-2021';
nsamp = 1200;
ERparam = 0.05;

loaddir = [pwd '/Data/'];
loadname = ['Data_ER_' num2str(ERparam) '_nsamp_' num2str(nsamp) '_' Date];
load([loaddir loadname '.mat'])
%==========================================================================
%================ Now Estimate the Network Using Poisson ==================
nodes = n;
B = zeros(nodes); %Initialize the adjacency matrix of the estimated network
vars = [1:nodes];
for i = 1:nodes
    i
    NODE = i;
    Y = A(:,i); %The variable we are interested in
    X = A(:,vars(vars~=i)); %All other variables
    newVars = vars(vars~=i); %Need this for later
    %% Test CSE Poisson

    % parameters
    par.ns = 1000;
    par.alpha = 0.001;

    % CSE Poisson foward
    [S] = CSE_forward_Poisson_Est(Y,X,par);
    par.alpha = 0.001;
    [Snew] = CSE_backward_Poisson_Est(Y,X,S,par);
    B(i,newVars(Snew)) = 1;
end

csvwrite([loaddir 'PoissonoCSE_Results_' loadname '.csv'],B)
```

Of course you will want to change the date for when your file(s) was (were) run.
