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

Next we will run glasso, note again that the code for glasso is written in R, and so it is contained in the file Glasso_ForPoisson_WithBIC.R. In this file we boxcox transform the data first and then run Glasso and choose lambda based upon the BIC condition. 

The R code for this is given below:

```
library(glasso)
library(Matrix)
library(MASS)
library(EnvStats)
library(rdetools)



ReadFolder = "/Data/"
RData = "Data_ER_0.05_nsamp_1200_"
ReadDate = "08-Apr-2021.csv"

WD = dirname(rstudioapi::getSourceEditorContext()$path)
ReadName = paste(ReadFolder,RData,ReadDate,sep="")

WriteDate = Sys.Date()
WriteFront = "Glasso_Results_"
WriteName = paste(ReadFolder,WriteFront,RData,WriteDate,".csv",sep="")

CSV = read.csv(paste(WD,"/",ReadName,sep=""),head=FALSE)

Mat = as.matrix(CSV,nrows=nrow(CSV),ncols=50)
BC = boxcox(CSV[,1]+1,objective.name="Log-Likelihood", optimize=TRUE)
lambda = BC[1]
lambda = as.double(lambda)
M = as.matrix(CSV[,1])
M = M+1
M1 = (M^lambda)/lambda
LS = logspace(-2,0,1000)

for(i in 2:50){
  print(i)
  BC = boxcox(CSV[,i]+1,objective.name="Log-Likelihood", optimize=TRUE)
  lambda = BC[1]
  lambda = as.double(lambda)
  M = as.matrix(CSV[,i])
  M = M+1
  M1 = cbind(M1,(M^lambda)/lambda)
}

S = var(M1)

BIC = -1e100
for(j in 1:1000){
  alpha = LS[j]
  GL = glasso(S,alpha,nobs=nrow(CSV))
  LogLikeliNum = as.double(GL[3])
  Adj = matrix(unlist(GL[2]),nrow=50,byrow=TRUE)
  Wh = which(Adj!=0)
  Adj[Wh] = 1
  B2 = length(Wh)*log(nrow(CSV))/2
  print(paste(i,BIC))
  #BIC Maximization Here
  if(LogLikeliNum-B2>BIC){
    A = Adj
    BIC = LogLikeliNum-B2
    print(BIC)
  }
}

write.table(A,paste(WD,"/",WriteName,sep=""),sep=",",row.names=FALSE,col.names=FALSE)
```
Finally we are ready for the hybrid method. The hybrid method first applies Glasso in the place of the forward Poisson oCSE, and then applies the backward Poisson to weed out extraneous edges that were discovered. In the example case given with a sample size of 1200 and a network which only has an expected average degree of 2.5 the combo forward Poisson and backward Poisson are expected to find all of the true edges, this is not the case for smaller sample size and in the paper we show that the hybrid method performs better in this small data scenario, note of course that the engine for weeding out spurious edges is still oCSE though! The code for running the hybrid method is labeled test_Hybrid_GLASSO_PoissonoCSE.m and a code snippet is given below:

```
clear
clc
close all


Date = '08-Apr-2021';
Date2 = '2021-04-08';
nsamp = 1200;
ERparam = 0.05;

loaddir = [pwd '/Data/'];
loadname = ['Data_ER_' num2str(ERparam) '_nsamp_' num2str(nsamp) '_'];
load([loaddir loadname Date '.mat'])

%GLASSO inferred network
C = csvread([loaddir 'Glasso_Results_' loadname Date2 '.csv']);
C = C - eye(n);

nodes = n;
D = zeros(nodes); %Initialize the adjacency matrix of the estimated network
vars = [1:nodes];
for i = 1:nodes
    i
    NODE = i;
    Y = A(:,i); %The variable we are interested in
    X = A(:,vars(vars~=i)); %All other variables
    newVars = vars(vars~=i); %Need this for later
    %% Test CSE Poisson
    V = find(C(i,:));
    S = [];
    for l = 1:length(V)
        S = [S find(newVars==V(l))];
    end
    % parameters
    par.ns = 1000;
    par.alpha = 0.001;
   
    [Snew] = CSE_backward_Poisson_Est(Y,X,S,par);
    D(i,newVars(Snew)) = 1;
end

csvwrite([loaddir 'Hybrid_GLASSO_PoissonoCSE_Results_' loadname Date '.csv'],D)
```

Finally we examine the performance of the methods. As mentioned above in this case using either the Poisson oCSE method or the hybrid method should give us nearly identical results. However when the sample size is small (or the average degree is large) the hybrid method with oCSE as the spurious edge removal performs best. The code for calculating true and false positive rates is called PerfEval.m and the code snippet is given below:

```
clear
clc
close all


%Performance evaluation

%Load in the true matrix (Aa)
load([pwd '/Data/' 'Data_ER_0.05_nsamp_1200_08-Apr-2021.mat'])

%Load in the matrix inferred by the Poisson oCSE and calculate the FPR and
%TPR
PoissInferred = csvread([pwd '/Data/' 'PoissonoCSE_Results_Data_ER_0.05_nsamp_1200_08-Apr-2021.csv']);
PoissFP = sum(sum((Aa-PoissInferred)<0));
PoissFPR = PoissFP/sum(sum(Aa))
PoissTPR = (sum(sum(PoissInferred)) - PoissFP)/sum(sum(Aa))

%Load in the matrix inferred by Glasso
GlassoInferred = csvread([pwd '/Data/' 'Glasso_Results_Data_ER_0.05_nsamp_1200_2021-04-08.csv']);

%Glasso always returns diagonal entries by design
GlassoInferred = GlassoInferred - eye(n);
GlassoFP = sum(sum((Aa-GlassoInferred)<0));
GlassoFPR = GlassoFP/sum(sum(Aa))
GlassoTPR = (sum(sum(GlassoInferred)) - GlassoFP)/sum(sum(Aa))

%Hybrid method
HybridInferred = csvread([pwd '/Data/' 'Hybrid_GLASSO_PoissonoCSE_Results_Data_ER_0.05_nsamp_1200_08-Apr-2021.csv']);
HybridFP = sum(sum((Aa-HybridInferred)<0));
HybridFPR = HybridFP/sum(sum(Aa))
HybridTPR = (sum(sum(HybridInferred)) - HybridFP)/sum(sum(Aa))
```
