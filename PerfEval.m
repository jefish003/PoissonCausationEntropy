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