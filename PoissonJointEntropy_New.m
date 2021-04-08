function [I_est] = PoissonJointEntropy_New(Cov);

%For Calculation of Poisson Joint Entropy
%
%Inputs: Cov - The (n x n) covariance (or correlation) matrix, where n is
%             the number of random variables.
%
%Outputs: I_est - The estimated Joint Entropy between the random variables.
%
%==========================================================================
%==========================================================================
%========================= Written By: Jeremie Fish =======================
%====================== Last Edited on January 29th 2019 ==================
%==========================================================================
%==========================================================================

%Find lambda_12, lambda_13,... lambda_1n, lambda_23, lambda_24,...
T = triu(Cov,1);
T= T(:);

%Now find lambda_11,lambda_22,...,lambda_nn
H = diag(Cov);

%The Joint Entropy is approximately the sum of the single variate entropies
%added to the lambdas we found in the T term. So first find the single 
%variate entropies summed together which is done in Ent1 below.
Ent1 = sum(PoissEnt(H(:)));

%Now sum all of the lambdas contained in T
Ent2 = sum(T(:));

%Now calculate the approximate entropy.
I_est = Ent1+Ent2;