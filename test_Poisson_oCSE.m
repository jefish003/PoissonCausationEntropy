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
