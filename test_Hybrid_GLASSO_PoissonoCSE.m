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