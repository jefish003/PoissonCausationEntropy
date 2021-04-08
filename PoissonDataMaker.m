function [X] = PoissonDataMaker(A,n,lambdas,noiseLam);

%Inputs: A- The Adjacency Matrix of the system must be (p x p)
%        n- The number of observations (natural number)
%        lambdas - The rates (p+p(p-1)/2 x 1)
%        noiseLam - The rates of the noise (p x 1)
%
%Output: X- The Poisson distributed data
%
%==========================================================================
%==========================================================================
%====================== Written By Jeremie Fish ===========================
%======================= Last Edited 12-05-2018 ===========================
%==========================================================================
%==========================================================================

p = size(A,1);
Zp = p*(p-1)/2;
%Note A is assumed to be a simple graph (no self edges). Thus we force it
%to be a simple graph below.
A(1:p+1:end) = 0;

%Get the permutation matrix P as discussed in the paper: 'A Local Poisson 
%Graphical Model for Inferring Networks From Sequencing Data' by Genevera I
%Allen and Zhandong Liu.
P = PermutationMatrix(p);

%Find 1_p as discussed in the paper
One_p = ones(p,1);

%Now find tri(A) as discussed in the paper.
NK = nchoosek(1:p,2);
Inds = sub2ind([p p],NK(:,1),NK(:,2));
T = A(Inds);

%I_p
I_p = eye(p);

%Matrix that P is taken the Hadamard Product with
T_p = One_p*T';

%The matrix B from the above listed paper
B = [I_p P.*T_p]';

%Now get Y from the paper
Lambs = lambdas'.*ones(n,Zp+p);
Y = poissrnd(Lambs);

%Now for E
Lambs = noiseLam'.*ones(n,p);
E = poissrnd(Lambs);

%Finally we can calculate X
X = Y*B + E;
