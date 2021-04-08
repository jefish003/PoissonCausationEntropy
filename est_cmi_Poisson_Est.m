function [I_est] = est_cmi_Poisson_Est(X,Y,Z,S);
% (parametric) estimation of conditional mutual information, assuming Poisson distribution

% Input:
%  X: n-by-dx, n points in the dx-dimensional sample space
%  Y: n-by-dy, n points in the dy-dimensional sample space
%  Z: n-by-dz, n points in the dz-dimensional sample space
%  S: d-by-d covariance matrix of the joint variable (X,Y,Z), d = dx + dy + dz

% Output:
%  I(X;Y|Z) using the Gaussian formula of mutual information
%  [Ref: T. M. Cover and J. A. Thomas, Elements of Information Theory, 
%        2nd ed., John Wiley & Sons, Hoboken, NJ, 2006.]

% Last update/check: Jan 29th, 2019
% Reference: 
%   J. Sun, D. Taylor, E. M. Bollt, "Causal network inference by optimal causation entropy"
%   SIAM Journal on Applied Dynamical Systems vol.14, 73-106 (2015).


%% Main

% dimensions
dx = size(X,2); 
dy = size(Y,2);
dz = size(Z,2);

% subspace covariance matrices
indX = [1:dx]; indY = [dx+1:dx+dy]; indZ = [dx+dy+1:dx+dy+dz];


% computation of cmi
if dz==0 % no conditioning
    S_est = S([indX,indY],[indX,indY]);
    %Dcov = diag(S_est);
    l_est = S_est - diag(diag(S_est));
    S_est(1:length(S_est)+1:end) = S_est(1:length(S_est)+1:end)-sum(l_est,2)';
    Dcov = diag(S_est) + sum(l_est,2);
    I_est = sum(PoissEnt(Dcov))-PoissonJointEntropy_New(S_est);
else % with conditioning
    Sa = S - diag(diag(S));
    Sa = sum(Sa,2);
    SS = S;
    SS(1:length(SS)+1:end) = SS(1:length(SS)+1:end)- Sa';
    SS(indX,indX) = SS(indX,indX)+ S(indX,indY); 
    SS(indY,indY) = SS(indY,indY) + S(indY,indX);
    
    %For Finding H(Y,Z)
    S_est1 = SS([indY,indZ],[indY,indZ]);
    %For finding H(X,Z)
    S_est2 = SS([indX,indZ],[indX,indZ]);
    HYZ = PoissonJointEntropy_New(S_est1);
    %For finding H(Z) noting that Z may be a set of variables in which case
    %we may need the joint entropy of the variable Z.
    SindZ = SS(indZ,indZ);
    HZ = PoissonJointEntropy_New(SindZ);
    HXYZ = PoissonJointEntropy_New(S-diag(Sa));
    HXZ = PoissonJointEntropy_New(S_est2);
    H_YZ = HYZ - HZ;
    H_XYZ = HXYZ - HXZ;
    I_est = H_YZ - H_XYZ;
    
end

