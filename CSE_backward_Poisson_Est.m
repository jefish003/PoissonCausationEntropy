function [Snew] = CSE_backward_Poisson_Est(Y,X,S,par)

% Inputs:
%  X: T-by-nx time series
%  Y: T-by-ny time series
%  S: set of candidate "causal" indices
%  par: parameters used in the method
%    par.ns: number of shuffles in the permutation test
%    par.alpha: alpha level

% Output:
%  Snew: set of indices of the columns of X that "cause" y

% Last update/check: Jan 29th 2019


%% Main

% initialize
Snew = S;

% dimensions
dy = size(Y,2); dx = size(X,2);
indY = [dx+1:dx+dy]; indZ = [];
try
if par.opt == 'Covariance'
    Sigma = cov([X,Y]);%corr([X,Y]);
end
if par.opt == 'Correlation'
    Sigma = corr([X,Y]);
end
catch
    %As it turns out correlation is the best way, this is another way of
    %finding the correlation.
    Sigma = cov([X,Y])./diag(cov([X,Y]));
end
%Negative lambda values are not possible in the model, thus a calculated
%negative lambda is assumed to be 0.
Sigma(Sigma<0) = 0;

rp = randperm(length(Snew));
for i = 1 : length(S)
    indX = S(rp(i));
    if length(Snew)>1
        indZ = setdiff(Snew,indX);
        Z = X(:,indZ);
        ind = [indX,indY,indZ];
        data.I = est_cmi_Poisson_Est(X(:,indX),Y,Z,Sigma(ind,ind));
    else
        indZ = []; Z = [];
        ind = [indX,indY,indZ];
        
        data.I = est_cmi_Poisson_Est(X(:,indX),Y,Z,Sigma(ind,ind));
    end
    data.X = X(:,indX); data.Y = Y; data.Z = Z; data.ind = ind; data.S = Sigma;
    [result] = shuffletest_cmi_Poisson_Est(data,par);
    if result.ans>0
    else % remove indX
        Snew = setdiff(Snew,indX);   
    end
end


