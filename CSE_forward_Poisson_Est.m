function [S Covar] = CSE_forward_Poisson_Est(Y,X,par)


% Input:
%  X: T-by-nx time series
%  Y: T-by-ny time series
%  par: parameters used in the method
%    par.ns: number of shuffles in the permutation test
%    par.alpha: alpha level
%    par.opt: option of Covariance or Correlation

% Output:
%  S: set of indices of the columns of X that "cause" Y

% Last update/check: Jan 29th , 2019
% Reference: 


%% Main

% initialize
S = []; Z = [];

% dimensions
dy = size(Y,2); dx = size(X,2);
indY = [dx+1:dx+dy]; indZ = [];

% covariance matrix
try
if par.opt == 'Covariance';
    Sigmatot = cov([X,Y]);%corr([X,Y]);
end
if par.opt == 'Correlation'
    Sigmatot = corr([X,Y]);
end
catch
    Sigmatot = cov([X,Y])./diag(cov([X,Y]));
end

Sigmatot(Sigmatot<0) = 0;
Covar = cov([X,Y]);

% misc
T = size(X,1); n = size(X,2);
burned = zeros(1,n); % record whether or not a node has been included in S

% main loop
isstop = 0;
i = 0;
saver = zeros(1,n);
iter = 0;
%MY STUFF!
Iresults = [];
Tresults = [];
while isstop==0
    iter = iter+1;
    CMI = -Inf(1,n);
    for j = 1 : n
        if burned(j)==0
            indX = j;
            
            ind = [indX,indY,indZ];
            Sigma = Sigmatot;
           
            aba= est_cmi_Poisson_Est(X(:,j),Y,Z,Sigma(ind,ind));
            CMI(j) = aba;
        end
    end
    saver(iter,:) = CMI;
    [val,pos] = max(CMI);
    if val>-Inf 
        i = i+1;
        data.I = val;
        data.X = X(:,pos); data.Y = Y; data.Z = Z; data.ind = ind;data.S = Sigmatot;
        [result] = shuffletest_cmi_Poisson_Est(data,par);
        if result.ans>0
            disp(['Success'])
           
            S = [S,pos];
            Z = [Z,X(:,pos)]; indZ = [indZ,pos]; 
            burned(pos) = 1;
        else
            isstop = 1;
        end
    else
        isstop = 1;
    end
end

