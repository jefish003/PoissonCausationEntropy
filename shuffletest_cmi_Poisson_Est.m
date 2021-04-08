function [result] = shuffletest_cmi_Poisson_Est(data,par);
% Shuffle test of the cmi: I(X;Y|Z)
%  H0 (null): I(X;Y|Z)=0
%  H1 (alternative): I(X;Y|Z)>0

% Input
%  data: a struct variable encoding data information
%    data.I: (estimated) cmi value
%    data.X: T-by-nx, T samples of X
%    data.Y: T-by-ny, T samples of Y
%    data.Z: T-by-nz, T samples of Z
%  par: a struct variable encoding parameter information
%    par.alpha: alpha-level
%    par.ns: number of shuffles to perform

% Output
%  result.ans: binary, 1 (if accepting H1) or 0 (if accepting H0)
%  result.pval: scalar, p-value
%  result.Is: 1-by-ns, shuffled I values (sored from small to large)
%  result.threshold: threshold I value

% Last update/check: Jan 29th 2019
% Reference: 
%   J. Sun, D. Taylor, E. M. Bollt, "Causal network inference by optimal causation entropy"
%   SIAM Journal on Applied Dynamical Systems vol.14, 73-106 (2015).


%% Main
% solution place holders
result.ans = [];
result.pval = [];
result.Is = zeros(1,par.ns);
result.threshold = [];

% cmi values of the shuffled data series
T = size(data.X,1); % number of samples per variable
num = round(par.ns*(par.alpha));
for k = 1 : par.ns
    
    % obtain a shuffled data series for X
    rp = randperm(T);
    Xshuffle = data.X(rp,:);
    % cmi value after the shuffle
    try
    if par.opt == 'Covariance'
        Sshuffle = cov([Xshuffle, data.Y, data.Z]);%corr([Xshuffle, data.Y, data.Z]);
    end
    if par.opt == 'Correlation'
        Sshuffle = corr([Xshuffle, data.Y, data.Z]);
        
    end
    catch
        Sshuffle = cov([Xshuffle, data.Y,data.Z])./diag(cov([Xshuffle, data.Y,data.Z])); %default to covariance
    end
    Sshuffle(Sshuffle<0) = 0;
    
    aa = est_cmi_Poisson_Est(Xshuffle,data.Y,data.Z,Sshuffle);
    [result.Is(k)] = aa;
    
    if length(find(data.I<result.Is))> num
        result.threshold = max(result.Is);
        result.ans = 0;
        result.I = data.I;
        break
    end
end
% from the shuffled cmi values: compute p value etc., and decide outcome
result.Is = sort(result.Is,'ascend');
result.pval = length(find(result.Is>data.I,1))/par.ns;
if k == par.ns
    result.threshold = result.Is(round(par.ns*(1-par.alpha)));
    result.I = data.I;
    if data.I>result.threshold
        result.ans = 1;
    else
        result.ans = 0;
    end
end


