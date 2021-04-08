function [Hest P] = PoissEnt(X);

%Input: X- An (N x 1) vector of the (independent) entropies you wish to
%          calculate. MUST BE (N x 1)!!!!!
%
%Output: Hest- An Nx1 vector of the (independent) estimated entropies
%
%NOTE: This does NOT calculate the joint entropy of the variable!
%
%==========================================================================
%==========================================================================
%====================== Written By Jeremie Fish ===========================
%==========================================================================
%==================== Last Edited on 11-30-2018 ===========================
%==========================================================================

%Must have a positive rate, so ignoring any negative sign...
if length(X) == 0
    Hest = [];
end
if length(X)>0
Lambda = abs(X);

%Calculate the first probability
First = exp(-X);
Psum = First;
p = containers.Map(0, First);
counter = 0;
%==========================================================================
%====  small is to ensure that this terminates if the probabilities  ======
%====  get really tiny.                                              ======
%==========================================================================

small = 1;
i = 1;
while max(1-Psum)> 1e-16 & small > 1e-75 
    counter = counter+1;
    prob = poisspdf(i,Lambda);%Lambda.*p(counter-1)/counter
    Psum = Psum + prob;
    p(counter) = prob;
    if i >= Lambda
        small = prob;
    end
    i = i+1;
end
%Now that we have the probablities, calculate the entropy
P = cell2mat(values(p)); %convert MAP to matrix for operations
est_a = P.*log(P); %Entropy calculation
est_a(est_a == -Inf) = 0; %Define 0 ln(0) to be 0
est_a(isnan(est_a)) = 0;
est = -sum(est_a,2);
est = real(est);
Hest = est;
%Hest(isnan(Hest))=0; 
end