function [P] = PermutationMatrix(p);

%Input: p- This relates to the size of the Permutation matrix P (as 
%          discussed in the paper: 'A local Poisson Graphical Model For 
%          for Inferring Networks From Sequencing Data by Genevera I. Allen
%          and Zhandong Liu), see the paper for details of what P is.
%
%Output P- The permutation matrix P as discussed in the paper listed above
%
%==========================================================================
%==========================================================================
%======================= Written By Jeremie Fish ==========================
%======================== Last Edited 12-05-2018 ==========================
%==========================================================================
%==========================================================================

Zp = p*(p-1)/2; %This will be the number of columns  of P
P = zeros(p,Zp); %Initialize the Permutation matrix P.

%==========================================================================
%========================== Description ===================================
% We need to be able to map subscripts like 12, 13, 14, 23, 24, .... to the
% permutation matrix. The first column of the permutation matrix represents
% 12 (which is the coupling between variable 1 and variable 2) and so on,
% thus the first p columns will be 12, 13, ..., 1p. The next columns will
% be 23,24,...,2p and so on. To get the proper subsripts we use nchoosek
% below. 
%==========================================================================
NK = nchoosek(1:p,2);  
%Now that we have the proper subsripts we can map them to P in a loop.
for i = 1:p
    Z = zeros(1,Zp);
    %We need to know the row (r) that the subscript occured in to map to P
    [r c] = find(NK == i);
    Z(r) = 1;
    P(i,:) = Z;
end
    