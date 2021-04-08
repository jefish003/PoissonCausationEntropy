function [A] = ER_Adj(N,p);

%Inputs: N- The number of nodes in the ER graph
%        p- The probability of forming an edge between nodes
%
%Output: A- The adjacency matrix

%
%
A = rand(N);

A(A<p) = 1;
A(A<1) = 0;
A = triu(A,1);
A = A + A';