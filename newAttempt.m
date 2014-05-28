clear all;close all;clc;
%%
N = 3;

TN = 2*eye(N) - diag(ones(N-1,1),1) - diag(ones(N-1,1),-1);
TNxN = kron(eye(N),TN) + kron(TN,eye(N));
reshape(magic(2),[1,2])