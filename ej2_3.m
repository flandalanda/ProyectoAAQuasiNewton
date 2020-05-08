% Script to perform excercise 2.3
clear;   close all;   clc;

% Handle for DIXMAANA function
f=@(x) dixmaana(x);

% Parameters n (DIXMAANA) and m (LM BFGS)
n = [240,960];
m = [1,3,5,17,29];

% Aux variables
norma = [];
valor = [];
tiempo = [];
iter = [];

for i = 1:length(n) 
	x0 = 2*ones([n(i),1]);
	for j = 1:length(m)
        tic
		[xk, k] = lineLM_BFGS(f, x0, 10^(-5), 1000, m(j));
        tiempo = [tiempo, toc];        
        norma = [norma, norm(apGrad(f,xk))];
        valor = [valor, f(xk)];
        iter = [iter, k];
	end
end
