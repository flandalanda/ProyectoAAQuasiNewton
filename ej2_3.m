% Script to perform excercise 2.3
clear;   close all;   clc;

%Define la función DIXMAANA
f=@(x) dixmaana(x);

% n en {240,960}
% m en {1,3,5,17,29}

for n in [240,960]
	x0 = 2*ones([n,1]);

	for m in [1,3,5,17,29]
		sol() = lineLM_BFGS(f,x0,m,tol,maxiter)

	end
end
