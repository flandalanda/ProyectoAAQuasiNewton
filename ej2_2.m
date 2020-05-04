% Script to perform excercise 2.2
clear;   close all;   clc;

%Define extende Rosenbrock function (generalized banana function)
f=@(x) rosenbrock(x);

%Setting aux variables
TrcSR1=zeros(4,1);
TLM=zeros(4,1);
TBFGS=zeros(4,1);
iterrcSR1=zeros(4,1);
iterLM=zeros(4,1);
iterBFGS=zeros(4,1);


%Running all methods for n={2,8,32,128} and x_0 [-1.2,1,-1.2,....,1]

for i=1:4
    xsize=2^(2*i-1);
    x_0=ones(xsize,1);
    res=ones(xsize,1);
    
    for j=1:xsize
        if mod(j,2)==1
            x_0(j)=-1.2;
        end 
%     tic
%       [res,iterLM(i)]=lineLM_BFGS( f, x_0, 10^(-5), 100 )
%     TLM(i)=toc;
      tic    
        [res,iterBFGS(i)]=lineBGFS( f, x_0, 10^(-5), 500);
      TBFGS(i)=toc;
      tic 
        [res,iterrcSR1(i)]= rcSR1(f, x_0, 200);
      TrcSR1(i)=toc;
    end    
end    