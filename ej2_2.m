% Script to perform excercise 2.2
clear;   close all;   clc; format long;

%Define extended Rosenbrock function (generalized banana function)
f=@(x) rosenbrock(x);

%Setting aux variables
TrcSR1=zeros(4,1);
TLM=zeros(4,1);
TBFGS=zeros(4,1);
iterrcSR1=zeros(4,1);
iterLM=zeros(4,1);
iterBFGS=zeros(4,1);
gradfinBFGS=zeros(4,5);
gradfinLM=zeros(4,5);
gradfinrcSR1=zeros(4,5);
rosenBFGS=zeros(4,5);
rosenLM=zeros(4,5);
rosenrcSR1=zeros(4,5);
errBFGS=zeros(4,5);
errLM=zeros(4,5);
errrcSR1=zeros(4,5);


%Running all methods for n={2,8,32,128} and x_0 [-1.2,1,-1.2,....,1]
tic
for i=1:4
    xsize=2^(2*i-1);
    x_0=ones(xsize,1);
    aux=ones(xsize,1);
    aux2=0;
    for j=1:xsize
        if mod(j,2)==1
            x_0(j)=-1.2;
        end 
        
      tic 
        [aux,iterrcSR1(i)]= rcSR1(f, x_0, 200);
      TrcSR1(i)=toc;
      
      rosenrcSR1(i,5)=f(aux);
      errrcSR1(i,5)=norm(ones(xsize)-aux);
      gradfinrcSR1(i,5)=norm(apGrad(f,aux));
      for k=1:4
        [aux,aux2]=rcSR1( f, x_0,iterrcSR1(i)-k);  
        gradfinrcSR1(i,5-k)=norm(apGrad(f,aux));
        rosenrcSR1(i,5-k)=f(aux);
        errrcSR1(i,5-k)=norm(ones(xsize)-aux);
      end
        
        tic
          [aux,iterLM(i)]=lineLM_BFGS( f, x_0, 10^(-5), 100, 128)
        TLM(i)=toc;

          rosenLM(i,5)=f(aux);
          errLM(i,5)=norm(ones(xsize)-aux);  
          gradfinLM(i,5)=norm(apGrad(f,aux));
          for k=1:4
            [aux,aux2]=lineLM_BFGS( f, x_0, 10^(-5),iterLM(i)-k, 128);  
            gradfinLM(i,5-k)=norm(apGrad(f,aux));
            rosenLM(i,5-k)=f(aux);
            errLM(i,5-k)=norm(ones(xsize)-aux);
          end

          tic    
            [aux,iterBFGS(i)]=lineBGFS( f, x_0, 10^(-5), 1000);
          TBFGS(i)=toc;

          rosenBFGS(i,5)=f(aux);
          errBFGS(i,5)=norm(ones(xsize)-aux);
          gradfinBFGS(i,5)=norm(apGrad(f,aux));
          for k=1:4
            [aux,aux2]=lineBGFS( f, x_0, 10^(-5),iterBFGS(i)-k);  
            gradfinBFGS(i,5-k)=norm(apGrad(f,aux));
            rosenBFGS(i,5-k)=f(aux);
            errBFGS(i,5-k)=norm(ones(xsize)-aux);
          end
      
    end    
end
toc
