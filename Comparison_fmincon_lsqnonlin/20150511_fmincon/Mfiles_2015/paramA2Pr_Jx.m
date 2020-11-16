function [Pr,Jx]= paramA2Pr_Jx(A,PW,alpha1,tau1,tau2,alphaL,sigmaL,lambdaL,taus,Mmatrix,PSPnopsum,eralphaL,exptsum)

Qr=compute_Vr(A,PW,tau1,alpha1);
%% compute the PSP and instantaneous firing rate
PSP=PSPnopsum*Qr;                                          % x=Qr.*x^0
temp=-(PSP)/sigmaL;

%     GPUresults=GPUexp(temp);
CPUresults=exp(temp);
exps=min(eralphaL*CPUresults,1e9);
% To prevent numeric error, set the maximal exps as
%     exps=min(exps,1e9);
%
lambdas=1./(1+exps);
%
Jx(1,6)=lambdas*Mmatrix;
lambdaj=Jx(1,6)*lambdaL;
Pr=1-exp(-lambdaj);                 % psychometric function value
%
lambdas2exps=(lambdas.^2.*exps);
Jx(1,4)=-lambdaL/sigmaL*lambdas2exps*Mmatrix;
Jx(1,5)=-lambdaL/(sigmaL.^2)*(lambdas2exps.*(PSP-alphaL))*Mmatrix;
%%
if A*(1-exp(-PW/tau1))-alpha1>=0
    Jx(1,1)=-pi*lambdaL/sigmaL*(lambdas2exps.*PSPnopsum)*Mmatrix;
    Jx(1,2)=A*PW*exp(-PW/tau1)/tau1/tau1*Jx(1,1);
end;
pxtau2=Qr/(tau2^2*(tau2-taus))*exptsum-PSP./(tau2-taus);
Jx(1,3)=lambdaL/sigmaL*(pxtau2.*lambdas2exps)*Mmatrix;