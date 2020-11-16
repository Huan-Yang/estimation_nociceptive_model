function [Pr,Jpsi]=compute_Psi_hazardL_est_end_J_4_par(stochastic_param,deterministic_param,IPI,NoP,PW,A_list)
%% Description
% The Jacobian Jpsi is d(Psi)/d(param) in the hazard model
% Applying the chain rule, Jx is d(lambdaT)/dparam, and
% Jpsi= (exp(-lambdaT))*Jx=(1-Pr)*Jx

% Pr=zeros(46,1);
% deterministic_param=[alpha1,tau1,tau2,taur,taus];
alpha1=deterministic_param(1);         % threshold
tau1=deterministic_param(2);           % time constant in the afferent fibres
tau2=deterministic_param(3);           % time constant in the secondary neuron
taus=deterministic_param(4);           % time constant in the synapse (rising time constant)
T=2000;%
dt=.1;
t=dt:dt:T;%[time]
% parameters in the hazard/escape model
alphaL=stochastic_param(1);
sigmaL=stochastic_param(2);
lambdaL=stochastic_param(3);
%%
index=round(IPI/dt);
PSPu=compute_vt(t,taus,tau2);
%
PSPnop=zeros(NoP,length(t));
exptnop=zeros(NoP,length(t));
%
PSPnop(1,:)=PSPu;
exptnop(1,:)=t.*exp(-t/tau2);
%%
for k=2:NoP
    PSPnop(k,1+index*(k-1):end)=PSPu(1:end-index*(k-1));
    exptnop(k,1+index*(k-1):end)=exptnop(1,1:end-index*(k-1));
end;
PSPnopsum=sum(PSPnop,1);                                       % x^0 or x^s
exptsum=sum(exptnop,1);                                        % sum((t-kIPI) exp(-(t-kIPI)/tau2))
%%
lambdaj=zeros(size(A_list));
Pr=lambdaj;
% Jx: d(lambdaT)/d(param)
Jx=zeros(length(A_list),6);            % Jacobian is initialized to be zeros, also for alpha1 tau1 piecewise non-linear expression
% to facilitate a matrix multiplication for speedup,
Mmatrix=ones(length(t),1);
Mmatrix([1,end])=1/2;
Mmatrix=Mmatrix*dt;
eralphaL=exp(alphaL/sigmaL);
%
for k=1:length(A_list)
    [Pr(k,1),Jx(k,:)]=paramA2Pr_Jx(A_list(k),PW,alpha1,tau1,tau2,alphaL,sigmaL,lambdaL,taus,Mmatrix,PSPnopsum,eralphaL,exptsum);
end;
Jpsi=Jx.*repmat((1-Pr),1,6);