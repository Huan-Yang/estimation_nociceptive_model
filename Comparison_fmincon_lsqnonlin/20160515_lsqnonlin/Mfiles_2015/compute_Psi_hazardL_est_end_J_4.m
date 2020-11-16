function [Pr,Jpsi] = compute_Psi_hazardL_est_end_J_4(stochastic_param,deterministic_param,IPI,NoP,PW,A_list)
%% DESCRIPTION
% computate the psychometric function (electrocutaneous stimulus) and its Jacobian matrix (w.r.t six system parameters) for the
% hazard model
%% INPUT
% stochastic_param      stochatic parameters in the hazard model: alphaL sigmaL lambdaL
% deterministic_param   deterministic parameters in the HM: alpha1 tau1 tau2 taus
% IPI                   interpulse interval
% NoP                   number of pulses
% PW                    pulse width
% A_list                a set of amplutudes for stimuli with temporal properties (IPI,NoP,PW) 
%% OUTPUT
% Pr                    psychometric function
% Jpsi                  Jacobian matrix of Pr with respect to six parameters in the hM
% The Jacobian Jpsi is d(Psi)/d(param) in the hazard model
% Applying the chain rule, Jx is d(lambdaT)/dparam, and
% Jpsi = (exp(-lambdaT))*Jx = (1-Pr)*Jx
%% History of the version
% 2016-03-28 comments added by Huan Yang 
alpha1 = deterministic_param(1);         % threshold
tau1 = deterministic_param(2);           % time constant in the afferent fibres
tau2 = deterministic_param(3);           % time constant in the secondary neuron
taus = deterministic_param(4);           % time constant in the synapse (rising time constant)
T = 2000;%
dt = .1;
t = dt:dt:T;%[time]
% parameters in the hazard/escape model
alphaL = stochastic_param(1);
sigmaL = stochastic_param(2);
lambdaL = stochastic_param(3);
%%
index = round(IPI/dt);
PSPu = compute_vt(t,taus,tau2);
%
PSPnop = zeros(NoP,length(t));
exptnop = zeros(NoP,length(t));
%
PSPnop(1,:) = PSPu;
exptnop(1,:) = t.*exp(-t/tau2);
%%
for k = 2:NoP
    PSPnop(k,1+index*(k-1):end) = PSPu(1:end-index*(k-1));
    exptnop(k,1+index*(k-1):end) = exptnop(1,1:end-index*(k-1));
end;
PSPnopsum = sum(PSPnop,1);                                       % x^0 or x^s
exptsum = sum(exptnop,1);                                        % sum((t-kIPI) exp(-(t-kIPI)/tau2))
%%
lambdaj = zeros(size(A_list));
Pr = lambdaj;
% Jx: d(lambdaT)/d(param)
Jx = zeros(length(A_list),6);            % Jacobian is initialized to be zeros, also for alpha1 tau1 piecewise non-linear expression
% to facilitate a matrix multiplication for speedup,
Mmatrix = ones(length(t),1);
Mmatrix([1,end]) = 1/2;
Mmatrix = Mmatrix*dt;
%
for k = 1:length(A_list)
    A = A_list(k);
    Qr = compute_Vr(A,PW,tau1,alpha1);
    %% compute the PSP and instantaneous firing rate
    PSP = PSPnopsum*Qr;                                          % x=Qr.*x^0
    CPUresults = exp((alphaL-PSP)/sigmaL);
    exps = min(CPUresults,1e12);
    %
    lambdas = 1./(1+exps);
    %
    Jx(k,6) = lambdas*Mmatrix;
    lambdaj(k) = Jx(k,6)*lambdaL;
    Pr(k) = 1-exp(-lambdaj(k));                 % psychometric function value
    %
    lambdas2exps = (lambdas.^2.*exps);
    Jx(k,4) = -lambdaL/sigmaL*lambdas2exps*Mmatrix;
    Jx(k,5) = -lambdaL/(sigmaL.^2)*(lambdas2exps.*(PSP-alphaL))*Mmatrix;
    %%
    if A*(1-exp(-PW/tau1))-alpha1 >= 0
        Jx(k,1) = -pi*lambdaL/sigmaL*(lambdas2exps.*PSPnopsum)*Mmatrix;
        Jx(k,2) = A*PW*exp(-PW/tau1)/tau1/tau1*Jx(k,1);
    end;
    pxtau2 = Qr/(tau2^2*(tau2-taus))*exptsum-PSP./(tau2-taus);
    Jx(k,3) = lambdaL/sigmaL*(pxtau2.*lambdas2exps)*Mmatrix;
end;
Jpsi = Jx.*repmat((1-Pr),1,6);