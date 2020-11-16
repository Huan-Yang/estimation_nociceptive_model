function [Pr,lambdaj] = compute_Psi_hazardL_est_end(stochastic_param,deterministic_param,IPI,NoP,PW,A_list)
%% DESCRIPTION
% computate the psychometric function (electrocutaneous stimulus) for the hazard model
% the difference to compute_Psi_hazardL_est_end_J is that no Jacobian
% matrix returned
%% INPUT
% stochastic_param      stochatic parameters in the hazard model: alphaL sigmaL lambdaL
% deterministic_param   deterministic parameters in the HM: alpha1 tau1 tau2 taus
% IPI                   interpulse interval
% NoP                   number of pulses
% PW                    pulse width
% A_list                a set of amplutudes for stimuli with temporal properties (IPI,NoP,PW) 
%% OUTPUT
% Pr                    psychometric function
% lambdaj               expected number of spikes (lambdaT) as a function of Time
%% History of the version
% 2016-03-28 comments added by Huan Yang 
alpha1 = deterministic_param(1);         % threshold
tau1 = deterministic_param(2);           % time constant in the afferent fibres
tau2 = deterministic_param(3);           % time constant in the secondary neuron
taus = deterministic_param(4);           % time constant in the synapse (rising time constant)
T = 2000;% stimulus duration [ms]
dt = .1;
t = dt:dt:T; %[time]
alphaL = stochastic_param(1);
sigmaL = stochastic_param(2);
lambdaL = stochastic_param(3);
%%
index = round(IPI/dt);
PSPu = compute_vt(t,taus,tau2);
%
PSPnop = zeros(NoP,length(t));
%
PSPnop(1,:) = PSPu;
%%
for k = 2:NoP
    PSPnop(k,1+index*(k-1):end) = PSPu(1:end-index*(k-1));
end;
PSPnopsum = sum(PSPnop,1);
%%
lambdaj = zeros(size(A_list));
Pr = lambdaj;
for k = 1:length(A_list)
    A = A_list(k);
    Qr = compute_Vr(A,PW,tau1,alpha1);
    %% compute the PSP and instantaneous firing rate
    PSP = PSPnopsum*Qr;
    %%
    lambda = 1./(1+exp(-1/sigmaL.*(PSP-alphaL)));
    lambdaj(k) = dt*(sum(lambda(2:end-1))+sum(lambda([1,end]))/2);
    Pr(k) = 1-exp(-lambdaL.*lambdaj(k));
end;