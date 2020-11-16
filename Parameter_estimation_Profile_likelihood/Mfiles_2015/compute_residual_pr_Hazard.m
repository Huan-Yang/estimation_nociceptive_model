function F = compute_residual_pr_Hazard(A,stochastic_param,deterministic_param,IPI,NoP,PW,pr0)
%% DESCRIPTION
%  root-finding of detection threshold A50: F(A50): =  psi(A50)-pr0  =  0
%% INPUT
% A                     stimulus amplutude [mA]
% stochastic_param      stochastic parameters: alphaL sigmaL lambdaL
% deterministic_param   parameters: alpha1 tau1 tau2 taus
% IPI                   interpulse interval
% NoP                   number of pulses
% PW                    pulse width
% pr0                   correponsind probability for a detection threshold (0.50 by default)
%% OUTPUT
% F                     compute the difference between model-based psychometric functin Pr and Pr0
%% History of the version
% 2016-03-28 comments added by Huan Yang 
[Pr] = compute_Psi_hazardL_est_end(stochastic_param,deterministic_param,IPI,NoP,PW,A);
alpha1 = deterministic_param(1);
tau1 = deterministic_param(2);
%%
Ac = (1-exp(-PW./tau1))^-1*alpha1;
% check the existness of the perception thresholds
if compute_Vr(Ac,PW,deterministic_param(2),deterministic_param(1))-pr0<0
    %% exist
else
    %% non-exist
end;
%%
Qr = compute_Vr(A,PW,deterministic_param(2),deterministic_param(1));
%%
if Qr == 0
    F = (A-Ac)+compute_Vr(Ac,PW,deterministic_param(2),deterministic_param(1))-pr0;
else
F = Pr-pr0;
end;