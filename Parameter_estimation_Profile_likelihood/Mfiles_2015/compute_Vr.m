function Vr=compute_Vr(A,PW,tau1,alpha1)
%% DESCRIPTION
% compute the amount of afferent activation
%% INPUT
% A         amplitude 
% PW        pulse width
% tau1      time constant
% alpha1    threshold
%% OUTPUT
% Vr        amount of afferent activation
%% History of the version
% 2016-03-28 comments added by Huan Yang 
L=(A.*(1-exp(-PW./tau1)));
Vr=pi*(L-alpha1).*(L>=alpha1);