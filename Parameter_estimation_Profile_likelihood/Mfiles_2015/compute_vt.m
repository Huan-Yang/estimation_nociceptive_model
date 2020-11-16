function vt=compute_vt(t_list,taus,tau2)
%% DESCRIPTION
% compute the impulse response of PSP
%% INPUT
% t_list        set of time points during the stimulation
% taus          time constant of AMPA receptors
% tau2          time constant of
%% OUTPUT
% vt            post-synaptic potential
%% History of the version
% 2016-03-28 comments added by Huan Yang 
vt=1./(tau2-taus).*(exp(-t_list/tau2)-exp(-t_list/taus));