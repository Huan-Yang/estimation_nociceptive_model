function vt=compute_vt(t_list,taus,tau2)
% compute the impulse response of PSP
vt=1./(tau2-taus).*(exp(-t_list/tau2)-exp(-t_list/taus));