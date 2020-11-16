function F=compute_residual_pr_Hazard(A,stochastic_param,deterministic_param,IPI,NoP,PW,pr0)
% this function is used in root-finding of detection threshold:
% psi(A)-0.5=0
[Pr]=compute_Psi_hazardL_est_end(stochastic_param,deterministic_param,IPI,NoP,PW,A);
alpha1=deterministic_param(1);
tau1=deterministic_param(2);
%%
Ac=(1-exp(-PW./tau1))^-1*alpha1;
% check the existness of the perception thresholds
if compute_Vr(Ac,PW,deterministic_param(2),deterministic_param(1))-pr0<0
    %% exist
else
    %% non-exist
end;
%%
Qr=compute_Vr(A,PW,deterministic_param(2),deterministic_param(1));
%%
if Qr==0
    F=(A-Ac)+compute_Vr(Ac,PW,deterministic_param(2),deterministic_param(1))-pr0;
else
F=Pr-pr0;
end;