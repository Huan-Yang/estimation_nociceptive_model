function Vr=compute_Vr(A,PW,tau1,alpha1)
% the amount of afferent activation
L=(A.*(1-exp(-PW./tau1)));
Vr=pi*(L-alpha1).*(L>=alpha1);