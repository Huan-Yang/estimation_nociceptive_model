function [F,J]=nln_likelihood_component_sub(paramsub,parami,paramivalue,TS,A_list,R1,len_ind)
ind_sub=setdiff([1:length(paramsub)+1],parami);
param(ind_sub)=paramsub;
param(parami)=paramivalue;

stochastic_param=param(4:6);
deterministic_param=[param(1:3),1.5];
F=[];
if nargout > 1
    J=[];
end;
for TSi=1:size(TS,1)
    NoP=TS(TSi,1);
    IPI=TS(TSi,2);
    PW=TS(TSi,3);
    
    A=A_list{TSi};
    r1=R1{TSi};
    n10=len_ind{TSi};
    if nargout > 1
        [psi,J_temp]=compute_Psi_hazardL_est_end_J_4(stochastic_param,deterministic_param,IPI,NoP,PW,A);
        
        psi=max(psi,1e-9);
        psi=min(psi,1-1e-9);
        F_temp=sqrt((r1-n10).*log(1-psi)-r1.*log(psi));
        J_temp=1./(2.*repmat(F_temp,1,6)).*repmat(-r1./psi-(r1-n10)./(1-psi),1,6).*J_temp;
        J=[J;J_temp(:,ind_sub)];
        %
    else
        psi=compute_Psi_hazardL_est_end(stochastic_param,deterministic_param,IPI,NoP,PW,A);
        psi=min(psi,1-1e-4);
        psi=max(psi,1e-4);
        F_temp=sqrt((r1-n10).*log(1-psi)-r1.*log(psi));
    end;
    F=[F;F_temp(:)];
    
end;