function [FF,JJ]=nln_likelihood_component_J(param,TS,A_list,R1,len_ind, NoP_P,IPI_P,PW_P,A_P,Psi_P)
% Description:
% Both (1)evaluation of the function componnet in the cost function, i.e. the
% sqrt root of the minus log-likelihood and
% (2) the stranghtward compuation of the Jacobian terms
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
    %
    if nargout > 1
        [psi,J_temp]=compute_Psi_hazardL_est_end_J_4(stochastic_param,deterministic_param,IPI,NoP,PW,A);
        
        psi=max(psi,1e-9);
        psi=min(psi,1-1e-9);
        F_temp=sqrt((r1-n10).*log(1-psi)-r1.*log(psi));
        J_temp=1./(2.*repmat(F_temp,1,6)).*repmat(-r1./psi-(r1-n10)./(1-psi),1,6).*J_temp;
        J=[J;J_temp];
        %
    else
        psi=compute_Psi_hazardL_est_end(stochastic_param,deterministic_param,IPI,NoP,PW,A);
        psi=max(psi,1e-9);
        psi=min(psi,1-1e-9);
        F_temp=sqrt((r1-n10).*log(1-psi)-r1.*log(psi));
    end;
    F=[F;F_temp];
end;
FF=sum(F(:).^2);
if nargout > 1
    JJ=zeros(1,6);
    for parami=1:6
        JJ(:,parami) = 2*sum(F(:).*J(:,parami));
    end;
end;