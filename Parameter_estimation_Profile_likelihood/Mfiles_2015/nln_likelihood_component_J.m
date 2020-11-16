function [F,J] = nln_likelihood_component_J(param,TS,A_list,R1,len_ind)
%% DESCRIPTION
% Evaluation of the psychometric function in the hazard model (HM), and the evaluation of the Jacobian matrix (optional)
%% INPUT 
% param         six parameter of the HM
% TS            conbinations of the temporal stimulus parameter
% A_list        sets of stimulus amplitudes
% R1            sets of 'yes' responses corresonding to the stimulus properties 
% len_ind       sets of the numbers of stimuli with specific stimulus properties
%% OUTPUT   
% F             sets of evaluated components in the log-likelihood term
%               sqrt((r1-n10).*log(1-psi)-r1.*log(psi))
% J             Jacobian matrix of F with respection to subspace parameters
%% History of the version
% 2016-03-26 comments added by Huan Yang 
stochastic_param = param(4:6);
deterministic_param = [param(1:3),1.5];
F = [];
if nargout > 1
    J = [];
end;
for TSi = 1:size(TS,1)
    NoP = TS(TSi,1);
    IPI = TS(TSi,2);
    PW = TS(TSi,3);
    
    A = A_list{TSi};
    r1 = R1{TSi};
    n10 = len_ind{TSi};
    %
    if nargout > 1
        [psi,J_temp] = compute_Psi_hazardL_est_end_J_4(stochastic_param,deterministic_param,IPI,NoP,PW,A);
        psi = max(psi,1e-9);
        psi = min(psi,1-1e-9);
        F_temp = sqrt((r1-n10).*log(1-psi)-r1.*log(psi));
        J_temp = 1./(2.*repmat(F_temp,1,6)).*repmat(-r1./psi-(r1-n10)./(1-psi),1,6).*J_temp;
        J = [J;J_temp];
        %
    else
        psi = compute_Psi_hazardL_est_end(stochastic_param,deterministic_param,IPI,NoP,PW,A);
        psi = max(psi,1e-9);
        psi = min(psi,1-1e-9);
        F_temp = sqrt((r1-n10).*log(1-psi)-r1.*log(psi));
    end;
    F = [F;F_temp(:)];
end;