function [theta_est_temp,Resnorm_temp,exitflag,output]=par_estsubkernel(Data,theta0,lb,ub,options,parami,paramivalue,TS,A_list,R1,len_ind)
%% DESCRIPTION
% local optimization for the HM (five parameters) with yes-or-no response dataset
%% INPUT
% Data          Dataset of the stimuli and yes-or-no responses
% fileSRP       for store results purpose
% randseed_HY   for store results purpose
% k             for store results purpose    
% theta0        starting values of parameters (five parameters)
% lb            lower boundary of the parameters to be tuned
% ub            upper boundary of the parameters to be tuned 
% options       options for the optimization
% parami        index of parameter, which is a contant
% paramivalue   value of the parameter with the index parami 
% TS            conbinations of the temporal stimulus parameter
% A_list        sets of stimulus amplitudes
% R1            sets of 'yes' responses corresonding to the stimulus properties 
% len_ind       sets of the numbers of stimuli with specific stimulus properties
%% OUTPUT
% thetaset      set of estimates of parameters from multiple-starting-value optimization
% chi2set       set of chi2 statistics
% exitflagset   set of exit flag values
% output        optimization summary
%% History of the version
% 2016-03-24 comments added by Huan Yang 
%%
theta_est_temp = theta0;
Resnorm_old  = 1e8;% dummy value
Resnorm_temp = Resnorm_old-1;
%
while Resnorm_old-Resnorm_temp>(1e-7)/2 && Resnorm_old>Resnorm_temp% suggest a changing Chi-square statistics
    Resnorm_old = Resnorm_temp;
    flag_continue = 1;
    k_check = 0;
    krestart_check = 0;
    while flag_continue == 1         % paramsub,parami,paramivalue,TS,A_list,R1,len_ind
        [theta_est_temp,Resnorm_temp,~,exitflag,output] = lsqnonlin(@nln_likelihood_component_sub,theta_est_temp,lb,ub,options,parami,paramivalue,TS,A_list,R1,len_ind);
        if exitflag == 1      % already find the 'local optimal'
            flag_continue = 0; % stop
        elseif exitflag == 0  % the iter limite was reached in the optimization
            krestart_check=krestart_check+1;
            if krestart_check>2
                flag_continue=0;
            else
                flag_continue=1;
            end;
        elseif exitflag>1
            %             disp('due to reaching the tolerance');
            k_check=k_check+1;
            if k_check>2
                flag_continue=0;
            else
                flag_continue=1;
            end;
        end;
    end;
%     fprintf('%4.10f   ',[theta_est_temp,Resnorm_temp*2]);fprintf('\n');
end;