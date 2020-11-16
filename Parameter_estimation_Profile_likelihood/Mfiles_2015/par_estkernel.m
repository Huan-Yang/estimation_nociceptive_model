function [theta_est_temp,Resnorm_temp,exitflag,output]=par_estkernel(Data,fileSRP,randseed_HY,k,theta0,lball,uball,options,TS,A_list,R1,len_ind)
%% DESCRIPTION
% local optimization for the HM (six parameters) with yes-or-no response dataset
%% INPUT
% Data          Dataset of the stimuli and yes-or-no responses
% fileSRP       for store results purpose
% randseed_HY   for store results purpose
% k             for store results purpose    
% theta0        starting values of parameters (six parameters)
% lball         lower boundary of the parameters to be tuned
% uball         upper boundary of the parameters to be tuned 
% options       options for the optimization
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
theta_est_temp = theta0;
Resnorm_old  = 1e8;% dummy value
Resnorm_temp = Resnorm_old-1;
%
while Resnorm_old-Resnorm_temp>(1e-7)/2 && Resnorm_old>Resnorm_temp% suggest a changing Chi-square statistics
    Resnorm_old = Resnorm_temp;
    flag_continue = 1;
    k_check=0;
    krestart_check=0;
    while flag_continue==1
        [theta_est_temp,Resnorm_temp,~,exitflag,output] = lsqnonlin(@nln_likelihood_component_J,theta_est_temp,lball,uball,options,TS,A_list,R1,len_ind);
        if exitflag==1      % already find the 'local optimal'
            flag_continue=0; % stop
        elseif exitflag==0  % the iter limite was reached in the optimization
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
%     compute the log-likelihood using the estimated values of
%     parameters in the HM
% logL_fitted_HM = -Resnorm_temp;
% %% best-fit prediction using the 'local optimal values'
% for TSi=1:length(Data)
%     stochastic_param=theta_est_temp(4:6);
%     deterministic_param=[theta_est_temp(1:3),1.5];
%     NoP=TS(TSi,1);
%     IPI=TS(TSi,2);
%     PW=TS(TSi,3);
%     [Pr_HM{TSi}]=compute_Psi_hazardL_est_end(stochastic_param,deterministic_param,IPI,NoP,PW,A_list{TSi});
% end;
%% save the computation results
% c=clock;
% yy='0000';mm='00';dd='00';hh='00';minutes='00';seconds='00';
% yy=num2str(c(1));
% mm(end-length(num2str(c(2)))+1:end)=num2str(c(2));
% dd(end-length(num2str(c(3)))+1:end)=num2str(c(3));
% hh(end-length(num2str(c(4)))+1:end)=num2str(c(4));
% minutes(end-length(num2str(c(5)))+1:end)=num2str(c(5));
% seconds(end-length(num2str(round(c(6))))+1:end)=num2str(round(c(6)));
% dataandtime=[yy,mm,dd,hh,minutes,seconds];
% %%
% if ~isunix
%     save(['RESULT\',fileSRP,'\','s_',num2str(randseed_HY),'_k_',num2str(k),'_',dataandtime],...
%         'A_list','k','R1','len_ind','TS','theta_est_temp','theta0','k_check','output','exitflag','Resnorm_temp');
% else
%     save(['RESULT/',fileSRP,'/','s_',num2str(randseed_HY),'_k_',num2str(k),'_',dataandtime],...
%         'A_list','k','R1','len_ind','TS','theta_est_temp','theta0','k_check','output','exitflag','Resnorm_temp');
% end;