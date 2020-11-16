function [theta_est_temp,Resnorm_temp,exitflag,output]=par_estkernel(Data,fileSRP,randseed_HY,k,theta0,lball,uball,options,TS,A_list,R1,len_ind)
%
theta_est_temp=theta0;
Resnorm_old  = 1e8;% dummy value
Resnorm_temp = Resnorm_old-1;
%
while Resnorm_old-Resnorm_temp>(1e-7)/2 && Resnorm_old>Resnorm_temp% suggest a changing Chi-square statistics
    Resnorm_old = Resnorm_temp;
    flag_continue=1;
    k_check=0;
    krestart_check=0;
    while flag_continue==1
%      [theta_est_temp,Resnorm_temp,~,exitflag,output] = lsqnonlin(@nln_likelihood_component_J,theta_est_temp,lball,uball,options,TS,A_list,R1,len_ind);
       [theta_est_temp,Resnorm_temp,exitflag,output] = fmincon(@nln_likelihood_component_J,theta_est_temp,...
            [],[],[],[],lball,uball,[],options,TS,A_list,R1,len_ind);
        
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
    fprintf('%4.10f   ',[theta_est_temp,Resnorm_temp*2]);fprintf('\n');
end;