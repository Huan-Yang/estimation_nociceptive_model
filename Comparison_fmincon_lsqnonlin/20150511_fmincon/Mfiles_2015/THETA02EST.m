function [theta_est_temp,chi2PL,exitflag]=THETA02EST(theta0,Data,lball,uball)
options = optimset('Display','off','TolFun',1e-12,'TolX',1e-12,'MaxIter',1e4,'MaxFunEvals',1e4,'Jacobian','on');
%% no calibration for applied current amplitudes in Floor's dataset
for TSi=1:length(Data)
    A_list{TSi}=unique(Data(TSi).Stimuli);
    A_list{TSi}=A_list{TSi}(:);
    TS(TSi,:)=[Data(TSi).NoP,Data(TSi).IPI,Data(TSi).PW];
    len_ind{TSi}=zeros(size(A_list{TSi}));
    R1{TSi}=zeros(size(A_list{TSi}));
    for Aj=1:length(A_list{TSi})
        ind_list=find(abs(Data(TSi).Stimuli-A_list{TSi}(Aj))<0.001);
        len_ind{TSi}(Aj)=length(ind_list); % the number of trials at this amplitude
        R1{TSi}(Aj)=sum(Data(TSi).Responses(ind_list)); % the number of trials with 'detected' response
%         [phat{TSi}(Aj),pci{TSi}(Aj,:)] = binofit(R1{TSi}(Aj),length(ind_list));
    end;
    % all the stimulus response pairs are stored in TS, R1 len_ind A_list
end;
%% fitting the stimulus-response pairs using the hazard model
try
    % generate a random set of the initial guess values
    flag_continue=1;
    k_check=0;
    krestart_check=0;
    theta0_temp=theta0;
    while flag_continue==1
        [theta_est_temp,Resnorm_temp,~,exitflag,output] = lsqnonlin(@nln_likelihood_component_J,theta0_temp,lball,uball,options,TS,A_list,R1,len_ind);
        theta0_temp=theta_est_temp;
        if exitflag==1      % already find the 'local optimal'
            flag_continue=0; % stop
        elseif exitflag==0  % the iteration limit was reached in the optimization
            krestart_check=krestart_check+1;
            if krestart_check>2
                flag_continue=0;
            else
                theta0_temp=theta_est_temp;
                flag_continue=1;
            end;
        elseif exitflag>1
            k_check=k_check+1;
            if k_check>2
                flag_continue=0;
            else
                theta0_temp=theta_est_temp;
                flag_continue=1;
            end;
        end;
    end;
    chi2PL=Resnorm_temp*2;
end;