function [current_pm,current_pp,current_PLm,current_PLp] = Data2PL(Data,paramH_index,optimized_p,lball,uball,fileSRP,PLk,minchi2)
%% DESCRIPTION
%% INPUT 
% Data              Dataset with stimuli and binary responses
% paramH_index      index of the focal parameter
% optimized_p       the maximum likelihood estimate
% lball             lower boundary of the system parameters
% uball             upper boundary of the system parameters
% fileSRP           name of the Dataset
% PLk               profile likelihood at k-th iteration
% minchi2           value of the cost function with the parameter thetaG
%% OUTPUT   
% current_pm        current value of the parameter in the negative direction
% current_pp        current value of the parameter in the positive direction
% current_PLm       current value of the profile likelihood in the negative direction
% current_PLp       current value of the profile likelihood in the positive direction
%% History of the version
% 2016-03-25 comments added by Huan Yang 

% a step constant to explore the profile likelihood
q 		 = .2;
% the total number of maximum steps
stepN	 = 500;
param_HM = 6;
%
paramset = 1:param_HM;
%
% get the optimized values of parameters
%
funchPL = @param2PL;
PLthres = chi2inv(0.95,1);
%
% first search for larger values of the focal parameter
current_pp = cell(param_HM,1);
current_PLp = cell(param_HM,1);
flag_s = 'p';
parfor parami = paramH_index
    subparam_ind    =   setdiff(paramset,parami);
    paramf          =   optimized_p(parami);
    %
    stepi           =   1;
    flag_c          =   1;
    % step 1 is at the optimal value
    current_pp{parami}(stepi,subparam_ind) = optimized_p(subparam_ind);
    current_pp{parami}(stepi,parami) = paramf;
    current_PLp{parami}(stepi) = feval(funchPL,current_pp{parami}(stepi,:),Data);
    %
    while flag_c
        %% find the focal parameter in the next step
        stepi = stepi + 1;
        if (uball(parami)-current_pp{parami}(stepi-1,parami))./uball(parami)<1e-3
            break;
        end;
        paramf = current2paramf(current_pp{parami}(stepi-1,:),parami,current_PLp{parami}(stepi-1),lball,uball,q,PLthres,Data,flag_s);
        if paramf > uball(parami) || paramf < lball(parami)
            flag_c = 0;
            break;
        end;
        %%
        %         now fix paramf at parami, and re-optimized the other parameters
        seed = 20142015;% dummy value
        PL2est1(Data,parami,paramf,seed,stepi-1,flag_s,current_pp{parami}(stepi-1,:),fileSRP,lball,uball,PLk,minchi2);
        %         check the PL at this point
        [current_osub] = checkpsub(parami,stepi-1,flag_s,fileSRP,PLk);
        % update the sub-parameters
        current_pp{parami}(stepi,[subparam_ind]) =   current_osub;
        current_pp{parami}(stepi,parami)         =   paramf;
        current_PLp{parami}(stepi)               =   feval(funchPL,current_pp{parami}(stepi,:),Data);
        %
        if stepi > stepN
            flag_c = 0;
        elseif current_PLp{parami}(stepi)>minchi2+chi2inv(0.95,param_HM)
            flag_c = 0;
        end;
    end;
end;
%% now search for smaller values of the focal parameter
current_pm = cell(param_HM,1);
current_PLm = cell(param_HM,1);
flag_s = 'm';
parfor parami=paramH_index
    subparam_ind    =   setdiff(paramset,parami);
    paramf          =   optimized_p(parami);
    %
    stepi           =   1;
    flag_c          =   1;
    %
    current_pm{parami}(stepi,[subparam_ind]) = optimized_p(subparam_ind);
    current_pm{parami}(stepi,parami) = paramf;
    current_PLm{parami}(stepi)  =   feval(funchPL,current_pm{parami}(stepi,:),Data);
    %
    while flag_c
        %% find the focal parameter in the next step
        stepi = stepi + 1; % increase in the iteration variable
        %
        if (current_pm{parami}(stepi-1,parami)-lball(parami))./lball(parami)<1e-3 
            break;
        end;
        paramf = current2paramf(current_pm{parami}(stepi-1,:),parami,current_PLm{parami}(stepi-1),lball,uball,q,PLthres,Data,flag_s);
        if paramf > uball(parami) || paramf < lball(parami)
            flag_c = 0;
            break;
        end;
        %%
        %  now fix paramf at parami, and re-optimized the other parameters
        seed = 20142015;
        PL2est1(Data,parami,paramf,seed,stepi-1,flag_s,current_pm{parami}(stepi-1,:),fileSRP,lball,uball,PLk,minchi2);
        %         check the PL at this point
        [current_osub] = checkpsub(parami,stepi-1,flag_s,fileSRP,PLk);
        % update the sub-parameters
        current_pm{parami}(stepi,[subparam_ind]) = current_osub;
        current_pm{parami}(stepi,parami)         =   paramf;
        current_PLm{parami}(stepi)               =   feval(funchPL,current_pm{parami}(stepi,:),Data);
        %
        if stepi > stepN
            flag_c = 0;
        elseif current_PLm{parami}(stepi)>minchi2+chi2inv(0.95,param_HM)
            flag_c = 0;
        end;
    end;
end;