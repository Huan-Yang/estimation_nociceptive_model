function [current_pm,current_pp,current_PLm,current_PLp]=Data2PL(Data,paramH_index,optimizaed_p,lball,uball,fileSRP,PLk,minchi2)
%%
q 		=.05;
stepN	=500;
param_HM=6;
%
paramset=[1:param_HM];
%
% get the optimized values of parameters
%
funchPL=@param2PL;
PLthres=chi2inv(0.95,1);

% first we search for larger values of the focal parameter
flag_s='p';
for parami=paramH_index
    subparam_ind    =   setdiff(paramset,parami);
    paramf          =   optimizaed_p(parami);
    %
    stepi           =   1;
    flag_c          =   1;
    % step 1 is at the optimal value
    current_pp{parami}(stepi,[subparam_ind])   =optimizaed_p(subparam_ind);
    current_pp{parami}(stepi,parami)= paramf;
    current_PLp{parami}(stepi)  =   feval(funchPL,current_pp{parami}(stepi,:),Data);
    %
    while flag_c
        %% find the focal parameter in the next step
        stepi = stepi + 1;
        if (uball(parami)-current_pp{parami}(stepi-1,parami))./uball(parami)<1e-3
            break;
        end;
        paramf=current2paramf(current_pp{parami}(stepi-1,:),parami,current_PLp{parami}(stepi-1),lball,uball,q,PLthres,Data,flag_s);
        if paramf > uball(parami) || paramf < lball(parami)
            flag_c = 0;
            break;
        end;
        %%
        %         now fix paramf at parami, and re-optimized the other parameters
        seed=20142015;% dummy value
        PL2est1(Data,parami,paramf,seed,stepi-1,flag_s,current_pp{parami}(stepi-1,:),fileSRP,lball,uball,PLk,minchi2);
        %         check the PL at this point
        [current_osub]=checkpsub(parami,stepi-1,flag_s,fileSRP,PLk);
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
%% now we search for smaller values of the focal parameter
flag_s='m';
for parami=paramH_index
    subparam_ind    =   setdiff(paramset,parami);
    paramf          =   optimizaed_p(parami);
    %
    stepi           =   1;
    flag_c          =   1;
    %
    current_pm{parami}(stepi,[subparam_ind])   =optimizaed_p(subparam_ind);
    current_pm{parami}(stepi,parami)= paramf;
    current_PLm{parami}(stepi)  =   feval(funchPL,current_pm{parami}(stepi,:),Data);
    %
    while flag_c
        %% find the focal parameter in the next step
        stepi = stepi + 1;
        %
        if (current_pm{parami}(stepi-1,parami)-lball(parami))./lball(parami)<1e-3 
            break;
        end;
        paramf=current2paramf(current_pm{parami}(stepi-1,:),parami,current_PLm{parami}(stepi-1),lball,uball,q,PLthres,Data,flag_s);
        if paramf > uball(parami) || paramf < lball(parami)
            flag_c = 0;
            break;
        end;
        %%
        %         now fix paramf at parami, and re-optimized the other parameters
        seed=20142015;
        PL2est1(Data,parami,paramf,seed,stepi-1,flag_s,current_pm{parami}(stepi-1,:),fileSRP,lball,uball,PLk,minchi2);
        %         check the PL at this point
        [current_osub]=checkpsub(parami,stepi-1,flag_s,fileSRP,PLk);
        % update the sub-parameters
        current_pm{parami}(stepi,[subparam_ind]) =   current_osub;
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