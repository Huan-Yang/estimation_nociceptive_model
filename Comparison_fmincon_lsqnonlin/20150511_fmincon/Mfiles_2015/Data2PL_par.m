function [current_pm,current_pp,current_PLm,current_PLp]=Data2PL_par(Data,paramH_index,optimizaed_p,lball,uball,fileSRP,PLk,minPL)
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
parami=paramH_index
subparam_ind    =   setdiff(paramset,parami);
paramf          =   optimizaed_p(parami);
%
stepi           =   1;
flag_c          =   1;
% step 1 is at the optimal value
current_pp(stepi,[subparam_ind])   =optimizaed_p(subparam_ind);
current_pp(stepi,parami)= paramf;
current_PLp(stepi)  =   feval(funchPL,current_pp(stepi,:),Data);
%
while flag_c
    %% find the focal parameter in the next step
    stepi = stepi + 1;
    paramf=current2paramf(current_pp(stepi-1,:),parami,current_PLp(stepi-1),lball,uball,q,PLthres,Data,flag_s);
    if paramf > uball(parami) || paramf < lball(parami)
        flag_c = 0;
        break;
    end;
    %%
    %         now fix paramf at parami, and re-optimized the other parameters
    seed=20142015;
    PL2est1(Data,parami,paramf,seed,stepi-1,flag_s,current_pp(stepi-1,:),fileSRP,lball,uball,PLk);
    %         check the PL at this point
    [current_osub]=checkpsub(parami,stepi-1,flag_s,fileSRP,PLk);
    % update the sub-parameters
    current_pp(stepi,[subparam_ind]) =   current_osub;
    current_pp(stepi,parami)         =   paramf;
    current_PLp(stepi)               =   feval(funchPL,current_pp(stepi,:),Data);
    %
    if stepi > stepN
        flag_c = 0;
    elseif current_PLp(stepi)>minPL+chi2inv(0.95,param_HM)
        flag_c = 0;
    end;
end;
%% now we search for smaller values of the focal parameter
flag_s='m';
parami=paramH_index
subparam_ind    =   setdiff(paramset,parami);
paramf          =   optimizaed_p(parami);
%
stepi           =   1;
flag_c          =   1;
%
current_pm(stepi,[subparam_ind])   =optimizaed_p(subparam_ind);
current_pm(stepi,parami)= paramf;
current_PLm(stepi)  =   feval(funchPL,current_pm(stepi,:),Data);
%
while flag_c
    %% find the focal parameter in the next step
    stepi = stepi + 1;
    paramf=current2paramf(current_pm(stepi-1,:),parami,current_PLm(stepi-1),lball,uball,q,PLthres,Data,flag_s);
    if paramf > uball(parami) || paramf < lball(parami)
        flag_c = 0;
        break;
    end;
    %%
    %         now fix paramf at parami, and re-optimized the other parameters
    seed=20142015;
    PL2est1(Data,parami,paramf,seed,stepi-1,flag_s,current_pm(stepi-1,:),fileSRP,lball,uball,PLk);
    %         check the PL at this point
    [current_osub]=checkpsub(parami,stepi-1,flag_s,fileSRP,PLk);
    % update the sub-parameters
    current_pm(stepi,[subparam_ind]) =   current_osub;
    current_pm(stepi,parami)         =   paramf;
    current_PLm(stepi)               =   feval(funchPL,current_pm(stepi,:),Data);
    %
    if stepi > stepN
        flag_c = 0;
    elseif current_PLm(stepi)>minPL+chi2inv(0.95,param_HM)
        flag_c = 0;
    end;
end;
