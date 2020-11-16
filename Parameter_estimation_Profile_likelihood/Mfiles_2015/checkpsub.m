function current_osub = checkpsub(parami,stepi,flag_s,Datamat,PLk)
%% DESCRIPTION
% obtain the optimized values of parameters in a subspace (5 in case of the hazard model) 
%% INPUT
% parami            index for the parameter profile likelihod
% stepi             index for the step in exploring the parameter of interest
% flag_s            to check the results with the information on the direction
%                   to explore the profile likelihood (either position or negative)
% Datamat           to check the results with the information on a Data file 
% PLk               the iteration index to compute the profle likelihood
%% OUTPUT
% current_osub      optimized values for the parameters in a subspace
%% History of the version
% 2016-03-28 comments added by Huan Yang 
if ~isunix
    Folder = ['RESULTPL\',Datamat(1:end-4),'_PL',num2str(PLk),'\P',num2str(parami),flag_s,'_s',num2str(stepi),'\'];
else
    Folder = ['RESULTPL/',Datamat(1:end-4),'_PL',num2str(PLk),'/P',num2str(parami),flag_s,'_s',num2str(stepi),'/'];
end;
mat_list = dir([Folder,'*.mat']);
Resnorm_all = [];
theta_est = [];
for mati = 1:length(mat_list)
    loadsave = mat_list(mati).name;
    load([Folder,loadsave],...
        'A_list','k',...
        'R1','len_ind','TS','theta_est_temp','theta0','k_check','output','exitflag','Resnorm_temp');
    Resnorm_all = [Resnorm_all;Resnorm_temp]; %#ok<AGROW>
    theta_est = [theta_est;theta_est_temp]; %#ok<AGROW>
end;
% find the optimal parameter values
ind =  find(Resnorm_all == min(Resnorm_all));
current_osub = theta_est(ind(1),:);