function current_osub=checkpsub(parami,stepi,flag_s,Datamat,PLk)
if ~isunix
    Folder=['RESULTPL\',Datamat(1:end-4),'_PL',num2str(PLk),'\P',num2str(parami),flag_s,'_s',num2str(stepi),'\'];
else
    Folder=['RESULTPL/',Datamat(1:end-4),'_PL',num2str(PLk),'/P',num2str(parami),flag_s,'_s',num2str(stepi),'/'];
end;
mat_list=dir([Folder,'*.mat']);
Resnorm_all=[];
theta_est=[];
for mati=1:length(mat_list)
    loadsave=mat_list(mati).name;
    
    load([Folder,loadsave],...
        'A_list','k',...
        'R1','len_ind','TS','theta_est_temp','theta0','k_check','output','exitflag','Resnorm_temp');
    
    Resnorm_all=[Resnorm_all;Resnorm_temp];
    theta_est=[theta_est;theta_est_temp];
end;
% find the optimal parameter values
ind = find(Resnorm_all==min(Resnorm_all));
current_osub=theta_est(ind(1),:);
