clear all; %#ok<CLALL>
close all;
clc;
%% to check the validaty of the profile-likelihood based confidence interval
warning off; %#ok<WNOFF>
% add the path of matlab scripts
if ~isunix
    addpath('Mfiles_2015\');
else
    addpath('Mfiles_2015/');
end;
%%
if ~isunix
    dirSRPfiles     = 'DATA\Floor\Day00_Floor\';
else
    dirSRPfiles     = 'DATA/Floor/Day00_Floor/';
end;
%% pick the data file
filename = uigetfile([dirSRPfiles,'*.mat']);
load([dirSRPfiles,filename]);
Datatemplate = Data;
clear Data;
% number of the multiple starting values for optimization
DatasetN = 200;
% pick the maximum likelihood estimate of parameters
% NOTE that the file should correspond to the data file
if ~isunix
    [thetaGname,thetaGpath] = uigetfile(['BIC\','*.mat']);
    load([thetaGpath,'\',thetaGname],'thetaG');
else
    [thetaGname,thetaGpath] = uigetfile(['BIC/','*.mat']);
    load([thetaGpath,'/',thetaGname],'thetaG');
end;
theta_true = thetaG;
try
    delete(gcp('nocreate'));
catch
end;
% Likelihood, profile likelihood for each true value of system parameter
LL      = zeros(DatasetN,1);
PLc     = zeros(DatasetN,1);
% difference between PLc and LL
dist    = zeros(DatasetN,1);
timecost= zeros(DatasetN,1);
%
Dataset_interval = 1:50;
DATA = cell(DatasetN,1);
chi2set = cell(DatasetN,1);
chi2setPL = cell(DatasetN,1);
thetaset = cell(DatasetN,6);
thetasetPL = cell(DatasetN,5);

for Dataseti = Dataset_interval
    % generate the dataset using the hazard model and true values
    DATA{Dataseti} = D2Dgenerate(theta_true,Datatemplate,Dataseti);
    %%
    disp(['compute the LL: ',num2str(Dataseti)]);
    configure_EST_PL;
    [thetaset{Dataseti},chi2set{Dataseti}] = feval(fun4est,lball,uball,DATA{Dataseti},Nest,num2str(Dataseti),theta_true);
    LL(Dataseti) =  min(chi2set{Dataseti});
    %% for PL of each of six parameters in the HM
    for parami = 1:6
        disp(['compute the LL: ', num2str(Dataseti),'  ',num2str(parami)]);
        tic;
        configure_EST_PL;
        fun4est = @DATA2ESTSUBSET;
        [thetasetPL{Dataseti,parami},chi2setPL{Dataseti,parami}] = feval(fun4est,lball,uball,DATA{Dataseti},Nest,num2str(Dataseti),parami,theta_true(parami),theta_true);
        %
        PLc(Dataseti,parami) 	=   min(chi2setPL{Dataseti,parami});
        %     dist(Dataseti) = PLc(Dataseti) - LL(Dataseti);
        timecost(Dataseti,parami) = toc;
    end;
    if ~isunix
        save(['validity_PL\_tmp_MSO_',num2str(Dataseti)],'PLc','DATA','timecost','chi2setPL','theta_true','Datatemplate','LL');
    else
        save(['validity_PL/_tmp_MSO_',num2str(Dataseti)],'PLc','DATA','timecost','chi2setPL','theta_true','Datatemplate');
    end;
end;