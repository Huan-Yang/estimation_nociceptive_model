clear all;
close all;
clc;
%%
% load the boundary values of the estimates
configure_EST_PL;
Nest = 10;
% pick one measurement SRPs dataset
if ~isunix
    pathSRP ='..\';
    fileSRP = 'D9450_1_0.mat';
    load([pathSRP,'\',fileSRP]);
else
    pathSRP ='../';
    fileSRP = 'D9450_1_0.mat';
    load([pathSRP,'/',fileSRP]);
end;
% STEP 1: MULTIPLE STARTING VALUE OPTIMIZATION
[thetaset,chi2set,exitflagset] = feval(fun4est,lball,uball,Data,Nest,fileSRP);
t_lsqnonlin = toc;
save('Mat_lsqnonlin.mat','thetaset','chi2set','exitflagset','t_lsqnonlin');
