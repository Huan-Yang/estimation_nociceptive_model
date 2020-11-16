clear all;
close all;
clc;
% Here we compare the fmincon and lsqnonlin to this MLE purpose
%
run(['20150511_fmincon','\','AUTO_EST_PL_SRPs.m']);
pause(10);
run(['20160515_lsqnonlin','\','AUTO_EST_PL_SRPs.m']);

%%
close all;
load(['20150511_fmincon','\Mat_fmincon.mat'],'chi2set');
chi2set_fmincon = chi2set;
load(['20160515_lsqnonlin','\Mat_lsqnonlin.mat']);
chi2set_lsqnonlin = chi2set;
%% check whether the sets of optimizated soluation can return identical global minimal for the cost function.
display(['chi2set_fmincon:   ', num2str(min(chi2set_fmincon))]);
display(['chi2set_lsqnonlin: ', num2str(min(chi2set_lsqnonlin))]);
%% show the results
figure('color','w')
plot(chi2set_fmincon,'o');
hold on;
plot(chi2set_lsqnonlin,'x');
set(gca,'yscale','log');