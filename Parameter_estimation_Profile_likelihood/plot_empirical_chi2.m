clear all;
close all;
clc;
if isunix
    load('validity_PL/_tmp_MSO_150.mat');
else
    load('validity_PL\_tmp_MSO_150.mat');
end;
newM = PLc(1:150,:)-repmat(LL(1:150)',[6,1])';
figure;
param_list = {'\alpha_1','\tau_1','\tau_2','\alpha_L','\sigma_L','\lambda_L'};
for parami = 1:6
    clear sample x f;
    sample= newM(:,parami);
    tmp = sample(find(sample>0));
    [f,x] = ecdf(tmp);
    
    subplot(2,3,parami);
    plot(x,f)
    hold on;
    xlist = [0.001:0.001:6];
    plot(xlist,chi2cdf(xlist,1));
    title(param_list(parami));
    xlabel('LPL(\theta^*)- 2LL');
    ylabel('ECDF');
end;
hl = legend('Empirical','\chi_2');