function [theta_est_temp,chi2PL,exitflag] = THETA02EST(theta0,Data,lball,uball)
%% DESCRIPTION
% estimate the parameter with a set of starting values for the system parameters in the hazard model 
%% INPUT
% theta0            parameters to be tuned for an optimization purpose
% Data              Dataset with stimuli and binary responses
% lball             lower boundaries of the system parameters
% uball             upper boundaries of the system parameters
%% OUTPUT
% theta_est_temp    estimated values of parameters
% chi2PL            chi2 statistics for the fitted values
% exitflag          exit flag of the optimizer
%% History of the version
% 2016-03-28 comments added by Huan Yang 
options  =  optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxIter',1e2,'MaxFunEvals',7*1e2,'Jacobian','on');
%% no calibration for applied current amplitudes in Floor's dataset
A_list = cell(length(Data),1);
TS = zeros(length(Data),3);
len_ind = cell(length(Data),1);
R1 = cell(length(Data),1);
for TSi = 1:length(Data)
    A_list{TSi} = unique(Data(TSi).Stimuli);
    A_list{TSi} = A_list{TSi}(:);
    TS(TSi,:) = [Data(TSi).NoP,Data(TSi).IPI,Data(TSi).PW];
    len_ind{TSi} = zeros(size(A_list{TSi}));
    R1{TSi} = zeros(size(A_list{TSi}));
    for Aj = 1:length(A_list{TSi})
        ind_list = find(abs(Data(TSi).Stimuli-A_list{TSi}(Aj))<0.001);
        len_ind{TSi}(Aj) = length(ind_list); % the number of trials at this amplitude
        R1{TSi}(Aj) = sum(Data(TSi).Responses(ind_list)); % the number of trials with 'detected' response
        %         [phat{TSi}(Aj),pci{TSi}(Aj,:)]  =  binofit(R1{TSi}(Aj),length(ind_list));
    end;
    % all the stimulus response pairs are stored in TS, R1 len_ind A_list
end;
%
randseed_HY = 0;
fileSRP = [];
k = 0;
[theta_est_temp, Resnorm_temp, exitflag, output] = par_estkernel(Data,fileSRP,randseed_HY,k,theta0,lball,uball,options,TS,A_list,R1,len_ind);
chi2PL = 2*Resnorm_temp;