function [thetaset,chi2set,exitflagset] = DATA2ESTSET(lball,uball,Data,Nest,fileSRP,varargin)
%% DESCRIPTION
% mulitple-starting-value optimization for six system parameters in the HM with yes-or-no datasets
%% INPUT
% lball         lower boundary of the parameters to be tuned
% uball         upper boundary of the parameters to be tuned 
% Data          Dataset of the stimuli and yes-or-no responses
% Nest          number of , indeed it is Nest+1, due to one set of 'true parameters'  
% fileSRP       file name containing the measurements, which is used to
%               creat folder to store reults from muliple-starting-value 
%               optimization
% varargin      true parameter (optional)
%% OUTPUT
% thetaset      set of estimates of parameters from multiple-starting-value optimization
% chi2set       set of chi2 statistics
% exitflagset   set of exit flag values
%% History of the version
% 2016-03-25 comments added by Huan Yang 
if length(varargin)==1
    theta_true = varargin{1};
else
    theta_true = [];
end;
%%
randseed_HY = 0;
rand('state',randseed_HY); %#ok<*RAND>
options  =  optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxIter',1e2,'MaxFunEvals',7*1e2,'Jacobian','on',...
    'Algorithm','trust-region-reflective');
lhsmatrix = lhsdesign(Nest,length(lball));
lbmatrix = repmat(lball,Nest,1);
ubmatrix = repmat(uball,Nest,1);
%% generate the Latin hypercubic random numbers for the starting values, in a logscale: for magnitudes
paramatrix = exp(log(lbmatrix)+(log(ubmatrix)-log(lbmatrix)).*lhsmatrix);
%% no calibration for applied current amplitudes in Floor's dataset
A_list = cell(length(Data));
TS = zeros(length(Data),3);
len_ind = cell(length(Data));
R1 = cell(length(Data));

for TSi = 1:length(Data)
    A_list{TSi} = unique(Data(TSi).Stimuli);
    A_list{TSi} = A_list{TSi}(:);
    TS(TSi,:) = [Data(TSi).NoP,Data(TSi).IPI,Data(TSi).PW];
    len_ind{TSi} = zeros(size(A_list{TSi}));
    R1{TSi} = zeros(size(A_list{TSi}));
    for Aj = 1:length(A_list{TSi})
        ind_list = find(abs(Data(TSi).Stimuli-A_list{TSi}(Aj))<0.001);
        len_ind{TSi}(Aj) = length(ind_list); % the number of trials at this amplitude
        R1{TSi}(Aj) = sum(Data(TSi).Responses(ind_list)); % the number of trils with 'detected' response
%         [phat{TSi}(Aj),pci{TSi}(Aj,:)]  =  binofit(R1{TSi}(Aj),length(ind_list));
    end;
    % all the stimulus response pairs are stored in TS, R1 len_ind A_list
end;
%% fitting the stimulus-response pairs using the hazard model
%%
N_SRP = zeros(length(Data),1);
for TSi = 1:length(Data)
   N_SRP(TSi) = sum(length(Data(TSi).Stimuli));
end;
N_total_SRP = sum(N_SRP);
%
paramatrix  =  [paramatrix;theta_true];
% a parallized implementation of multiple-starting-value optimization
for kindex = 1:size(paramatrix,1)
    disp((kindex));
    [thetaset(kindex,:),Resnorm_temp(kindex),exitflagset(kindex,1)] = par_estkernel(Data,fileSRP,randseed_HY,kindex,paramatrix(kindex,:),lball,uball,options,TS,A_list,R1,len_ind);
end;
% chi2set contains (Nest+1) chi2 statistics.
chi2set(:,1)     =  2*Resnorm_temp;