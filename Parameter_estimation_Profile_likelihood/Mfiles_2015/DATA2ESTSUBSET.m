function [thetaset,chi2set,exitflagset]=DATA2ESTSUBSET(lball,uball,Data,Nest,fileSRP,parami,paramf,varargin)
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
%% OUTPUT
% thetaset      set of estimates of parameters from multiple-starting-value optimization
% chi2set       set of chi2 statistics
% exitflagset   set of exit flag values
% varargin      true parameter (optional)
%% History of the version
% 2016-03-25 comments added by Huan Yang 
if length(varargin)==1
    theta_all_true = varargin{1};
else
    theta_all_true = [];
end;

%%
randseed_HY=0;
rand('state',randseed_HY); %#ok<*RAND>
options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxIter',1e2,'MaxFunEvals',7*1e2,'Jacobian','on');
numparam_HM = 6;
subind = setdiff([1:numparam_HM],parami);
lb=lball(subind);
ub=uball(subind);
lhsmatrix=lhsdesign(Nest,length(lb));
lbmatrix=repmat(lb,Nest,1);
ubmatrix=repmat(ub,Nest,1);
%% generate the Latin hypercubic random numbers for the starting values, in a logscale: for magnitudes
paramatrix=exp(log(lbmatrix)+(log(ubmatrix)-log(lbmatrix)).*lhsmatrix);
%% no calibration for applied current amplitudes in Floor's dataset
A_list = cell(length(Data));
TS = zeros(length(Data),3);
len_ind = cell(length(Data));
R1 = cell(length(Data));
for TSi=1:length(Data)
    A_list{TSi}=unique(Data(TSi).Stimuli);
    A_list{TSi}=A_list{TSi}(:);
    TS(TSi,:)=[Data(TSi).NoP,Data(TSi).IPI,Data(TSi).PW];
    len_ind{TSi}=zeros(size(A_list{TSi}));
    R1{TSi}=zeros(size(A_list{TSi}));
    for Aj = 1:length(A_list{TSi})
        ind_list = find(abs(Data(TSi).Stimuli-A_list{TSi}(Aj))<0.001);
        len_ind{TSi}(Aj) = length(ind_list); % the number of trials at this amplitude
        R1{TSi}(Aj) = sum(Data(TSi).Responses(ind_list)); % the number of trils with 'detected' response
        %         [phat{TSi}(Aj),pci{TSi}(Aj,:)]  =  binofit(R1{TSi}(Aj),length(ind_list));
    end;
end;
%% fitting the stimulus-response pairs using the hazard model
% boundary values
%%
for TSi = 1:length(Data)
    N_SRP(TSi) = sum(length(Data(TSi).Stimuli));
end;
N_total_SRP = sum(N_SRP);
theta_true = theta_all_true(setdiff(1:numparam_HM,parami));
paramatrix = [paramatrix;theta_true];
%
% multiple-starting-value optimization for a subspace with five parameters
parfor kindex = 1:size(paramatrix,1)
    disp(kindex)
    [thetaset(kindex,:),Resnorm_temp(kindex),exitflagset(kindex,1)] = par_estsubkernel(Data,paramatrix(kindex,:),lb,ub,options,parami,paramf,TS,A_list,R1,len_ind);
end;
% chi2set is the chi2 stat.
chi2set(:,1)    = 2*Resnorm_temp;