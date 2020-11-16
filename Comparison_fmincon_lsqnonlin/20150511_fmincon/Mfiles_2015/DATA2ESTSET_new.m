function [thetaset,chi2set,exitflagset]=DATA2ESTSET(lball,uball,Data,Nest,fileSRP)
%%
numparam_LC =8;
numparam_HM =6;
%%
randseed_HY=1;
rand('state',randseed_HY); %#ok<*RAND>
options = optimset('Display','off','TolFun',1e-10,'TolX',1e-10,'MaxIter',5e3,'MaxFunEvals',7*5e3,'Jacobian','on');
lhsmatrix=lhsdesign(Nest,length(lball));
lbmatrix=repmat(lball,Nest,1);
ubmatrix=repmat(uball,Nest,1);
%% generate the Latin hypercubic random numbers for the starting values, in a logscale: for magnitudes
paramatrix=exp(log(lbmatrix)+(log(ubmatrix)-log(lbmatrix)).*lhsmatrix);
% Nk is the number in the Monte Carlo like searching for the initial
% guess values of system parameters
if ~isunix
    mkdir (['RESULT\',fileSRP,'\']);
    check_mat = dir(['RESULT\',fileSRP,'\*.mat']);
else
    mkdir (['RESULT/',fileSRP,'/']);
    check_mat = dir(['RESULT/',fileSRP,'/*.mat']);
end;

results_ind_list=[];
for result_matfilei=1:length(check_mat)
    if ~isunix
        ind_temp_name=find(check_mat(result_matfilei).name=='_');
    else
        ind_temp_name=find(check_mat(result_matfilei).name=='_');
    end;
    results_ind_list=[results_ind_list;str2num( check_mat(result_matfilei).name(ind_temp_name(3)+1:ind_temp_name(4)-1) )];
end;
%
%% no calibration for applied current amplitudes in Floor's dataset
for TSi=1:length(Data)
    A_list{TSi}=unique(Data(TSi).Stimuli);
    A_list{TSi}=A_list{TSi}(:);
    TS(TSi,:)=[Data(TSi).NoP,Data(TSi).IPI,Data(TSi).PW];
    len_ind{TSi}=zeros(size(A_list{TSi}));
    R1{TSi}=zeros(size(A_list{TSi}));
    for Aj=1:length(A_list{TSi})
        ind_list=find(abs(Data(TSi).Stimuli-A_list{TSi}(Aj))<0.001);
        len_ind{TSi}(Aj)=length(ind_list); % the number of trials at this amplitude
        R1{TSi}(Aj)=sum(Data(TSi).Responses(ind_list)); % the number of trils with 'detected' response
        [phat{TSi}(Aj),pci{TSi}(Aj,:)] = binofit(R1{TSi}(Aj),length(ind_list));
    end;
    % all the stimulus response pairs are stored in TS, R1 len_ind A_list
end;
%% fitting the stimulus-response pairs using the hazard model
% boundary values
%%
for TSi=1:length(Data)
    [bglm(:,TSi),~,stats{TSi}] =glmfit(A_list{TSi},[R1{TSi} len_ind{TSi}],'binomial','link','logit');
    [bestGLMfit{TSi}] = glmval(bglm(:,TSi),A_list{TSi},'logit',stats{TSi},'size',len_ind{TSi});
    Pr_LC{TSi}=bestGLMfit{TSi}./len_ind{TSi};
    logL_LC_component(TSi)=sum(log(Pr_LC{TSi}).*R1{TSi}+(len_ind{TSi}-R1{TSi}).*log(1-Pr_LC{TSi}));
    N_SRP(TSi)=sum(length(Data(TSi).Stimuli));
end;
N_total_SRP=sum(N_SRP);
logL_fitted_LC=sum(logL_LC_component);

AIC_LC=-2*logL_fitted_LC+2*numparam_LC;
BIC_LC=-2*logL_fitted_LC+log(N_total_SRP)*numparam_LC;
for kindex=fliplr(setdiff(1:Nest,results_ind_list))
    disp((kindex));
    [thetaset(kindex,:),Resnorm_temp(kindex),exitflagset(kindex,1)]=par_estkernel(Data,fileSRP,randseed_HY,kindex,paramatrix(kindex,:),lball,uball,options,TS,A_list,R1,len_ind);
end;
% chi2set is the chi2 stat.
chi2set(:,1)    = 2*Resnorm_temp;