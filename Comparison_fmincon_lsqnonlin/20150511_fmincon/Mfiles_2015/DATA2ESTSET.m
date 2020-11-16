function [thetaset,chi2set,exitflagset]=DATA2ESTSET(lball,uball,Data,Nest,fileSRP)
%%
%%
randseed_HY=0;
rand('state',randseed_HY); %#ok<*RAND>
% options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxIter',1e2,'MaxFunEvals',7*1e2,'Jacobian','on');
% note that for fmincon, it is nessary to provide the alrithom
options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxIter',1e2,'MaxFunEvals',7*1e2,'algorithm','trust-region-reflective','GradObj','on');

lhsmatrix=lhsdesign(Nest,length(lball));
lbmatrix=repmat(lball,Nest,1);
ubmatrix=repmat(uball,Nest,1);
%% generate the Latin hypercubic random numbers for the starting values, in a logscale: for magnitudes
paramatrix = exp(log(lbmatrix)+(log(ubmatrix)-log(lbmatrix)).*lhsmatrix);
% Nk is the number in the Monte Carlo like searching for the initial
% guess values of system parameters
% if ~isunix
%     mkdir (['RESULT\',fileSRP,'\']);
% else
%     mkdir (['RESULT/',fileSRP,'/']);
% end;
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
    for Aj=1:length(A_list{TSi})
        ind_list=find(abs(Data(TSi).Stimuli-A_list{TSi}(Aj))<0.001);
        len_ind{TSi}(Aj)=length(ind_list); % the number of trials at this amplitude
        R1{TSi}(Aj)=sum(Data(TSi).Responses(ind_list)); % the number of trils with 'detected' response
%         [phat{TSi}(Aj),pci{TSi}(Aj,:)] = binofit(R1{TSi}(Aj),length(ind_list));
    end;
    % all the stimulus response pairs are stored in TS, R1 len_ind A_list
end;

for kindex=1:Nest
    disp((kindex));
    %     theta0=paramatrix(kindex,:);
    [thetaset(kindex,:),Resnorm_temp(kindex),exitflagset(kindex,1)]=par_estkernel(Data,fileSRP,randseed_HY,kindex,paramatrix(kindex,:),lball,uball,options,TS,A_list,R1,len_ind);
end;

% chi2set is the chi2 stat.
chi2set(:,1)    = 2*Resnorm_temp;