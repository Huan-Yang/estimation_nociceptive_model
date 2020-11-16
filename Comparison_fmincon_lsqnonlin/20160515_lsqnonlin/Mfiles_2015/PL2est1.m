function         current_osub = PL2est1(Data,parami,paramf,seed,stepi,flag_s,current_p,Datamat,lball,uball,PLk,minchi2,flag_store)
%% DESCRIPTION
% compute the profile likelihood and store the results of profile likelihood
%% INPUT
% Data          Dataset with stimuli and binary responses
% parami        index of parameter
% paramf        value of the fixed parameter
% seed          seed of the random generator
% stepi         step index for the iterated exploration of the profile likelihood
% flag_s        direction indictor of either positive or negative
% current_p     current values of the sytem parameters
% Datamat       filename of the Data file for a restore purpuse
% lball         lower boundary of the system parameters
% uball         upper boundary of the system parameters
% PLk           profile likeihood at the k-th iteration
% minchi2       the minimal value of the chi2 statistics
% flag_store    flag to store the result
%% OUTPUT
% no output, because all the results are restored in some folders
%% History of the version
% 2016-03-26 comments added by Huan Yang 
numparam_HM = 6;
%% load HM-simulated dataset with SRPs
rand('state',seed); %#ok<*RAND>
subind =  setdiff([1:numparam_HM],parami);
% obtain the boundaries of the supspace with indeices subind
lb = lball(subind);
ub = uball(subind);
% options for the profile likelihood
options  =  optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxIter',1e2,'MaxFunEvals',7*1e2,'Jacobian','on');
for TSi = 1:length(Data)
    A_list{TSi} = unique(Data(TSi).Stimuli);
    A_list{TSi} = A_list{TSi}(:);
    TS(TSi,:) = [Data(TSi).NoP,Data(TSi).IPI,Data(TSi).PW];
    %
    len_ind{TSi} = zeros(size(A_list{TSi}));
    R1{TSi} = zeros(size(A_list{TSi}));
    for Aj = 1:length(A_list{TSi})
        ind_list = find(abs(Data(TSi).Stimuli-A_list{TSi}(Aj))<0.001);
        len_ind{TSi}(Aj) = length(ind_list); % the number of trials at this amplitude
        R1{TSi}(Aj) = sum(Data(TSi).Responses(ind_list)); % the number of trils with 'detected' response
    end;
    % all the stimulus response pairs are stored in TS, R1 len_ind A_list
end;

if ~isunix
    Folder = ['RESULTPL\',Datamat(1:end-4),'_PL',num2str(PLk),'\P',num2str(parami),flag_s,'_s',num2str(stepi),'\'];
else
    Folder = ['RESULTPL/',Datamat(1:end-4),'_PL',num2str(PLk),'/P',num2str(parami),flag_s,'_s',num2str(stepi),'/'];
end;
if flag_store
    mkdir(Folder);
end;
%
k = 0;
try
    % generate a random set of the initial guess values
    theta0 = current_p(subind);
    flag_continue = 1;
    k_check = 0;
    krestart_check = 0;
    theta0_temp = theta0;
    Resnorm_temp_pre = inf;
    % a further optimization with the starting value (theta0_temp)
    Resnorm_temp = norm(feval(@nln_likelihood_component_sub,theta0_temp,parami, paramf,TS,A_list,R1,len_ind));
    while Resnorm_temp_pre-Resnorm_temp>1e-7/2 && Resnorm_temp_pre>Resnorm_temp % check whether the change in Resnorm_temp is small
        Resnorm_temp_pre = Resnorm_temp;
        while flag_continue == 1
            % optimize the likelihood with five free parameters together
            % with a fixed parameter
            [theta_est_temp,Resnorm_temp,~,exitflag,output]  =  lsqnonlin(@nln_likelihood_component_sub,theta0_temp,lb,ub,options,parami,paramf,TS,A_list,R1,len_ind);
            theta0_temp = theta_est_temp;
            %
            if exitflag == 1      % already find the 'local optimum'
                flag_continue = 0; % stop
            elseif exitflag == 0  % the iteration limit was reached in the optimization
                krestart_check = krestart_check+1;
                if krestart_check>2
                    flag_continue = 0;
                else
                    theta0_temp = theta_est_temp;
                    flag_continue = 1;
                end;
            elseif exitflag > 1  % due to small change in the function value or in the scearching variable
                k_check = k_check+1;
                if k_check > 2   % just for a double-check purpose
                    flag_continue = 0;
                    k_check = 0;
                else
                    theta0_temp = theta_est_temp;
                    flag_continue = 1;
                end;
            end;
        end;
    end;
    %% save the computation results
    if flag_store
        c = clock;
        yy = '0000';mm = '00';dd = '00';hh = '00';minutes = '00';seconds = '00';
        
        yy = num2str(c(1));
        mm(end-length(num2str(c(2)))+1:end) = num2str(c(2));
        dd(end-length(num2str(c(3)))+1:end) = num2str(c(3));
        
        hh(end-length(num2str(c(4)))+1:end) = num2str(c(4));
        minutes(end-length(num2str(c(5)))+1:end) = num2str(c(5));
        seconds(end-length(num2str(round(c(6))))+1:end) = num2str(round(c(6)));
        
        dataandtime = [yy,mm,dd,hh,minutes,seconds];
        %%
        save([Folder,'s_',num2str(seed),'_k_',num2str(k),'_',dataandtime],...
            'A_list','k',...
            'R1','len_ind','TS','theta_est_temp','theta0','k_check','output','exitflag','Resnorm_temp','paramf','parami');
    end;
    fprintf('p. index: %d, step: %d, p: %8.8f, chi2 stat.: %8.8f, min.chi2: %8.8f\n', parami,stepi+1,paramf,2*Resnorm_temp,minchi2);
    current_osub = theta_est_temp;
    %     [current_osub,paramf] =
    %     disp([num2str(parami),'    ',num2str(2*Resnorm_temp),'    ',num2str(paramf)]);
end;