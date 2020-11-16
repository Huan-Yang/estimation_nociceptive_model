clear all; %#ok<CLALL>
close all;
clc;
%%
warning off; %#ok<WNOFF>
%%
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
fileSRP = uigetfile([dirSRPfiles,'*.mat']);
load([dirSRPfiles,fileSRP]);
try
    delete(gcp('nocreate'));
catch
end;
%
configure_EST_PL;
Nest    = 100;%;
% step 1: multiple starting value optimization
[thetaset,chi2set] = feval(fun4est,lball,uball,Data,Nest,fileSRP);
%% compute the profile likelihood
% STEP 2: PICK THE BEST FIT from multiple-starting optimization
minchi2 	=   min(chi2set);
indmchi2	=   find(chi2set==minchi2);
thetaG  	=   thetaset(indmchi2(1),:);
%
flag_global = 1; %1: there is a change in the global optimization
flag_PLsub  = 1; %1: there is a set of parameter values with chi2 value below

%
try
    delete(gcp('nocreate'));
catch
end;
% determine the maximum logical core in the computer
n_c=str2double(getenv('NUMBER_OF_PROCESSORS'));
% start the matlab parallelization
try
    parpool(max(n_c,6));
catch
end;
%
k       = 1;	%  the step index of running profile likelihood
[PLp{k},PL{k},PLsubp{k}] = Theta2PL(Data,thetaG,lball,uball,fileSRP,k,minchi2);
mPL     = PL{k};
mPLp    = PLp{k};
mPLsubp = PLsubp{k};
figure(1)
for parami = 1:length(lball)
    subplot(1, length(lball), parami);
    hold on;
    plot(mPLp{parami},mPL{parami},'b-','linewidth',3);
    set(gca,'xscale','log');
    ylim([min(mPL{parami})-0.5,min(mPL{parami})+10]);
    xlim([lball(parami),uball(parami)]);
end;
%%
while (flag_global+flag_PLsub)>0
    % flag_global: the estimate in thetaG is the global optimum
    %
    disp(k);
    % STEP 3: COMPUTE THE (LOCAL) PROFILE LIKELIHOOD with thetaG
    chi2PL=inf;
    for parami=1:6
        if chi2PL>min(PL{k}{parami})
            chi2PL  = min(PL{k}{parami});
            ind_theta_a=find(PL{k}{parami}==chi2PL);
            theta_a(parami) 	= PLp{k}{parami}(ind_theta_a(1));
            theta_a(setdiff([1:6],parami)) =PLsubp{k}{parami}(ind_theta_a(1),:);
        end;
    end;
    % plot the (merged) profile likelihood
    % STEP 4: check whether there exists a even better set of values
    % check whether there is any set of values of parameters during the PL computation with a even better fit (a smaller chi2 statics)
    if minchi2-chi2PL > 1e-4
        thetaG0  = theta_a; % the even-better candidate
        % a further optimization from this even better candidate as the starting value
        % the thetaG and minchi2 are updated
        [thetaG,minchi2] = THETA02EST(thetaG0,Data,lball,uball);
        k       = k+1;
        %
        % compute the PL in a new round with the 'best' initial value;
        [PLp{k},PL{k},PLsubp{k}]=Theta2PL(Data,thetaG,lball,uball,fileSRP,k,minchi2);
        [mPL,mPLp,mPLsubp]     = mergePL_up(mPL,mPLp,mPLsubp,PL{k},PLp{k},PLsubp{k});
        % show the merged PL results
        figure(1)
        for parami=1:6
            subplot(1,6,parami);
            hold on;
            plot(mPLp{parami},mPL{parami},'bx--');
            set(gca,'xscale','log');
            ylim([min(mPL{parami})-0.5,min(mPL{parami})+10]);
            xlim([lball(parami),uball(parami)]);
        end;
    end;
    % after global optima is confirmed
    if ~isunix
        try
            mkdir(['_tmp\',fileSRP(1:end-4)]);
        catch
        end;
        save(['_tmp\',fileSRP(1:end-4),'\_tmp_mergedPL_',fileSRP(1:end-4),num2str(k)],'mPL','mPLp','mPLsubp','PLp','PL','PLsubp','thetaG','minchi2');
    else
        try
            mkdir(['_tmp/',fileSRP(1:end-4)]);
        catch
        end;
        save(['_tmp/',fileSRP(1:end-4),'/_tmp_mergedPL_',fileSRP(1:end-4),num2str(k)],'mPL','mPLp','mPLsubp','PLp','PL','PLsubp','thetaG','minchi2');
    end;
    flag_global = 0;
    flag_PLsub =1;
    while flag_PLsub
        % STEP 5: Check whether there exists a set of values in step 1 with a chi2 statics below the merged profile likelihood.
        [ind_PLb,thetasetb,chi2setb,diffCHI2] = findPLb(mPL,mPLp,thetaset(setdiff([1:size(thetaset,1)],indmchi2),:),chi2set(setdiff([1:size(thetaset,1)],indmchi2)));
        % only pick the index with the most detectable difference
        % between computed PL and local minima from step 1
        indPL=find(diffCHI2==max(diffCHI2)); %
        %
        fprintf('How many theta are left in the multiple-starting value optimization procedure: %d \n',length(ind_PLb)); %#ok<CTPCT
        if isempty(indPL)
            flag_PLsub  = 0;
            disp('PL computation finished!!!');
            % stored outcome are mPL and thetaG
        else
            % show the local minima together with the merged PL
            figure(1);
            for parami=1:6
                subplot(1,6,parami);
                hold on;
                scatter(thetasetb(:,parami),chi2setb,50,'fill');
                ylim([min(mPL{parami})-0.5,min(mPL{parami})+10]);
                xlim([lball(parami),uball(parami)]);
            end;
            %
            k       = k+1;
            theta_temp  = thetasetb((indPL(1)),:);
            [PLp{k},PL{k},PLsubp{k}]=Theta2PL(Data,theta_temp,lball,uball,fileSRP,k,minchi2);
            [mPL,mPLp,mPLsubp]     = mergePL_up(mPL,mPLp,mPLsubp,PL{k},PLp{k},PLsubp{k});
        end;
        if ~isunix
            save(['_tmp\',fileSRP(1:end-4),'\_tmp_mergedPL_',fileSRP(1:end-4),num2str(k)],'mPL','mPLp','mPLsubp','PLp','PL','PLsubp','thetaG','minchi2');
        else
            save(['_tmp/',fileSRP(1:end-4),'/_tmp_mergedPL_',fileSRP(1:end-4),num2str(k)],'mPL','mPLp','mPLsubp','PLp','PL','PLsubp','thetaG','minchi2');
        end;
    end;
end;