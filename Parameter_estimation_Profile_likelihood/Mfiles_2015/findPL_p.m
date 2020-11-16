function F = findPL_p(p,parami,psub,current_PL,qPLthres,Data)
%% DESCRITPION
% a subrounte to determine the value of the focal paraemeter
%% INPUT
% p                 parameter value (to be tuned)
% parami            parameter index
% psub              values of other parameters
% current_PL        value of current profile likelihood
% qPLthres          target threshold for the profile likelihood
% Data              Dataset with stimuli and binary responses
%% OUTPUT
% F                 objective function
%% History of the version
% 2016-03-26 comments added by Huan Yang 
pall(parami)=p;
pall(setdiff([1:6],parami))=psub;
funhPL=@param2PL;
F = feval(funhPL,pall,Data)-(current_PL+qPLthres);