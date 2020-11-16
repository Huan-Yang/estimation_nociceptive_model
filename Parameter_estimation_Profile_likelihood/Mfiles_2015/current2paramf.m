function  paramf = current2paramf(current_p,parami,current_PL,lb,ub,q,PLthres,Data,flag_s)
%% DESCRITPION
% determine the value of the focal parameter in computing the profile
% likelihood
%% INPUT
% current_p     current value of the parameter
% parami        index of the system parameter
% current_PL    current value of the profile likelihood
% lb            lower boundary of the six parameters
% ub            upper boundary of the six parameters
% q             a proportional factor
% PLthres       a pre-fixed threshold for profile likelihood
% Data          Dataset with stimuli and responses 
% flag_s        directional flag: 'p' for position, 'm' for negative 
%% OUTPUT
% paramf        value of the focal parameter
%% History of the version
% 2016-03-28 comments added by Huan Yang 
p       =   current_p(parami);
psub    =   current_p(setdiff([1:6],parami));
ubp     =   ub(parami);
lbp     =   lb(parami);
%
k=1;
kN=5;
%
options = optimset('Display','off');
while k<kN
    if flag_s == 'p'
        [pnext,fq] = lsqnonlin(@findPL_p,p,p,ubp,options,parami,psub,current_PL,PLthres*q,Data);
    elseif flag_s == 'm'
        [pnext,fq]=lsqnonlin(@findPL_p,p,lbp,p,options,parami,psub,current_PL,PLthres*q,Data);
    end;
    q = q/2;
    k = k+1;
    if fq < 1e-6		% check whether chi2(next)-chi2(now)=chi2thres*q
        paramf = pnext;
        return;
    end;
end;
%
if k >= kN
    if flag_s == 'p'
        paramf = p*1.02;
    elseif flag_s == 'm'
        paramf = p/1.02;
    end;
end;