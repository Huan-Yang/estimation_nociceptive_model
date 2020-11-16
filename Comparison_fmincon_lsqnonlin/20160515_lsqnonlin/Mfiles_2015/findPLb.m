function [ind_PLb,thetasetb,chi2setb,diffCHI2] = findPLb(mPL,mp,thetaset,chi2set)
%% DESCRITPION
% check whether any point below the merged profile likelihood
%% INPUT
% mPL               merged profile likelihood 
% mp                value of the focal parameter of the merged profile likelihood
% thetaset          a set of estimated parameter values from multiple-starting-value scheme
% chi2set           minimal value of the objective function
%% OUTPUT
% ind_PLb           index of parameter values being below 
% thetasetb         parameter values being below the merged profile likelihood
% chi2setb          chi2 statistics for thetasetb
% diffCHI2          mPL - chi2setb
%% History of the version
% 2016-03-29 comments added by Huan Yang 
flag_record = 0;
ind_PLb = []; thetasetb = []; chi2setb = []; diffCHI2 = [];
for thetaseti = 1:size(thetaset,1)
    tmp_theta = thetaset(thetaseti,:);
    chi2check = chi2set(thetaseti);
    if chi2check > min(mPL{1})+chi2inv(0.95,length(mPL))
        continue;
    end;
    for parami = 1:length(mPL)		% for each parameter
        clear tmp_mPL tmp_mp tmp_nPL tmp_np
        tmp_mPL     =   mPL{parami};
        tmp_mp      =   mp{parami};
        lb          =   tmp_mp(1);
        ub          =   tmp_mp(end);
        %
        tmp_thetai = tmp_theta(parami);
        if (tmp_thetai > ub ||tmp_thetai<lb) 
            flag_record = 1;
            break;
        else
            tmp_cPL     = interp1(tmp_mp,tmp_mPL,tmp_thetai,'pchip');
            if tmp_cPL-chi2check>1e-6	% check when the interpolated value from computed PL is really larger than the local minima value 
                flag_record = 1;
                diffCHI2 = [diffCHI2;tmp_cPL-chi2check];		% store the difference of chi2 statics in the diffCHI2
                break;
            end;
        end;     
    end;
    if flag_record
        ind_PLb     = [ind_PLb;thetaseti];
        thetasetb   = [thetasetb;tmp_theta];
        chi2setb    = [chi2setb;chi2check];
		% reset the flag back to 0
        flag_record=0;
    end;
end;