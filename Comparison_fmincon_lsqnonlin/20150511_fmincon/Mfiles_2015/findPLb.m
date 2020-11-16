function [ind_PLb,thetasetb,chi2setb,diffCHI2] = findPLb(mPL,mp,thetaset,chi2set)
% check whether any point below the merged profile likelihood
flag_record=0;
ind_PLb=[];thetasetb=[];chi2setb=[];diffCHI2=[];
for thetaseti=1:size(thetaset,1)
    tmp_theta=thetaset(thetaseti,:);
    chi2check=chi2set(thetaseti);
    if chi2check>min(mPL{1})+ chi2inv(0.95,length(mPL))
        continue;
    end;
    for parami=1:length(mPL)		% for each parameter
        clear tmp_mPL tmp_mp tmp_nPL tmp_np tmp_mpold tmp_mPLold
        tmp_mPLold     =   mPL{parami};
        tmp_mpold      =   mp{parami};
        %
        [tmp_mp]=   unique(tmp_mpold); 
        tmp_mPLp    =   [];
        for unique_tmpi=1:length(tmp_mp)
            tmp_mPLp(unique_tmpi)=min(tmp_mPLold(find(tmp_mpold==tmp_mp(unique_tmpi))));
        end;
        tmp_mPL=tmp_mPLp;
        %
        lb          =   tmp_mp(1);
        ub          =   tmp_mp(end);
        %
        tmp_thetai = tmp_theta(parami);
        if (tmp_thetai>ub ||tmp_thetai<lb) 
            flag_record=1;
            break;
        else
            tmp_cPL     = interp1(tmp_mp,tmp_mPL,tmp_thetai,'pchip');
            if tmp_cPL-chi2check>1e-6	% check when the interpolated value from computed PL is really larger than the local minima value 
                flag_record=1;
                diffCHI2=[diffCHI2;tmp_cPL-chi2check];		% store the difference of chi2 statics in the diffCHI2
                break;
            end;
        end;     
    end;
    if flag_record
        ind_PLb     =[ind_PLb;thetaseti];
        thetasetb   =[thetasetb;tmp_theta];
        chi2setb    =[chi2setb;chi2check];
		% reset the flag back to 0
        flag_record=0;
    end;
end;