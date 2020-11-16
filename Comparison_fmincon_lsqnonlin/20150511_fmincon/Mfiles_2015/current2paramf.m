function  paramf=current2paramf(current_p,parami,current_PL,lb,ub,q,PLthres,Data,flag_s)
p       =   current_p(parami);
psub    =   current_p(setdiff([1:6],parami));
ubp     =   ub(parami);
lbp     =   lb(parami);
%
k=1;
kN=5;
%
options=optimset('Display','off');
while k<kN
    if flag_s=='p'
        [pnext,fq]=lsqnonlin(@findPL_p,p,p,ubp,options,parami,psub,current_PL,PLthres*q,Data);
    elseif flag_s=='m'
        [pnext,fq]=lsqnonlin(@findPL_p,p,lbp,p,options,parami,psub,current_PL,PLthres*q,Data);
    end;
    q=q/2;
    k=k+1;
    if fq<1e-6		% check whether chi2(next)-chi2(now)=chi2thres*q
        paramf=pnext;
        return;
    end;
end;
%
if k>=kN
    if flag_s=='p'
        paramf=p*1.02;
    elseif flag_s=='m'
        paramf=p/1.02;
    end;
end;