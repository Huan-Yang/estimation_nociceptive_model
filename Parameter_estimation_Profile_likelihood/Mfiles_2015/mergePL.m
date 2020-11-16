function [mPL,mp,msubp,uPL,up,usubp] = mergePL(mPL,mp,msubp,nPL,np,nsubp)
%% DESCRIPTION:
% merge the profile likelihood, i.e. to obtain the lower values of the profile likelihood
%% INPUT
% mPL 	     currently merged profile likelihood
% mp  	     currently merged values of focal parameters
% msubp	     currently merged computed values of sub-parameters
%
% nPL 	     new computed profile likelihood
% np  	     new computed values of focal parameters
% nsubp	     new computed values of sub-parameters
%% OUTPUT
% mPL        merged profile likelihod
% mp         merged parameter value
% msubp	     values of other parameters for merged profile likelihood
% uPL        updated profile likelihood
% up         updated parameter for the profile likelihood
% usubp      updated values of other parameters for the profile likelihood
%% History of the version
% 2016-03-28 comments added by Huan Yang 
% 2016-05-09 remove the bug on the sub dimensional problem by Huan Yang
%%
for parami=1:length(mPL)
    clear tmp_mPL tmp_mp tmp_msubp tmp_nPL tmp_np tmp_nsubp tmp_cPL A B indtemp
    tmp_mPL     =   mPL{parami};  tmp_mPL = tmp_mPL(:);
    tmp_mp      =   mp{parami};   tmp_mp      =   tmp_mp(:);
    tmp_m       =   length(msubp{parami}(:));
    tmp_msubp   =   reshape(msubp{parami},[round(tmp_m/(length(mPL)-1)),5]);
    msubp{parami}=  tmp_msubp;
    
    tmp_nPL     =   nPL{parami};      tmp_nPL = tmp_nPL(:);
    tmp_np      =   np{parami};       tmp_np      =   tmp_np(:);
    tmp_n       =   length(nsubp{parami}(:));
    tmp_nsubp   =   reshape(nsubp{parami},[round(tmp_n/(length(nPL)-1)),5]);    
    nsubp{parami}=  tmp_nsubp;

    msubp{parami}=tmp_msubp;
    if isempty(tmp_mPL)
        mPL{parami}		=nPL{parami}(:);
        mp{parami}		=np{parami}(:);
        msubp{parami}	=nsubp{parami}(:);
        continue;
    end;
    lb          =   tmp_mp(1);
    ub          =   tmp_mp(end);
    %

    indub       =   find(tmp_np > ub);
    indlb       =   find(tmp_np < lb);
    % extend the profile likelihood on the possible boundary values
    parami
    tmp_mPL     = [  tmp_nPL(indlb); tmp_mPL(:); tmp_nPL(indub)];
    tmp_mp      = [  tmp_np(indlb); tmp_mp(:); tmp_np(indub)];
    tmp_msubp    = [  tmp_nsubp(indlb,:); tmp_msubp; tmp_nsubp(indub,:)];
    % which indices need to be updated and check the PL values within the
    % boundary
    ind_interp1 = setdiff(setdiff([1:length(tmp_np)],indub),indlb);
    tmp_cPL     = interp1(tmp_mp, tmp_mPL, tmp_np(ind_interp1), 'cubic');
    %
    updateind{parami} = find(tmp_cPL >= tmp_nPL(ind_interp1)) + length(indlb);
    %
    tmp_mPL     = [tmp_nPL(updateind{parami}); tmp_mPL(:)];
    tmp_mp      = [tmp_np(updateind{parami}); tmp_mp(:)];
    tmp_msubp   = [tmp_nsubp(updateind{parami},:); tmp_msubp];
	
    up{parami}  = tmp_np(updateind{parami});
    uPL{parami} = tmp_nPL(updateind{parami});
    usubp{parami} = tmp_nsubp(updateind{parami},:);
	%
    up{parami}  = [tmp_np(indlb); up{parami}; tmp_np(indub)];
    uPL{parami} = [tmp_nPL(indlb); uPL{parami}; tmp_nPL(indub)];
    usubp{parami} = [tmp_nsubp(indlb,:); usubp{parami}; tmp_nsubp(indub,:)];
    %
    A           = [tmp_mPL, tmp_mp, tmp_msubp];
    B           = sortrows(A,2);
    mp{parami}  = unique(B(:,2));
    for tempi = 1:length(mp{parami})
        indtemp 			= find(B(:,2)==mp{parami}(tempi));
        mPL{parami}(tempi)	= min(B(indtemp,1));
		ind_sub				= find(B(indtemp,1)==min(B(indtemp,1)));
        msubp{parami}(tempi,:) = B((ind_sub(1)),3:end);
    end;
end;