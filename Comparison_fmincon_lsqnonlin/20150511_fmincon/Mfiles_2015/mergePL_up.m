function [mPL,mp,msubp]     = mergePL_up(mPL,mp,msubp,nPL,np,nsubp)
% history:  bugs are fixed on Feb 6. Feb 16 March 19.
%

%
[~,~,~,uPL1,up1,usubp1] 	= mergePL(mPL,mp,msubp,nPL,np,nsubp);
[~,~,~,uPL2,up2,usubp2]     = mergePL(nPL,np,nsubp,mPL,mp,msubp);
%
paramN=length(mPL);
for parami=1:length(mPL)
    clear A B;
    A(:,1)          =  [up1{parami};up2{parami}];
    if paramN>1
        A(:,2:paramN)   =  [usubp1{parami};usubp2{parami}];
    end;
    A(:,paramN+1)   =  [uPL1{parami};uPL2{parami}];
    %
    B               =  sortrows(A,1);
    mp{parami}    =  B(:,1);
    if paramN>1
        msubp{parami} =  B(:,2:paramN);
    end;
    
    mPL{parami}     =  B(:,paramN+1);
end;