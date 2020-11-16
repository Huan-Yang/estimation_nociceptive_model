function [mPL,mp,msubp]     = mergePL_up(mPL,mp,msubp,nPL,np,nsubp)
%% DESCRITPION
% update the profile likelihood
%% INPUT
% mPL               parameter profile likelihood
% mp                value of the focal parameter
% msubp             values of parameters in the subspace
% nPL
% np
% nsubp
%% OUTPUT
% mPL               merged profile likelihood
% mp                merged value of the focal parameter
% msubp             valus of other parameters after the merge
%% Example
% p1{1} = rand(100,1)*2*pi;
% p1{1} = sort(p1{1});
% y1{1} = sin(p1{1});
% msubp{1} = rand(100,2);
% 
% p2{1} = rand(100,1)*2*pi;
% p2{1} = sort(p2{1});
% y2{1} = cos(p2{1});
% nsubp{1} = rand(100,2);
% 
% figure;
% plot(p1{1},y1{1});
% hold on;
% plot(p2{1},y2{1});
% 
% [y3,p3,usubp]     = mergePL_up(y1,p1,msubp,y2,p2,nsubp);
% plot(p3{1},y3{1},'r','linewidth',3);
%% History of the version
% 2016-03-28 comments added by Huan Yang 

[~,~,~,uPL,up,usubp]  = mergePL(mPL,mp,msubp,nPL,np,nsubp);
[mPL,mp,msubp]        = mergePL(uPL,up,usubp,mPL,mp,msubp);