function [PLp,PL,PLsubp]=Theta2PL(Data,thetaG,lball,uball,fileSRP,k,minchi2)
% Description
% 1. compute the profile likelihood based on the given values of thetaG 
% 2. also the parameter in their subspace 
paramHM=6;
[current_pm,current_pp,current_PLm,current_PLp]   = Data2PL(Data,1:paramHM,thetaG,lball,uball,fileSRP,k,minchi2);

PLp=cell(paramHM,1);
PL=cell(paramHM,1);
PLsubp=cell(paramHM,1);

for parami=1:paramHM
    % focal parameter values
    PLp{parami} = [flipud(current_pm{parami}(2:end,parami));current_pp{parami}(1:end,parami)]; 
    % -2 times of the profile likelihood, PL is just a place-holder
    PL{parami}  = [fliplr(current_PLm{parami}(2:end)),current_PLp{parami}(1:end)]';
    %
    PLsubp{parami}=[flipud(current_pm{parami}(2:end,setdiff([1:6],parami)));current_pp{parami}(1:end,setdiff([1:6],parami))];
end;