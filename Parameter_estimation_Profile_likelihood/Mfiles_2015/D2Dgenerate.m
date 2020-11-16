function Data = D2Dgenerate(theta,Datatemplate,Dataseti)
%% DESCRITPION
% generate the modeled based dataset with a template dataset from real
% measurements
% note that: the stimulus properpties are the same (adaptive 
% procedeure is not considered here.) 
%% INPUT
% theta             true values of system parameters
% Datatemplate      input and output dataset
% Dataseti          seed of random generator
%% OUTPUT
% Data              generated Dataset
%% History of the version
% 2016-03-25 comments added by Huan Yang 
rand('state',Dataseti);
for TSi = 1:length(Datatemplate)
    Data(TSi).Time = Datatemplate(TSi).Time;
    Data(TSi).Stimuli = Datatemplate(TSi).Stimuli;
    %
    Data(TSi).NoP = Datatemplate(TSi).NoP;
    Data(TSi).IPI = Datatemplate(TSi).IPI;
    Data(TSi).PW = Datatemplate(TSi).PW;
    xi = rand(length(Datatemplate(TSi).Stimuli),1);
    %
    deterministic_param = [theta(1:3),1.5];
    stochastic_param = [theta(4:6)];
    %
    Psi{TSi} = compute_Psi_hazardL_est_end(stochastic_param,deterministic_param,Data(TSi).IPI,Data(TSi).NoP,Data(TSi).PW,Data(TSi).Stimuli);
    Data(TSi).Responses = (Psi{TSi}>xi);
end;