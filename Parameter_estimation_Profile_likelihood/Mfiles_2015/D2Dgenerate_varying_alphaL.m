function Data = D2Dgenerate_varying_alphaL(theta,Datatemplate,Dataseti,flag)
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
    %
    deterministic_param = [theta(1:3),1.5];
    for k = 1:length(Data(TSi).Stimuli)
        stochastic_param = [theta(4)*(1+0.01/60*Data(TSi).Time(k)),theta(5:6)];
        Psi{TSi} = compute_Psi_hazardL_est_end(stochastic_param,deterministic_param,Data(TSi).IPI,Data(TSi).NoP,Data(TSi).PW,Data(TSi).Stimuli(k));
        xi = rand(1,1);
        Data(TSi).Responses(k) = (Psi{TSi}>xi);
    end;
end;