warning off; %#ok<WNOFF>
%
if ~isunix
    addpath('Mfiles_2015\');
else
    addpath('Mfiles_2015/');
end;
fun4est = @DATA2ESTSET;
% number of the multiple starting values
Nest    = 10;%;
% boundary of the system parameters
lball = [1e-6,    1e-2,   2,      1e-5,   1e-8,   log(2)/2000];
uball = [1,       3,      1e3,    1,      1e-1,   100];             %tau2 and sigmaL are changed
try
    delete(gcp('nocreate'));
catch
end;
% determine the maximum logical core in the computer
n_c=str2double(getenv('NUMBER_OF_PROCESSORS'));
% start the matlab parallelization
if 0
try
    parpool(n_c);
catch
end;
end;