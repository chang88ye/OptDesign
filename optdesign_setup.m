function optdesign_setup()

% launch cobratoolbox
try 
    initCobraToolbox
catch
    error(' > cobratoolbox is not available.\n');
end

warning on
% use Gurobi (if installed) as the default solver for LP, QP and MILP problems
solOK=changeCobraSolver('gurobi', 'ALL', 0);
if ~solOK
    warning('>GUROBI is not fully accessible, please check gurobi installation or switch to another MILP solver./n')
end

% make sure the current folder is OptDesign
str = which('optdesign.m');
if isempty(str)
    error('Current directory is not NIHBA, please change it to OptDesign.\n')
end

warning off
% add OptDesign to matlab search path
curPath=fileparts(str);
addpath(genpath(curPath));
end