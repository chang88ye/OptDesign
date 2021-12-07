warning('off','all')

%% Preprossing
% clear;
close all;


% define target product
targetRxn='EX_succ_e';
% targetRxn='EX_lyco_e';
% targetRxn='EX_narg_e';

if strcmp(targetRxn,'EX_lyco_e')
    load('iML1515_lyco.mat')
end
 
if strcmp(targetRxn,'EX_narg_e')
    load('iML1515_narg.mat')
end

if strcmp(targetRxn,'EX_succ_e')
    sel_model = input('Which model do you want to use, Ecoli_core_model (1) or iML1515 (2) ?');
    if sel_model==1
        % use E. coli core model
        load('ecoli_core_model.mat');     
    elseif sel_model==2
        load('iML1515.mat')
    else
        error('wrong input!')
    end   
end

% get biomass reaction
if isfield(model,'csense')
    if ~iscolumn(model.csense)
        model.csense=columnVector(model.csense);
    end
else
    model.csense=char('E'*ones(length(model.b),1));
end
biomassRxn=model.rxns{model.c==1};

% set uptake rate of oxygen and glucose
oxygenRxn='EX_o2_e';
substrate='EX_glc__D_e';
% model = changeRxnBounds(model,{substrate,oxygenRxn},[-20,-20],{'l','l'});

if strcmpi(model.description,'ecoli_core_model') % identifier is different in ecoli_core model
    targetRxn='EX_succ(e)';
    oxygenRxn='EX_o2(e)';
    substrate='EX_glc(e)';
end

% limit reaction rate in realistic range
model.lb(model.lb<-100)=-100;
model.ub(model.ub>100)=100;

orimodel=model;

% when the model size is big, it is better to compress the model so that
% the linear reactions can be reduced
if size(model.S,1)>100
    % compress the model and get compressed candidate reactions for knockout
    [model,candidate]=preprocessing(orimodel,substrate,oxygenRxn,biomassRxn,targetRxn);
else
    candidate.rxns=model.rxns; 
end

% make sure target reaction is in the reaction list of compressed model
if ~strcmp(model.rxns, targetRxn)
    newTargetRxn=model.rxns{contains(model.rxns, targetRxn)};
else
    newTargetRxn=targetRxn;
end

% remove unreasonable reactions from candidate knockout set
selectedRxns=setdiff(candidate.rxns, {'ATPM', biomassRxn, newTargetRxn});

disp(['Preprocessing step completed, ', num2str(length(selectedRxns)),' candidate reactions are selected for knockout.'])

%% -----------Set up algorithm--------------
% set options for OptDesign algorithm
regBlackList={substrate,'ATPM',biomassRxn, newTargetRxn, 'SUCCtex','SUCCt1pp'};
options.regBlackList=struct('rxnList', {[regBlackList,regBlackList]}, 'typeReg','UUUUUUDDDDDD'); %'EX_h2o','EX_o2','EX_h_e','SUCCtex','SUCCt1pp','GLCt2pp'};
options.selectedRxns=selectedRxns(~contains(selectedRxns,'EX_')); % exclude transport reactions
options.targetRxn=newTargetRxn;
options.minGrowth=0.1; % minimal growth in production mutants
options.minProductRatio=1; % anticipated ratio of production
options.maxKO=5; % at most the number of knockouts
options.maxM=10; % maximal number of manipulations
options.changeThresh={'flux',1}; % noticeable flux change threshold
options.timeLimit=10000;   % a couple of mins


%% ----------- Run algorithm and save results--------------
[solution] = optdesign(model, options,1);

% save results
save(['sol_', targetRxn, '_K',num2str(options.maxKO),'_M', num2str(options.maxM),'.mat'],...
    'solution');

%% visualization
% Draw the flux values on the map "target.svg" which can be opened in FireFox
if strcmpi(orimodel.description,'ecoli_core_model')
    
    % draw the first solution in the pool
    wt_v=solution.pool(1).wt_v;
    mt_v=solution.pool(1).mt_v;
    
    map=readCbMap('ecoli_Textbook_ExportMap');
    options.lb = -10;
    options.ub = 10;
    options.zeroFluxWidth = 0.1;
    options.rxnDirMultiplier = 10;
    options.fileName='WildType.svg';
    options.textsize = 12;
    
    % wild type
    WTsolution=optimizeCbModel(orimodel);
    drawFlux(map, orimodel, wt_v, options);
    web(options.fileName);
    
    % mutant strain
    options.fileName='Mutant.svg';
    drawFlux(map, orimodel, mt_v, options);
    web(options.fileName, '-new');
end

