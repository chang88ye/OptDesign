% an example of buiding butanol biosythetic pathway based on iML1515 and then using optDesign to find design strategies.

load('iML1515.mat')

But_1 = 'nad_c + btcoa_c <=> nadh_c + b2coa_c + h_c'; % called ACOAD1 in B. subtilis (iYO844)
But_2 = 'nad_c + coa_c + btal_c <=> nadh_c + btcoa_c + h_c'; % 
But_3 = 'nadh_c + btal_c + h_c -> nad_c + btOH_c'; % make this reaction irreversible, otherwise btOH will be a carbon source for growth.
But_4 = 'btOH_c <=>';

model = addReaction(model,'R01171',But_1);
model = addReaction(model,'R01172',But_2);
model = addReaction(model,'R03544',But_3);
model = addReaction(model,'EX_ButOH',But_4);

model.description='iML1515_ButOH';
model.csense(end+1:length(model.b))='E';
save('iML1515_ButOH.mat', 'model');


warning('on','all')

%% Preprossing
% clear;
close all;


% define target product
targetRxn='EX_ButOH';
load('iML1515_ButOH.mat')

% check model structure (essential components should be consistent in size)
assert(length(model.csense)==length(model.b) && length(model.b)==size(model.S,1), 'Number of mass balance constraints should be same as number of metabolites.');
assert(length(model.lb)==length(model.ub) && length(model.lb)== length(model.rxns) && length(model.lb)==size(model.S,2), 'number of reactions should be same as number of flux bounds.')

% name unspecified subsystems as 'Unassigned'
unassigned=cell2mat(cellfun(@(x) ~ischar(x), model.subSystems, 'UniformOutput',false));
model.subSystems(unassigned)={'Unassigned'};

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

% change carbon source
model=changeRxnBounds(model,{substrate,oxygenRxn}, [-20, 0], {'l', 'b'}); 

% knokout nanK/nanA
model = changeRxnBounds(model,{'AMANK','ACMNL'},[0, 0],['b','b']);
% model = changeRxnBounds(model,'ACNAMt2pp',0,'b');

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
    candidate=model.rxns; 
end

% make sure target reaction is in the reaction list of compressed model
if ~strcmp(model.rxns, targetRxn)
    newTargetRxn=model.rxns{contains(model.rxns, targetRxn)};
else
    newTargetRxn=targetRxn;
end

% remove unreasonable reactions from candidate knockout set
selectedRxns=setdiff(candidate, {'ATPM', biomassRxn, newTargetRxn});

disp(['Preprocessing step completed, ', num2str(length(selectedRxns)),' candidate reactions are selected for knockout.'])

%% -----------Set up algorithm--------------
% set options for OptDesign algorithm
regBlackList={substrate,'ATPM',biomassRxn, newTargetRxn, 'R03488' 'R03544','GRTT', 'crtE' 'EX_h2o','EX_o2','EX_h_e'};
options.regBlackList=struct('rxnList', {[regBlackList, regBlackList]}, 'typeReg','UUUUUUUUUUUDDDDDDDDDDD'); %'EX_h2o','EX_o2','EX_h_e','SUCCtex','SUCCt1pp','GLCt2pp'};
options.selectedRxns=selectedRxns(~contains(selectedRxns,'EX_')); % exclude transport reactions
options.targetRxn=newTargetRxn;
options.minGrowth=0.05; % minimal growth in production mutants
options.minProductRatio=0.99; % anticipated ratio of production
options.maxKO=10; % at most the number of knockouts
options.maxM=20; % maximal number of manipulations
options.changeThresh={'flux',0.5}; % noticeable flux change threshold
options.timeLimit=5000;   % a couple of mins


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

