warning('off','all')

%% Preprossing
% clear;
% close all;
hold on

[parentdir,~,~]=fileparts(pwd);

% define target product
% targetRxn='EX_succ_e';
targetRxn='EX_lyco_e';
% targetRxn='EX_narg_e';

if strcmp(targetRxn,'EX_lyco_e')
    load(fullfile(parentdir,'iML1515_lyco.mat'))
end
 
if strcmp(targetRxn,'EX_narg_e')
    load(fullfile(parentdir,'iML1515_narg.mat'))
end

if strcmp(targetRxn,'EX_succ_e')
    sel_model = input('Which model do you want to use, Ecoli_core_model (1) or iML1515 (2) ?');
    if sel_model==1
        % use E. coli core model
        load(fullfile(parentdir,'ecoli_core_model.mat'));     
    elseif sel_model==2
        load(fullfile(parentdir,'iML1515.mat'))
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

% Attention: make sure you use reaction names in compressed models.
regList={{},''}; % wild type

% regList ={{'PGI','PTAr/ACKr','GND','ATPS4rpp', 'SUCDi', 'HEX1','MALS','PPC'}, 'kkkkkuuu'}; % succ
% regList ={{'TALA','ALLTN/ALLTAMH2/UGCIAMH','DHAPT','PFK', 'SUCDi', 'FRD3','MALS','PGM'},'kkkkkuud'}; % succ-ref

% regList = {{'CPPPGO','OHPBAT/PERD/E4PD','ALCD19','TKT2','AI2K/PAI2I/PAI2T','EX_pi_e/PItex','TPI'}, 'kkkkkdu'}; % lyco
% regList = {{'CPPPGO';'ALCD19';'AI2K/PAI2I/PAI2T';'EX_pi_e/PItex';'TPI'}, 'kkkdu'}; % lyco-ref
% regList ={{'ALCD19';'AI2K/PAI2I/PAI2T';'EX_pi_e/PItex';'TPI'}, 'kkdu'};
regList={{'ALCD19';'FDH4pp';'AI2K/PAI2I/PAI2T'; 'DXPRIi/CDPMEK/MEPCT/MECDPS/MECDPDH5';'EX_pi_e/PItex';'TPI'}, 'kkkudu'};

% regList ={{'IMPD';'TRPAS2';'ASNS2';'ATPS4rpp';'AI2K/PAI2I/PAI2T';'EX_pi_e/PItex';'H2Otpp';'POR5';'FUM';'RPE'}, 'kkkkkddudd'}; % narg
% regList ={{'EDA/EDD';'GHMT2r';'ATPS4rpp';'GLUDy';'CYTBO3_4pp'; 'EX_pi_e/PItex';'AKGDH';'PGCD';'PPA';'TKT1'}, 'kkkkkddddd'};% narg-ref
% regList ={{'PPM';'ASNS2';'CBMKr';'ATPS4rpp';'AI2K/PAI2I/PAI2T';'EX_pi_e/PItex';'RPI'; 'ENO'}, 'kkkkkddd'}; 
%  regList ={{'IMPD';'TRPAS2';'ASNS2';'ATPS4rpp';'AI2K/PAI2I/PAI2T';'EX_pi_e/PItex';'H2Otpp';'PGI';'SUCDi'}, 'kkkkkddud'};


% regList ={{'ALCD19','TKT2','AI2K','PItex', 'TPI'}, 'kkkdu'};
% regList ={{'R1PK','FUM', 'ADK3', 'PItex', 'ADK1'}, 'kkddd'};

[minGrowthRate, delta] =deal(0.1, 1);

LP=buildLPFromStrategies(model, regList, newTargetRxn, minGrowthRate, delta);


newProductionEnvelope(model,LP, 'b',newTargetRxn,biomassRxn,20);

