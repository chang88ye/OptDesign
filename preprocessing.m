function [redModel,candidate]=preprocessing(model,substrateRxn,oxygenRxn,biomassRxn,targetRxn)

% collect a priori knowledge about undoable/essential genes/reactions/subsystems.
options=blackList(model);
gene_essential_data=readtable('ecoli_essential_genes.csv');
options.expEssentialGenes=gene_essential_data{:,2};

options.substrateRxn=substrateRxn;
options.oxygenRxn=oxygenRxn;
options.biomassRxn=biomassRxn;
options.targetRxn=targetRxn;
options.loadRedModel=true;

[redModel,candidate]=selectCandidates(model,options);
redModel.substrateRxns=findSubstrateInReducedModel(model,substrateRxn);
end


%% subfunctions
function [reducedModel,candidate]=selectCandidates(model,options)
tol=1e-4;

if ~exist('reduced', 'dir')
       mkdir('reduced');
end
%% step0: compress model,including removing zeroflux rxns and lumping rxns
redModelFile=['reduced', filesep, 'reduced_', model.description '.mat'];
if isfield(options, 'loadRedModel') && options.loadRedModel && exist(redModelFile,'file')
    load(redModelFile);
else
    disp('reduction process starts:')
    [reducedModel,biomassRxn,targetRxn] = compressor(model,options);
    reducedModel.biomassRxn=biomassRxn;
    reducedModel.targetRxn=targetRxn;
    disp('reduction process ends!')
    
    % test if reduced model is consistent with original model
    WTGsol=optimizeCbModel(model);
    RWTGsoL = optimizeCbModel(reducedModel,'max');
    if (abs(RWTGsoL.f-WTGsol.f)/WTGsol.f>5e-2)
        error('the reduced model is inconsistent with the original model!')
    end
    
    % save reduced model
    save(redModelFile,'reducedModel');
end

selectedRxns=reducedModel.rxns;

%% step1: removal of rxns acting on molecules containing more than 7/10 carbons
if (length(model.rxns)>500) % avoid reduction for small models.
    deleteRxns = findCarbonRxns(reducedModel,200);
    selectedRxns=setdiff(selectedRxns,deleteRxns);
end

%% step2: removal of non-gene associated rxns
nonGeneAssocRxnID=cellfun(@isempty,reducedModel.rules);
deleteRxns=reducedModel.rxns(nonGeneAssocRxnID);
selectedRxns=setdiff(selectedRxns,deleteRxns);

%% step3: removel of rxns related essential genes (both in vivo and in silico genes)
essentialGenes=union(reducedModel.essentialGenes,options.expEssentialGenes);
deleteRxns=findEssentialRxnsFromGenes(reducedModel,essentialGenes);
selectedRxns=setdiff(selectedRxns,deleteRxns);

%% step4: exclusion of centain subsystems
tmpSubSystems=cellfun(@(x) strsplit(x,'/'), reducedModel.subSystems, 'UniformOutput',false);
isSubset = cellfun(@(a)all(ismember(a,options.selSubSystems)),tmpSubSystems);
deleteRxns = reducedModel.rxns(isSubset);
selectedRxns=setdiff(selectedRxns,deleteRxns);

%% step5: removal of transport reactions with non genes
transRxns= findTransRxns(reducedModel,true); % Remove transport reactions
[~,transIds]=ismember(transRxns,reducedModel.rxns);
nonGeneIds=cellfun(@isempty,reducedModel.rules(transIds));
selectedRxns=setdiff(selectedRxns,transRxns(nonGeneIds));

%% remove exchange reactions
% transRxns = findTransRxns(reducedModel,true);
% selectedRxns=setdiff(selectedRxns,transRxns);
% excRxns=selectedRxns(contains(selectedRxns, 'EX_'));
% selectedRxns=setdiff(selectedRxns,excRxns);

%% step6: removal of special reactions: glycogen, thioredoxin & flavodoxin
tmpRxns=cellfun(@(x) strsplit(x,'/'), selectedRxns, 'UniformOutput',false);
isSubset = cellfun(@(x)all(ismember(x,options.excludedRxns)),tmpRxns);
selectedRxns=selectedRxns(~isSubset);

%% step7: removal of computationally essential reactions
% fva with 5% growth rate of wild type.
% [minflux,maxflux]= fluxVariability(reducedModel,5);
minflux=reducedModel.mins;maxflux=reducedModel.maxs;
zeroFluxRxns=reducedModel.rxns(abs(minflux)<tol & abs(maxflux)<tol);
mustFluxRxns=reducedModel.rxns(maxflux<-tol | minflux>tol);
deleteRxns=union(zeroFluxRxns,mustFluxRxns);
selectedRxns=setdiff(selectedRxns,deleteRxns);

%% Change target rxn bounds
reducedModel=changeRxnBounds(reducedModel,{options.biomassRxn,options.targetRxn},[0 0],{'l','l'});


%% get genes from reactions
selectedGenes=findGenesFromRxns(reducedModel,selectedRxns);
selectedGenes = cat(1,selectedGenes{:});
selectedGenes=setdiff(selectedGenes,reducedModel.essentialGenes);

%% get geneSets
isSelected=ismember(reducedModel.rxns,selectedRxns);
getGeneSets=any(reducedModel.geneSetRxnMat(:,isSelected),2);
selectedGeneSets=setdiff(reducedModel.geneSets(getGeneSets),reducedModel.geneSets(reducedModel.essentialGeneSets));
selectedGeneSets=setdiff(selectedGeneSets,reducedModel.essentialGenes);

%% return selected candidates
candidate=struct();
candidate.rxns=selectedRxns;
candidate.genes=selectedGenes;
candidate.geneSets=selectedGeneSets;
end

function options=blackList(model)
% return options struct: [expEssentialGenes, excludedRxns, excludedSubSystems] 

%  special reactions: glycogen, thioredoxin & flavodoxin
specialRxns=...
    {'GLCP','GLCP2','GLCS1','GLCS2',... % glycogen
    'FTR','TRDR','TRDRm','ASR2','THIORDXi','THIOGabc',... %thioredoxin
    'FLNDPR2r','FLDR2','FLDR','PFOR'};

% Remove some reactions from allowable knockouts
switch model.description
    case {'Ecoli_core_model','e_coli_core','iAF1260','iAF1260b', 'iJO1366', 'iML1515'}
        preserveRxns=...
            {'ENO','GAPD','PGK','PGM',...       % Essential reactions in actual cell.
            'H2Ot','H2Otex','EX_h2o_e',...
            'EX_acald_e','ACALDtex','ACALDt',...
            'EX_co2_e','CO2tex','CO2t','EX_h_e','Htex'... % Hard to knockout only these reactions.
            'CAT','DHPTDNR','DHPTDNRN',...
            'FHL', 'SPODM', 'SPODMpp',...
            'SUCASPtpp','SUCFUMtpp',...
            'SUCMALtpp','SUCTARTtpp'};         % thermodynamically infeasible reactions
    otherwise
        error('model description incorrect.')
end

options=struct();

options.excludedRxns = union(preserveRxns,specialRxns);
options.selSubSystems = {'Cell Envelope Biosynthesis', 'Glycerophospholipid Metabolism','Inorganic Ion Transport and Metabolism','Lipopolysaccharide Biosynthesis / Recycling',...
    'Membrane Lipid Metabolism','Murein Biosynthesis', 'Murein Recycling','Transport Inner Membrane', 'Transport, Inner Membrane','Transport Outer Membrane', 'Transport, Outer Membrane', ...
    'Transport Outer Membrane Porin','tRNA Charging', 'Nucleotide Salvage Pathway','Unassigned',%''
    };
end

function substrateRxns=findSubstrateInReducedModel(model, subtrate)
if ~iscell(subtrate)
    substrateRxns={};
    for j=1:length(model.rxns)
        tmp=strfind(model.rxns{j},subtrate);
        if ~isempty(tmp)
            substrateRxns=model.rxns(j);
            return;
        end
    end
else
    substrateRxns={};
    for i=1:length(subtrate)
        for j=1:length(model.rxns)
            tmp=strfind(model.rxns{j},subtrate);
            if ~isempty(tmp)
                substrateRxns= [substrateRxns,model.rxns(j)];
                break;
            end
        end
    end
end
end

function essentialRxns=findEssentialRxnsFromGenes(model,essentialGenes)
[~,geneIds]=ismember(essentialGenes,model.genes);
geneIds=geneIds(geneIds>0);
if isempty(geneIds)
    warning('target genes not found in the model');
    essentialRxns={};
end

essentialRxns={};
x0 = true(size(model.genes));
for i=1:length(geneIds)
    x=x0; x(geneIds(i))=false;
    rxnIds=find(model.rxnGeneMat(:,geneIds(i)));
    for j=1:length(rxnIds)
        
        if (~isempty(model.rules{rxnIds(j)})) %To avoid errors if the rule is empty
            if (~eval(model.rules{rxnIds(j)}))
                essentialRxns{end+1}=model.rxns{rxnIds(j)};
            end
        end
    end
end
end