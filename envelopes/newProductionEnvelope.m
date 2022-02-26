function [biomassValues, targetValues, lineHandle] = newProductionEnvelope(model, LP, lineColor, targetRxn, biomassRxn, nPts)
% Calculates the byproduct secretion envelope
%
% USAGE:
%
%    [biomassValues, targetValues, lineHandle] = productionEnvelope(model, deletions, lineColor, targetRxn, biomassRxn, geneDelFlag, nPts)
%
% INPUTS:
%    model:            COBRA model structure
%
% OPTIONAL INPUTS:
%    deletions:        List of reaction or gene deletions (empty if wild type)
%    lineColor:        Line color for plotting (see help plot for color codes)
%    targetRxn:        Target metabolite production reaction name
%    biomassRxn:       Biomass rxn name
%    geneDelFlag:      Perform gene and not reaction deletions
%    nPts:             Number of points in the plot
%
% OUTPUTS:
%    biomassValues:    Biomass values for plotting
%    targetValues:     Target upper and lower bounds for plotting
%    lineHandle:       Handle to lineseries object
%
% .. Author: - Markus Herrgard 8/28/06

if (nargin < 2)
    lineColor = 'k';
end
if (nargin < 3)
    % Target flux
    targetRxn = 'EX_etoh(e)';
end
if (nargin < 4)
    % Biomass flux
    biomassRxn = 'biomass_SC4_bal';
end
if (nargin < 5)
    nPts = 20;
end


% vector
vecBiomass=ismember(model.rxns,biomassRxn);
vecTarget=ismember(model.rxns,targetRxn);

% Run FBA to get upper bound for biomass
LP.c=[vecBiomass;vecBiomass];
solMax = optimizeCbModel(LP,'max');
solMin = optimizeCbModel(LP,'min');

% Create biomass range vector
biomassValues = linspace(solMin.f,solMax.f,nPts);

% Max/min for target production
LP.C(end+1,:)=-LP.c';
LP.csense(end+1)=char('L');
LP.d(end+1)=0;
LP.c=[vecTarget;vecTarget];
for i = 1:length(biomassValues)
    LP.d(end)=-biomassValues(i); % change growth value
    sol = optimizeCbModel(LP,'max');
    if (sol.stat > 0)
        targetUpperBound(i) = sol.f;
    else
        targetUpperBound(i) = NaN;
    end
    sol = optimizeCbModel(LP,'min');
    if (sol.stat > 0)
        targetLowerBound(i) = sol.f;
    else
        targetLowerBound(i) = NaN;
    end
end

% Plot results
lineHandle=plot([biomassValues fliplr(biomassValues)],[targetUpperBound fliplr(targetLowerBound)],lineColor,'LineWidth',2);
axis tight;
%ylabel([strrep(targetRxn,'_','-') ' (mmol/gDW h)']);
%xlabel('Growth rate (1/h)');

biomassValues = biomassValues';
targetValues = [targetLowerBound' targetUpperBound'];