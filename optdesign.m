function [solution] = optdesign(model, options, display)
% This is the core procedure of OptDesign consisting of two steps: step 1,
% identify potential regulation targets; step 2, find optimal combination of 
% knockout and regulation targets leading to the best lower bound of production 

% USAGE:
%
%    [solution, bilevelMILPProblem, gdlsSolutionStructs] = optdesign(model, varargin)
%
% INPUTS:
%    model:                  Cobra model structure
%    targetRxn:              Reaction(s) to be maximized (Cell array of strings)
%
% OPTIONAL INPUTS:           settings for OptDesign, default values will be
%                            used if not provided
%    varargin:               parameters entered using either a structure or list of
%                            parameter, parameter value. List of optional parameters:
%                              *  `maxM' - Maximum number of manipulations (default: 10)
%                              *  `maxKO` - Maximum number of knockouts (default: 5)
%                              *  `selectedRxns` - List of reactions that can be knocked out
%                              *  `regBlackList` - List of reactions that cannot be regulated (default: {})
%                              *  `targetRxn' - Target production reaction 
%                              *  `minProdRatio` - Minimum ratio of production to theoretical maxmimum value (default: 0.5)
%                              *  `minGrowth` - Minimum growth rate (default: 0.1)
%                              *  'changeThresold' - minimum noticeable flux change, which can be either absolute flux ({'flux',1}) 
%                                  or flux ratio ({'ratio',0.1}).   (default: {'flux',1})
%                              *  `timeLimit` - Maximum run time in seconds (default: 252000)
%
% OutputFlag:                display intermediate results and solver progress

%
% OUTPUTS:
%    solution:           solution structure (solution.pool(i) is the i-th solution solution returned)
%

%% Parse input paramets
MAXFLUX = 1000;
MAXDUAL = 1000;

if nargin<3
    display=0;
end

if (~exist('options','var') || isempty(options) )
    error('No target reaction specified')
else 
    %Default Values
    if isfield(options,'targetRxn')
        selTargetRxns = logical(ismember(model.rxns,options.targetRxn));
        if ~any(selTargetRxns)
            error([options.targetRxns ' not found. Double check spelling.']);
        end
    else
        error('No target reaction specified');
    end
    
    if ~isfield(options,'x0'), options.x0=[]; end
    if ~isfield(options,'selectedRxns'), options.selectedRxns=model.rxns;end
    if ~isfield(options,'maxKO'), options.maxKO=5; end
    if ~isfield(options,'maxM'), options.maxM=2; end
    if ~isfield(options,'minGrowth'), options.minGrowth=0; end
    if ~isfield(options,'minProdRatio'), options.minProdRatio=0.5; end
    if ~isfield(options,'regBlackList'), options.regBlackList={}; end
    if ~isfield(options,'changeThresh'), options.changeThresh={'ratio',0.01};end

    % import options for MILP solvers
    if ~isfield(options,'timeLimit'), options.timeLimit=252000; end
    if ~isfield(options,'Heuristics'), options.Heuristics = 1.0; end
    if ~isfield(options,'MIPFocus'), options.MIPFocus = 1; end
    if ~isfield(options,'ImproveStartGap'), options.ImproveStartGap = Inf; end
end

% solver options
solverOpt.timeLimit=options.timeLimit;
solverOpt.Heuristics=options.Heuristics;
solverOpt.MIPFocus=options.MIPFocus;
solverOpt.ImproveStartGap=options.ImproveStartGap;

if isfield(options,'selectedRxns')
    selSelectedRxns = logical(ismember(model.rxns,options.selectedRxns));
end

%% Generate selection reaction matrix
model.selRxnMatrix = selMatrix(selSelectedRxns)';
possibleKOList = model.rxns(selSelectedRxns);


objectiveRxn = model.c;

%% Setup model
model.ub(isinf(model.ub)) = MAXFLUX;
model.ub(model.ub>MAXFLUX) = MAXFLUX;
model.lb(isinf(model.lb)) = -MAXFLUX;
model.lb(model.lb<-MAXFLUX) = -MAXFLUX;
model.rxnsPresent = ones(length(model.rxns),1);

%% Create bi-level MILP problem
[nMets, nRxns] = size(model.S);

reg_up=ones(nRxns,1);
reg_down=reg_up;

% Use y indicating 
% vj=0,       if   yj=1
% lbj<vj<ubj  else yj=0


lb=model.lb;
ub=model.ub;

lbA=lb;
ubA=ub;
ubB=ub-lb;
lbB=-ubB;

ubB(~reg_up)=0;
lbB(~reg_down)=0;

% test model
model.lb=lb; model.ub=ub;
solWT=optimizeCbModel(model);
obj_c=model.c;
lbA(obj_c==1)=floor(solWT.f*1e6)*1e-6;

model.c=double(selTargetRxns);
solWP=optimizeCbModel(model);
lb(selTargetRxns)=options.minProductRatio*solWP.f; % minimal production rate

model.lb=lb;
model.c=obj_c;
solWT=optimizeCbModel(model);
lb(obj_c==1)=min(options.minGrowth,0.9*solWT.f); % minimal growth rate

% Calculate sigma: minimal noticeable flux change
if strcmp(options.changeThresh{1},'flux')
    sigma=options.changeThresh{2}*ones(nRxns,1);
elseif strcmp(options.changeThresh{1},'ratio')
    sigma=options.changeThresh{2}*ubB;
    sigma(sigma<0.1)=0.1; % assume flux change is unnoticeable if it is less than 0.1
else
    error('Noticeable flux threshold not specified correctly');
end

% objective
c = [zeros(2*nRxns,1); (1)*ones(2*nRxns,1)];


%% %%%%%%%%%% STEP 1: Identify a list of candidates for up/down regulation %%%%
%x= [ v                     u=detlav                             x+                                x-                                ]
A = [ 
    % Sv=0   //mass balance
    model.S                 sparse(nMets,nRxns)          sparse(nMets, 2*nRxns); 
    % Su=0   //mass balance
    sparse(nMets,nRxns)     model.S                      sparse(nMets, 2*nRxns);
    % lb<=u+v<=ub
    speye(nRxns)            speye(nRxns)                 sparse(nRxns, 2*nRxns);
    speye(nRxns)            speye(nRxns)                 sparse(nRxns, 2*nRxns);
    % u<=x+*ub+sigma
    sparse(nRxns,nRxns)     speye(nRxns)                -speye(nRxns).*(ubB)                  sparse(nRxns, nRxns);
    % u>=lb(1-x+)+sigma*x+
    sparse(nRxns,nRxns)     speye(nRxns)                 speye(nRxns).*(lbB-sigma)            sparse(nRxns, nRxns);
    % u>=lbx- -sigma
    sparse(nRxns,nRxns)     speye(nRxns)                 sparse(nRxns, nRxns)               -speye(nRxns).*(lbB)  ;
    % u<=ub(1-x-)-sigma*x-
    sparse(nRxns,nRxns)     speye(nRxns)                 sparse(nRxns, nRxns)                speye(nRxns).*(ubB+sigma);
    % v+u<=ub
    speye(nRxns)            speye(nRxns)                 sparse(nRxns, 2*nRxns)              ;
    % v+u>=lb
    speye(nRxns)            speye(nRxns)                 sparse(nRxns, 2*nRxns)              ;
    % x+ + x- <=1
    sparse(nRxns,nRxns)     sparse(nRxns,nRxns)          speye(nRxns)                        speye(nRxns) ;
    ];
b = [ zeros(nMets, 1);
    zeros(nMets, 1);
    ub;
    lb;
    sigma;
    lbB;
    -sigma;
    ubB;
    ub;
    lb;
    ones(nRxns, 1);
    ];
csense = char(['=' * ones(nMets, 1);
    '=' * ones(nMets, 1);
    '<' * ones(nRxns, 1);
    '>' * ones(nRxns, 1);
    '<' * ones(nRxns, 1);
    '>' * ones(nRxns, 1);
    '>' * ones(nRxns, 1);
    '<' * ones(nRxns, 1);
    '<' * ones(nRxns, 1);
    '>' * ones(nRxns, 1);
    '<' * ones(nRxns, 1);
    ]);
p_lb = [ lbA;  
    lbB;
    zeros(2*nRxns,1)];
p_ub = [ ubA;
    ubB;
    reg_up;
    reg_down;
    ];
vartype = char(['C' * ones(nRxns, 1);
    'C' * ones(nRxns, 1);
    'B' * ones(2*nRxns, 1); ]);
osense = 1; %minimize

[bilevelMILPProblem.c, bilevelMILPProblem.A,...
    bilevelMILPProblem.b, bilevelMILPProblem.lb,...
    bilevelMILPProblem.ub, bilevelMILPProblem.csense,...
    bilevelMILPProblem.vartype, bilevelMILPProblem.osense] = ...
    deal(c, A, b, p_lb, p_ub, csense, vartype, osense);

% Gurobi Solver settings
solOpt.OutputFlag=1;
solOpt.timeLimit=200;
%solOpt.MIPGap=0.05;
sol= solveCobraMILP(bilevelMILPProblem,solOpt);

if isempty(sol.full)
    error('no feasible solution can be found.')
end

%solution.x=sol.full;
y= sol.int> 1e-4;

clear solOpt A c b p_lb p_ub csense vartype osense


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% production requirement
lb(selTargetRxns)=0;
lb(obj_c==1)=options.minGrowth;

y_up=y(1:nRxns);
y_down=y(nRxns+1:2*nRxns);

% Set upper bound for regulation variables in black list
if isfield(options,'regBlackList')
    regBlackList=contains(model.rxns,options.regBlackList);
    y_up(regBlackList)=0;
    y_down(regBlackList)=0;
end

upCandidateRxns=model.rxns(y_up);
downCandidateRxns=model.rxns(y_down);

%% Get regulation candidates
if display
    v=sol.full(1:nRxns);
    u=sol.full(nRxns+1:2*nRxns);
    regUpDn=(y_up|y_down);
    up=abs(v(regUpDn)+u(regUpDn))>abs(v(regUpDn)); % absolute flux in mutant is larger than in wild type
    dn=abs(v(regUpDn)+u(regUpDn))<abs(v(regUpDn));
    possibleRegRxns=model.rxns(y_up|y_down);
    possibleRegRxns(up),possibleRegRxns(dn)
end

%% Create auxilary paramters
selUpMatrix = selMatrix(y_up)';
selDnMatrix = selMatrix(y_down)';
selKOMatrix = selMatrix(selSelectedRxns)';

nInt_u=size(selUpMatrix,2);
nInt_d=size(selDnMatrix,2);
nInt_x=sum(selSelectedRxns);

%% Create MILP 
% Objective
c=[-1*selTargetRxns; -1*selTargetRxns; zeros(2*nMets+6*nRxns+(nInt_u+nInt_d)+nInt_x,1);1e-5*ones(nInt_u+nInt_d+nInt_x,1)];
 
%[v, u, lmdv, lmdu a^, av, b^, bv, c^, cv, z+, z-, pi, y+, y-, y*]
B=[
    % Sv=0   //mass balance
    model.S                 sparse(nMets,nRxns)          sparse(nMets, 2*nMets+6*nRxns+2*(nInt_u+nInt_d+nInt_x));
    % Su=0   //mass balance
    sparse(nMets,nRxns)     model.S                      sparse(nMets, 2*nMets+6*nRxns+2*(nInt_u+nInt_d+nInt_x));
    % lb<=u+v<=ub
    speye(nRxns)            speye(nRxns)                 sparse(nRxns, 2*nMets+6*nRxns+2*(nInt_u+nInt_d+nInt_x));
    speye(nRxns)            speye(nRxns)                 sparse(nRxns, 2*nMets+6*nRxns+2*(nInt_u+nInt_d+nInt_x));
    % u>=lb(1-x+)+sigma*x+
    sparse(nRxns,nRxns)     speye(nRxns)                 sparse(nRxns, 2*nMets+6*nRxns+nInt_u+nInt_d+nInt_x)    (lbB-sigma).*selUpMatrix    sparse(nRxns, nInt_d+nInt_x);
    % u<=ub(1-x-)-sigma*x-
    sparse(nRxns,nRxns)     speye(nRxns)                 sparse(nRxns, 2*nMets+6*nRxns+nInt_u+nInt_d+nInt_x)    sparse(nRxns, nInt_u)     (sigma+ubB).*selDnMatrix  sparse(nRxns, nInt_x);
    % v+u<=ub(1-x^)
    speye(nRxns)            speye(nRxns)                 sparse(nRxns, 2*nMets+6*nRxns+2*(nInt_u+nInt_d)+nInt_x)     ub.*selKOMatrix;
    % v+u>=lb(1-x^)
    speye(nRxns)            speye(nRxns)                 sparse(nRxns, 2*nMets+6*nRxns+2*(nInt_u+nInt_d)+nInt_x)     lb.*selKOMatrix;
    % duality
    selTargetRxns'         selTargetRxns'               sparse(1, 2*nMets)     -lb'     ub'    -lbA'   ubA'  -lbB'  ubB' -(sigma-lbB)'*selUpMatrix   -(sigma+ubB)'*selDnMatrix sparse(1,nInt_x) sparse(1,nInt_u+nInt_d+nInt_x);
    % dual for v
    sparse(nRxns,2*nRxns)                                model.S'        sparse(nRxns,nMets)     speye(nRxns)    -speye(nRxns)    speye(nRxns)    -speye(nRxns) sparse(nRxns, (2*nRxns+nInt_u+nInt_d))  selKOMatrix   sparse(nRxns,nInt_u+nInt_d+nInt_x);
    % dual for u
    sparse(nRxns,2*nRxns)                                sparse(nRxns,nMets)      model.S'       speye(nRxns)    -speye(nRxns)    sparse(nRxns, 2*nRxns)   speye(nRxns)    -speye(nRxns)   sparse(nRxns, (nInt_u+nInt_d))  selKOMatrix   sparse(nRxns,nInt_u+nInt_d+nInt_x);
    % -H*y*<pi<Hy*
    sparse(nInt_x, 2*nRxns+2*nMets+6*nRxns+(nInt_u+nInt_d))  speye(nInt_x)   sparse(nInt_x,nInt_u+nInt_d) -speye(nInt_x);
    sparse(nInt_x, 2*nRxns+2*nMets+6*nRxns+(nInt_u+nInt_d))  speye(nInt_x)   sparse(nInt_x,nInt_u+nInt_d) speye(nInt_x);
    % z<c
     sparse(nInt_u+nInt_d, 2*nRxns+2*nMets+4*nRxns)  -sparse(blkdiag(selUpMatrix',selDnMatrix'))   speye(nInt_u+ nInt_d) sparse(nInt_u+ nInt_d,nInt_x+nInt_u+nInt_d+nInt_x);
    % z<Hy
    sparse(nInt_u+nInt_d, 2*nRxns+2*nMets+6*nRxns)   speye(nInt_u+ nInt_d)  sparse(nInt_u+nInt_d,nInt_x) -MAXDUAL*speye(nInt_u+ nInt_d) sparse(nInt_u+nInt_d,nInt_x);
    % z>c-H(1-y)
    sparse(nInt_u+nInt_d, 2*nRxns+2*nMets+4*nRxns)  -sparse(blkdiag(selUpMatrix',selDnMatrix'))   speye(nInt_u+ nInt_d)  sparse(nInt_u+nInt_d,nInt_x) -MAXDUAL*speye(nInt_u+ nInt_d) sparse(nInt_u+nInt_d,nInt_x);
    
    %%upper level constraints
    % y+  + y-   + y*<=1
    sparse(nRxns,2*nRxns+2*nMets+6*nRxns+(nInt_u+nInt_d+nInt_x))    selUpMatrix       selDnMatrix         selKOMatrix;
    % limit number of knockouts
    sparse(1,2*nRxns+2*nMets+6*nRxns+2*(nInt_u+nInt_d)+nInt_x)     ones(1, nInt_x);
    % limit total number of manipulations
    sparse(1,2*nRxns+2*nMets+6*nRxns+nInt_u+nInt_d+nInt_x)         ones(1, nInt_u+nInt_d+nInt_x);
    
    % target production constraint
    -1*selTargetRxns' -1*selTargetRxns'  sparse(1,2*nMets+6*nRxns+2*(nInt_u+nInt_d+nInt_x));
];

b = [ zeros(nMets, 1);
    zeros(nMets, 1);
    ub;
    lb;
    lbB;
    ubB;
    ub;
    lb;
    0;
    selTargetRxns;
    selTargetRxns;
    zeros(nInt_x, 1);
    zeros(nInt_x, 1);
    zeros(nInt_u+nInt_d, 1);
    zeros(nInt_u+nInt_d, 1);
    -MAXDUAL*ones(nInt_u+nInt_d, 1);
    
    ones(nRxns,1);
    options.maxKO;
    options.maxM;
    
    -0.1;
    ];
csense = char(['=' * ones(nMets, 1);
    '=' * ones(nMets, 1);
    '<' * ones(nRxns, 1);
    '>' * ones(nRxns, 1);
    '>' * ones(nRxns, 1);
    '<' * ones(nRxns, 1);
    '<' * ones(nRxns, 1);
    '>' * ones(nRxns, 1);
    '=';
    '=' * ones(nRxns, 1);
    '=' * ones(nRxns, 1);
    '<' * ones(nInt_x, 1);
    '>' * ones(nInt_x, 1);
    '<' * ones(nInt_u+nInt_d,1);
    '<' * ones(nInt_u+nInt_d,1);
    '>' * ones(nInt_u+nInt_d,1);

    '<' * ones(nRxns,1);
    '<';
    '<';
    
    '<';
    ]);
p_lb = [ lbA;  
    lbB;
    -Inf * ones(nMets, 1);
    -Inf * ones(nMets, 1);
    zeros(6*nRxns,1);
    zeros((nInt_u+nInt_d),1);
    -Inf * ones(nInt_x, 1);
    zeros(nInt_u+nInt_d,1);
    zeros(nInt_x,1);
    ];
p_ub = [ ubA;
    ubB;
    Inf * ones(nMets, 1);
    Inf * ones(nMets, 1);
    Inf * ones(6*nRxns,1);
    Inf * ones((nInt_u+nInt_d),1);
    Inf * ones(nInt_x, 1);
    ones(nInt_u+nInt_d,1);
    ones(nInt_x,1);
    ];
vartype = char(['C' * ones(2*nRxns, 1);
    'C' * ones(2*nMets+6*nRxns+(nInt_u+nInt_d+nInt_x), 1);
    'B' * ones(nInt_u+nInt_d+nInt_x, 1); 
    ]);

osense = 1; %minimize

[bilevelMILPProblem.c, bilevelMILPProblem.A,...
    bilevelMILPProblem.b, bilevelMILPProblem.lb,...
    bilevelMILPProblem.ub, bilevelMILPProblem.csense,...
    bilevelMILPProblem.vartype, bilevelMILPProblem.osense,...
    bilevelMILPProblem.x0] = ...
    deal(c, B, b, p_lb, p_ub, csense, vartype, osense, options.x0);

solverOpt.OutputFlag=1;
%solverOpt.MIPFocus=3;
solverOpt.Heuristics=1;
solverOpt.PoolSolutions=10;
solverOpt.PoolSearchMode=2;
sol= solveCobraMILP(bilevelMILPProblem,solverOpt);

clear selUpMatrix selDnMatrix selKOMatrix bilevelMILPProblem solverOpt 

if isempty(sol.full)
    error('no feasible solution can be found.')
end

lb(objectiveRxn==1)=0;
biomassRxn=model.rxns{objectiveRxn==1};

if ~isfield(sol,'pool')|| isempty(fieldnames(sol.pool)) % if gurobi 7.5 or older verision used
    sol.pool(1).xn=sol.full;
end


for i=1:length(sol.pool)
    % iterate all pool solutions and get values for integer variables
    y=sol.pool(i).xn(end+1-(nInt_u+nInt_d+nInt_x):end)>1e-4;
    [yu, yd, yx]=deal(y(1:nInt_u),y(1+nInt_u:nInt_u+nInt_d), y(1+nInt_u+nInt_d:end));
    
    % get flux
    v=sol.pool(i).xn(1:nRxns);
    u=sol.pool(i).xn(nRxns+1:2*nRxns);    
    
    %Generate Solution Structure
    fprintf('\nGenerating Output\n');
    
    %solution.y=y;
    solution.pool(i).KOset = possibleKOList(logical(yx));
    
    regUpDnRxns=union(upCandidateRxns(logical(yu)),downCandidateRxns(logical(yd)));
    [~,regUpDn]=ismember(regUpDnRxns,model.rxns);
    
    up=abs(v(regUpDn)+u(regUpDn))>abs(v(regUpDn)); % absolute flux in mutant is larger than in wild type
    dn=abs(v(regUpDn)+u(regUpDn))<abs(v(regUpDn));
    
    solution.pool(i).UPset=regUpDnRxns(up);
    solution.pool(i).DOWNset=regUpDnRxns(dn);
    solution.pool(i).wt_v=v;
    solution.pool(i).mt_v=u+v;
end
% can use as a start solution in other runs
solution.x=sol.full;

clear sol y yu yd yx

end
