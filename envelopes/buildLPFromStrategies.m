function LP=buildLPFromStrategies(model, regList,targetRxn, minGrowthRate, delta)

% regList must be in the form of {{'FUM','ATPS4rpp','DXPS'}, 'kku'};

% wild typle flux bounds
[w_lb, w_ub] =deal(model.lb, model.ub);
solWT=optimizeCbModel(model);
w_lb(model.c==1)=floor(solWT.f*1e6)*1e-6;

% mutant flux bounds
[m_lb, m_ub]=deal(model.lb, model.ub);
m_lb(model.c==1)=minGrowthRate;

[~, nRxns]=size(model.S);

koIDs=strfind(regList{2}, 'k');
ko_vector=zeros(1, nRxns);
%iscontains=cellfun(@(c)contains(c,regList{1}(koIDs)),model.rxns);
iscontains=findRxnIDs(model,regList{1}(koIDs));
ko_vector(iscontains)=1;

% Delta-v boudns
d_lb=-Inf*ones(nRxns,1);
d_ub= Inf*ones(nRxns,1);


for i=1:length(regList{2})
    if regList{2}(i) == 'u'
%         d_lb(contains(model.rxns, regList{1}{i}))=delta;
        d_lb(findRxnIDs(model, regList{1}{i}))=delta;
    elseif regList{2}(i) == 'd'
        d_ub(findRxnIDs(model, regList{1}{i}))=-delta;
    end
end

% set objective
c = zeros(size(model.c));
c(findRxnIDs(model, targetRxn))=1;
obj=[c;c];

% create LP constraints: Ax=b, Cx<=d
A=[blkdiag(model.S, model.S); % sv=0, su=0
    selMatrix(ko_vector), selMatrix(ko_vector)]; %v+u=0 for knockouts

b=zeros(size(A,1),1); %b=0

C=[-eye(nRxns),-eye(nRxns); %-(v+u)<=-lb
    eye(nRxns), eye(nRxns)];

d=[-model.lb; model.ub];

csense=char(['E'*ones(size(A,1),1); 'L'*ones(size(d))]);

lb=[w_lb; d_lb]; 
ub=[w_ub; d_ub];

[LP.c, LP.S, LP.dxdt, LP.C, LP.d, LP.csense, LP.lb, LP.ub] = deal(obj, A, b, C, d, csense, lb, ub);


end