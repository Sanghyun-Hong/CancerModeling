function [ foundRxns ] = findRxnsWithGenesFromMModel( mmodel, geneNames, essGeneNames )
% return a list of reactions that neither blocked nor essential

geneIndexes = getGeneIndexFromMModel(mmodel, geneNames);

foundRxns = [];
for i = 1:length(geneIndexes)
  gi = geneIndexes(i);
  foundRxns = [foundRxns; find(mmodel.rxnGeneMat(:,gi))];
end

unique(foundRxns);

% find blocked reactions using FVA
er = cellfun('length', strfind(mmodel.rxns,'EX_'));
modelC = mmodel;
modelC.ub(er==1) = 1000;
modelC.lb(er==1) = -1000;

[minFlux, maxFlux, Vmin, Vmax] = fluxVariability(modelC);
blockedRxns = find((abs(minFlux) < 1e-10) & (abs(maxFlux) < 1e-10));

% get reactions that are not blocked
foundRxns = foundRxns(~ismember(foundRxns, blockedRxns)); 

% now find essential genes
essGeneIndexes = getGeneIndexFromMModel(mmodel, essGeneNames);
% find reactions related to the essential genes
foundEssRxns = [];
for i = 1:length(essGeneIndexes)
  gi = essGeneIndexes(i);
  foundEssRxns = [foundEssRxns; find(mmodel.rxnGeneMat(:,gi))];
end

unique(foundEssRxns);

% get reactions that are not essential
foundRxns = foundRxns(~ismember(foundRxns, foundEssRxns));

end
