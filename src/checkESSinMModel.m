function [ essGeneIndexes, sol ] = checkESSinMModel( mmodel, method, essGenesNames )
% Check the essentiality of a metabolic model with essential genes extracted from
% cacner model
% mmodel: metabolic model
% method: either FBA or MOMA
% essGenesNames: a cell array of the name of essential genes from a cancer model
% ess: output --- essentiality

usingMOMA = -1;
if strcmp(method, 'MOMA')
  usingMOMA = 1;
elseif strcmp(method, 'FBA')
  usingMOMA = 0;
else
  % exit the program
  fprintf('Usage: either MOMA or FBA needs to be specified in the parameter\n');
  return
end

% Initialize variables
geneNames = model.genes_unique_names;
% Get an array of the index of essential genes
essGeneIndexes = find(ismember(geneNames, essGenesNames));

% Knockout a reaction associated to each reaction in the array
sol = zeros(length(essGeneIndexes));
output = zeros(length(essGeneIndexes));
for i = 1:length(essGeneIndexes)
  modelC = mmodel;
  gi = essGeneIndexes(i);
  r = find(mmodel.rxnGeneMat(:,gi));
  modelC.ub(r) = 0;
  modelC.lb(r) = 0;
  if usingMOMA == 1
    % MOMA
    output(i) = MOMA(mmodel, modelC);
    sol(i) = output(i).f;
  elseif not(usingMOMA)
    % FBA
    output(i) = optimizeCbModel(modelC);
    sol(i) = output(i).f;
  else
    % should not reach here
  end
end

end
