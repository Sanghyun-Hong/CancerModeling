function [ essGeneIndexes, sol, optSol ] = checkESSinMModel( mmodel, method, essGenesNames )
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

isRecon1 = -1;
if isfield(mmodel, 'genes_unique_names')
  isRecon1 = 1;
elseif isfield(mmodel, 'genes_name')
  isRecon1 = 0;
else
  % exit the program
  fprintf('Usage: either RECON1 or RECON2 needs to be provided in the parameter\n');
  return
end


% Initialize variables
if isRecon1 == 1
  geneNames = mmodel.genes_unique_names;
elseif not(isRecon1)
  geneNames = mmodel.genes_name;
else
  % should not reach here
end

% Get an array of the index of essential genes
essGeneIndexes = find(ismember(geneNames, essGenesNames));

% Knockout a reaction associated to each reaction in the array
sol = zeros(length(essGeneIndexes), 1);
optSol = -1;
for i = 1:length(essGeneIndexes)
  modelC = mmodel;
  gi = essGeneIndexes(i);
  r = find(mmodel.rxnGeneMat(:,gi));
  modelC.ub(r) = 0;
  modelC.lb(r) = 0;
  if usingMOMA == 1
    % MOMA
    [output, outputWT] = MOMA(mmodel, modelC);
    sol(i) = output.f;
    if optSol == -1
      optSol = outputWT.f;
    end
  elseif not(usingMOMA)
    % FBA
    output = optimizeCbModel(modelC);
    sol(i) = output.f;
    if optSol == -1
      output = optimizeCbModel(mmodel);
      optSol = output.f;
    end
  else
    % should not reach here
  end
end

end
