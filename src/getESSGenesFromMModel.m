function [ essGenes, notEssGenes ] = getESSGenesFromMModel( mmodel, method, foundGenesInds, bmDropRatio )

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

% get a base biomass value
if usingMOMA == 1
  [output, outputWT] = MOMA(mmodel, mmodel);
  optSol= outputWT.f;
else
  output = optimizeCbModel(mmodel);
  optSol = output.f;
end

essGenes = cell(length(foundGenesInds), 1);
notEssGenes = cell(length(foundGenesInds), 1);

for i = 1:length(foundGenesInds)
  modelC = mmodel;
  gi = foundGenesInds(i);
  geneExp = sprintf('x(%d)', gi);
  
  rxns = find(~cellfun(@isempty, strfind(mmodel.rules, geneExp)) > 0);
  x = ones(size(mmodel.genes));
  x(gi) = 0;
  
  for j = 1:length(rxns)
    rxnId = rxns(j);
    if(~eval(mmodel.rules{rxnId}))
      % knock-out
      modelC.ub(rxnId) = 0;
      modelC.lb(rxnId) = 0;
    end
  end
  
  if usingMOMA == 1
    [output, outputWT] = MOMA(mmodel, modelC);
    sol = output.f;
  else
    FBAsol = optimizeCbModel(modelC);
    sol = FBAsol.f;
  end
  
  if sol < optSol * (1 - bmDropRatio)
    essGenes(i) = mmodel.genes_unique_names(mmodel.genes_unique_map(gi));
  else
    notEssGenes(i) = mmodel.genes_unique_names(mmodel.genes_unique_map(gi));
  end
end

% compact cell arrays
essGenes = essGenes(~cellfun('isempty', essGenes));
notEssGenes = notEssGenes(~cellfun('isempty', notEssGenes));

% for i = 1:length(mmodel.genes)
%   modelC = mmodel;
%   inds = ones(size(mmodel.genes));
%   inds(i) = 0;
%   for j = 1:length(mmodel.rxns)
%     vv(j) = evalExpRule2(model.rules{j}, inds);
%   end
%   modelC.ub(vv==0) = 0;
%   modelC.lb(vv==0) = 0;
%     
%   if usingMOMA == 1
%     [output, outputWT] = MOMA(mmodel, modelC);
%     sol(i) = output.f;
%   else
%     FBAsol(i) = optimizeCbModel(modelC);
%     sol(i) = FBAsol(i).f;
%   end
% end

end

