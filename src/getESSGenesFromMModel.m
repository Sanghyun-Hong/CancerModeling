function [ essGenesOfModel ] = getESSGenesFromMModel( mmodel, method, bmDropRatio )

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

for i = 1:length(mmodel.genes)
  modelC = mmodel;
  geneExp = sprintf('x(%d)', i);
  rxns = find(~cellfun(@isempty, strfind(mmodel.rules, geneExp)) > 0);
  
  x = ones(size(mmodel.genes));
  x(i) = 0;
  
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
    sol(i) = output.f;
  else
    FBAsol(i) = optimizeCbModel(modelC);
    sol(i) = FBAsol(i).f;
  end
end

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

essGenesIndexes = find(sol < optSol * (1 - bmDropRatio));
essGenesOfModel = cell(size(essGenesIndexes), 1);

for i = 1:length(essGenesIndexes)
  gi = essGenesIndexes(i);
  % store a name of each essential gene
  essGenesOfModel(i) = mmodel.genes_unique_names(mmodel.genes_unique_map(gi));
end

end

