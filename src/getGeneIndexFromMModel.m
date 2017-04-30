function [ geneIndexes ] = getGeneIndexFromMModel( mmodel, geneNames )
% check whether a given mmodel is RECON1 or RECON2
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
  modelGeneNames = mmodel.genes_unique_names;
elseif not(isRecon1)
  modelGeneNames = mmodel.genes_name;
else
  % should not reach here
end

% Get an array of the index of essential genes
geneIndexes = find(ismember(modelGeneNames, geneNames));


end

