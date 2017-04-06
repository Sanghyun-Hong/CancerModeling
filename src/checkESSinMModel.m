function [ ess ] = checkESSinMModel( mmodel, essGenesNames )
% Check the essentiality of a metabolic model with essential genes extracted from
% cacner model
% mmodel : metabolic model
% essGenesNames: a cell array of the name of essential genes from a cancer model
% ess: output --- essentiality

% Initialize variables
geneNames = model.genes_unique_names;
geneMap = model.genes_unique_map;
rxnMat = model.rxnGeneMat;
% Get an array of the index of essential genes
essGeneIndexes = find(ismember(geneNames, essGenesNames));
% TODO: Get an array of the index of associated reaction
% rxnIndexes = rxnMat(:, essGeneIndexes); % this code is not working

ess = None;

end
