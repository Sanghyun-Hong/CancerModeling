function [ foundRxns ] = findRxnsWithGenesFromMModel( mmodel, geneNames )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

geneIndexes = getGeneIndexFromMModel(mmodel, geneNames);

sol = zeros(length(geneIndexes), 1);
for i = 1:length(geneIndexes)
  gi = geneIndexes(i);
  foundRxns = [foundRxns; find(mmodel.rxnGeneMat(:,gi))];
end

unique(foundRxns);

end
