function [ genesToRem ] = findGenesToRemove( her2data, others, lexprCompMethod, diffRatio )

if strcmp(lexprCompMethod, 'median')
  % choose the median value of data expression in each cancer type
  her2Expr = median(her2data.ESS, 2);
  otherExpr = median(others.ESS, 2);
elseif strcmp(lexprCompMethod, 'average')
  % get the number of genes and cellines
  [her2GeneNum, her2CellineNum] = size(her2data.ESS);
  [otherGeneNum, otherCellineNum] = size(others.ESS);
  % calculate the average of data expression in each cancer type
  her2Expr = sum(her2data.ESS, 2) / her2CellineNum;
  otherExpr = sum(others.ESS, 2) / otherCellineNum;
else
  % Unknown case
  genesToRem = NaN;
  return
end

% first compare the positive value 
her2Positive = her2Expr;
otherPositive = otherExpr;
her2Positive(her2Expr < 0) = 0;
otherPositive(otherExpr < 0) = 0;
pExprGenes = her2data.genes_ess(her2Positive < otherPositive * diffRatio);
% then compare the negative value 
her2Negative = her2Expr;
otherNegative = otherExpr;
her2Negative(her2Expr > 0) = 0;
otherNegative(otherExpr > 0) = 0;
nExprGenes = her2data.genes_ess(her2Negative * diffRatio < otherNegative);
% extract genes to remoeve
genesToRem = union(pExprGenes, nExprGenes);

end

