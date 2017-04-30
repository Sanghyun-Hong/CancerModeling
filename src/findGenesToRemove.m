function [ genesToRem ] = findGenesToRemove( her2data, others, usingMedian )

% get the number of genes and cellines
[her2GeneNum, her2CellineNum] = size(her2data.ESS);
[otherGeneNum, othrCellineNum] = size(others.ESS);

if usingMedian == true
  % choose the median value of data expression in each cancer type
  her2Median = median(her2data.ESS, 2);
  otherMedian = median(others.ESS, 2);
  % extract genes to remoeve
  genesToRem = her2data.genes_ess(her2Median < otherMedian);
elseif usingMedian == false
  % calculate the average of data expression in each cancer type
  her2Expr = sum(her2data.ESS, 2) / her2CellineNum;
  otherExpr = sum(others.ESS, 2) / otherCellineNum;
  % extract genes to remoeve
  genesToRem = her2data.genes_ess(her2Expr < otherExpr);
else
  % Unknown case
  genesToRem = NaN;
end

end

