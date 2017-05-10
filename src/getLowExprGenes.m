function [ interGenes, pvals ] = getLowExprGenes( mmodel, her2, others )
% return p value of each gene

% find intersected genes between cancer data and metabolic model
[interGenes, cidx, midx] = intersect(her2.genes, mmodel.genes_unique_names);

X = her2.GEcell(cidx,:);
Y = others.GEcell(cidx,:);

pvals = zeros(size(interGenes));
for i = 1:length(interGenes)
  % calculate p value by using ranksum
  [pvals(i), h(i)] = ranksum(X(i,:), Y(i,:), 'tail', 'left');
end

end

