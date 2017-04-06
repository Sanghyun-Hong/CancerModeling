function [ essGenes ] = findESSGenesFromCModel( cmodel, threshold, clRatio )
% cmodel: Cacner model (for now, assuming ESSbrca.mat)
% threshold: threhsold to determine essential genes
%            generally less than -0.5 will be ESS genes
% clRatio: the ratio of essential ones through each gene's celllines to
% determine whether that gene is considered as essential gene
% essGenes: output --- a list of indexes of essential genes

[geneNum, cellLines] = size(cmodel.ESS);

essGenes = cmodel.genes_ess(sum((cmodel.ESS < threshold) == 1, 2) >= cellLines * clRatio);

end

