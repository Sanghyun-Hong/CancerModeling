function [ essGenes ] = findESSGenesFromCancerData( cdata, threshold, clRatio )
% cmodel: Cacner model (for now, assuming ESSbrca.mat)
% threshold: threhsold to determine essential genes
%            generally less than -0.5 will be ESS genes
% clRatio: the ratio of essential ones through each gene's celllines to
% determine whether that gene is considered as essential gene
% essGenes: output --- a list of the name of essential genes

[geneNum, cellLines] = size(cdata.ESS);

essGenes = cdata.genes_ess(sum((cdata.ESS < threshold) == 1, 2) >= cellLines * clRatio);

end
