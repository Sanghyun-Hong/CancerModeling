function runProjectScript( cellineRatio, essThreshold, bmDropRatio, pValThreshold )

% prefix direcotry for the datasets folder
datasets_dir = '';

% load ESSbrca.mat
load(strcat(datasets_dir, 'ESSbrca.mat'));
% load recon1.mat
load(strcat(datasets_dir, 'recon1.mat'));
% load recon2.mat
load(strcat(datasets_dir, 'recon2.mat'));

% initialize variables
cdata = ESSbrca;
r1model = defineHumanMediaNCI60(recon1, '');
%r2model = recon2;
method = 'FBA'; % FBA or MOMA
%lexprCompMethod = exprCompMethod; % median or average
%exprDiffRatio = 0.01;

% name a csv file
deli = '_';
ext = '.csv';
csvfile = char(strcat(method, deli, string(essThreshold), deli, ...
             string(cellineRatio), deli, string(bmDropRatio), deli, ...
             pValThreshold, ext));
           
% open a csv file to write the method
fid = fopen(csvfile, 'w');
% first write parameters in the first row
fprintf(fid, 'method, %s\n', method);
fprintf(fid, 'ess threshold, %f\n', essThreshold);
fprintf(fid, 'celline ratio, %f\n', cellineRatio);
fprintf(fid, 'biomass drop ratio, %f\n', bmDropRatio);
fprintf(fid, 'p value threshold, %f\n', pValThreshold);

% get HER2 model
[her2data, others] = getHER2fromCancerData(cdata);

% get essential genes from HER2 model
essGenes = findESSGenesFromCancerData(cdata, essThreshold, cellineRatio);
% write a list of essential genes from HER cancer data into csv file
toCSV = essGenes';
fprintf(fid, '<HER2 ESSgenes>,');
fprintf(fid, '%s,', toCSV{1:end-1});
fprintf(fid, '%s\n', toCSV{end});

% get a list of essGenes found in RECON1
foundGeneInds = find(ismember(r1model.genes_unique_names, essGenes));
foundGeneNames = r1model.genes_unique_names(foundGeneInds);
% write a list of essential genes from HER cancer data into csv file
toCSV = foundGeneNames';
fprintf(fid, '<HER2 ESSgenes found in RECON1>,');
fprintf(fid, '%s,', toCSV{1:end-1});
fprintf(fid, '%s\n', toCSV{end});

% get essential genes of a metabolic model
essGeneNames = findESSGenesFromMModel(r1model, bmDropRatio);
fprintf(fid, '[DEFAULT RECON1]\n');
% write a list of essential genes found in default metabolic model into csv file
toCSV = essGeneNames';
fprintf(fid, '<Essential>,');
if ~isempty(toCSV)
  fprintf(fid, '%s,', toCSV{1:end-1});
  fprintf(fid, '%s\n', toCSV{end});
else
  fprintf(fid, '\n');
end

% choose the low expression genes from HER2 --- deprecated 5/9/2017
%genesToRem = findGenesToRemove(her2data, others, lexprCompMethod, exprDiffRatio);
% Instead, we use p value to choose low expression gene data
[interGenes, pvals] = getLowExprGenes(r1model, her2data, others);
genesToRem = interGenes(pvals < pValThreshold);

% find rxns from genes
rxnsToRem = findRxnsWithGenesFromMModel(r1model, genesToRem, essGeneNames);
rxnNames = r1model.rxns(rxnsToRem);
% write a list of reactions to remove from metabolic model
toCSV = rxnNames';
fprintf(fid, 'Rxns to Remove,');
fprintf(fid, '%s,', toCSV{1:end-1});
fprintf(fid, '%s\n', toCSV{end});

for i = 1:length(rxnNames)
  mmodel = r1model;
  % modify a model by removing each reaction from the list
  mmodel = removeRxns(mmodel, rxnNames(i));
  
  % write the name of a removed reaction
  fprintf(fid, '[%s]\n', char(rxnNames(i)));
  % get essential genes of a MODIFIED metabolic model
  essGeneNames = findESSGenesFromMModel(mmodel, bmDropRatio);
  
  % write a list of essential genes from default metabolic model into csv file
  toCSV = essGeneNames';
  fprintf(fid, '<Essential>,');
  if ~isempty(toCSV)
    fprintf(fid, '%s,', toCSV{1:end-1});
    fprintf(fid, '%s\n', toCSV{end});
  else
    fprintf(fid, '\n');
  end
end

fclose(fid);

end