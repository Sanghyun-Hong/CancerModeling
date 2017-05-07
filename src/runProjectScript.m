function runProjectScript( clRatio, dropRatio, exprCompMethod )

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
threshold = -0.5; 
cellineRatio = clRatio;
bmDropRatio = dropRatio;
lexprCompMethod = exprCompMethod; % median or average
exprDiffRatio = 0.01;

% name a csv file
deli = '_';
ext = '.csv';
csvfile = char(strcat(method, deli, string(threshold), deli, ...
             string(cellineRatio), deli, lexprCompMethod, deli, ...
             string(bmDropRatio), deli, string(exprDiffRatio), ext));
           
% open a csv file to write the method
fid = fopen(csvfile, 'w');
% first write parameters in the first row
fprintf(fid, 'method, %s\n', method);
fprintf(fid, 'threshold, %f\n', threshold);
fprintf(fid, 'celline ratio, %f\n', cellineRatio);
fprintf(fid, 'geneExpr compare method, %s\n', lexprCompMethod);
fprintf(fid, 'biomass drop ratio, %f\n', bmDropRatio);
fprintf(fid, 'gene express diff ratio, %f\n', exprDiffRatio);

% get HER2 model
[her2data, others] = getHER2fromCancerData(cdata);

% get essential genes from HER2 model
essGenes = findESSGenesFromCancerData(cdata, threshold, cellineRatio);
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

% choose the low expression genes from HER2
genesToRem = findGenesToRemove(her2data, others, lexprCompMethod, exprDiffRatio);

% find rxns from genes
rxnsToRem = findRxnsWithGenesFromMModel(r1model, genesToRem);
rxnNames = r1model.rxns(rxnsToRem);
% write a list of reactions to remove from metabolic model
toCSV = rxnNames';
fprintf(fid, 'Rxns to Remove,');
fprintf(fid, '%s,', toCSV{1:end-1});
fprintf(fid, '%s\n', toCSV{end});

% get essential genes of a metabolic model
[essGenesOfModel, nonESSGenesOfModel] = getESSGenesFromMModel(r1model, method, foundGeneInds, bmDropRatio);
fprintf(fid, '[DEFAULT RECON1]\n');
% write a list of essential genes found in default metabolic model into csv file
toCSV = essGenesOfModel';
fprintf(fid, '<Essential>,');
if ~isempty(toCSV)
  fprintf(fid, '%s,', toCSV{1:end-1});
  fprintf(fid, '%s\n', toCSV{end});
else
  fprintf(fid, '\n');
end
% write a list of non-essential genes found in default metabolic model into csv file
toCSV = nonESSGenesOfModel';
fprintf(fid, '<Non-Essential>,');
if ~isempty(toCSV)
  fprintf(fid, '%s,', toCSV{1:end-1});
  fprintf(fid, '%s\n', toCSV{end});
else
  fprintf(fid, '\n');
end

% TODO: first calculate the base accuracy (by Mr. Hong)

for i = 1:length(rxnNames)
  mmodel = r1model;
  % modify a model by removing each reaction from the list
  mmodel = removeRxns(mmodel, rxnNames(i));
  % write the name of a removed reaction
  fprintf(fid, '[%s]\n', char(rxnNames(i)));
  % get essential genes of a MODIFIED metabolic model
  [essGenesOfModel, nonESSGenesOfModel] = getESSGenesFromMModel(mmodel, method, foundGeneInds, bmDropRatio);
  % write a list of essential genes from default metabolic model into csv file
  toCSV = essGenesOfModel';
  fprintf(fid, '<Essential>,');
  if ~isempty(toCSV)
    fprintf(fid, '%s,', toCSV{1:end-1});
    fprintf(fid, '%s\n', toCSV{end});
  else
    fprintf(fid, '\n');
  end
  % write a list of non-essential genes found in default metabolic model into csv file
  toCSV = nonESSGenesOfModel';
  fprintf(fid, '<Non-Essential>,');
  if ~isempty(toCSV)
    fprintf(fid, '%s,', toCSV{1:end-1});
    fprintf(fid, '%s\n', toCSV{end});
  else
    fprintf(fid, '\n');
  end

  % TODO: calculate the accuracy of the modified model (by Mr. Hong)
  
end

% remove all reactions from metabolic model
mmodel = r1model;
for i = 1:length(rxnNames)
  % modify a model by removing each reaction from the list
  mmodel = removeRxns(mmodel, rxnNames(i));
end

fprintf(fid, '[All reactions removed]\n');
% get essential genes of a MODIFIED metabolic model
[essGenesOfModel, nonESSGenesOfModel] = getESSGenesFromMModel(mmodel, method, foundGeneInds, bmDropRatio);
% write a list of essential genes from default metabolic model into csv file
toCSV = essGenesOfModel';
fprintf(fid, '<Essential>,');
if ~isempty(toCSV)
  fprintf(fid, '%s,', toCSV{1:end-1});
  fprintf(fid, '%s\n', toCSV{end});
else
  fprintf(fid, '\n');
end
% write a list of non-essential genes found in default metabolic model into csv file
toCSV = nonESSGenesOfModel';
fprintf(fid, '<Non-Essential>,');
if ~isempty(toCSV)
  fprintf(fid, '%s,', toCSV{1:end-1});
  fprintf(fid, '%s\n', toCSV{end});
else
  fprintf(fid, '\n');
end

% TODO: calculate the accuracy of the modified model (by Mr. Hong)

fclose(fid);

end