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
cellineRatio = 0.9;
bmDropRatio = 0.9;
lexprCompMethod = 'median'; % median or average
exprDiffRatio = 0.01;

% name a csv file
deli = '_';
ext = '.csv';
csvfile = char(strcat(method, deli, string(threshold), deli, ...
             string(cellineRatio), deli, lexprCompMethod, deli, ...
             string(bmDropRatio), ext));
           
% open a csv file to write the method
fid = fopen(csvfile, 'w');
% first write parameters in the first row
fprintf(fid, '%s, %f, %f, %s, %f\n', method, threshold, cellineRatio, ...
  lexprCompMethod, bmDropRatio);

% get HER2 model
[her2data, others] = getHER2fromCancerData(cdata);

% get essential genes from HER2 model
essGenes = findESSGenesFromCancerData(cdata, threshold, cellineRatio);
% write a list of essential genes from HER cancer data into csv file
toCSV = essGenes';
fprintf(fid, '[ESSgenes of HER2],');
fprintf(fid, '%s,', toCSV{1:end-1});
fprintf(fid, '%s\n', toCSV{end});

% get essential genes of a metabolic model
essGenesOfModel = getESSGenesFromMModel(r1model, method, bmDropRatio);
% write a list of essential genes from default metabolic model into csv file
toCSV = essGenesOfModel';
fprintf(fid, '[ESSgenes of RECON1],');
fprintf(fid, '%s,', toCSV{1:end-1});
fprintf(fid, '%s\n', toCSV{end});

% TODO: first calculate the base accuracy (by Mr. Hong)

% choose the lower expression genes from HER2
% 3rd param:  using median = true
%             using average = false
genesToRem = findGenesToRemove(her2data, others, lexprCompMethod, exprDiffRatio);

% find rxns from genes
rxnsToRem = findRxnsWithGenesFromMModel(r1model, genesToRem);
rxnNames = mmodel.rxns(rxnsToRem);
% write a list of reactions to remove from metabolic model
toCSV = rxnNames';
fprintf(fid, '[Rxns to Remove],');
fprintf(fid, '%s,', toCSV{1:end-1});
fprintf(fid, '%s\n', toCSV{end});

for i = 1:length(rxnNames)
  mmodel = r1model;
  % modify a model by removing each reaction from the list
  mmodel = removeRxns(mmodel, rxnNames(i));
  
  % get essential genes of a MODIFIED metabolic model
  essGenesOfModel = getESSGenesFromMModel(mmodel, method, bmDropRatio);
  % write a list of essential genes from default metabolic model into csv file
  toCSV = essGenesOfModel';
  fprintf(fid, '[%s],', char(rxnNames(i)));
  fprintf(fid, '%s,', toCSV{1:end-1});
  fprintf(fid, '%s\n', toCSV{end});

  % TODO: calculate the accuracy of the modified model (by Mr. Hong)
  
end

fclose(fid);

