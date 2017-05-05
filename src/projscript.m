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
r1model = recon1;
%r2model = recon2;
method = 'FBA';
threshold = -0.5;
cellineRatio = 0.9;
bmDropRatio = 0.5;

% name a csv file
deli = '_';
ext = '.csv';
csvfile = char(strcat(method, deli, string(threshold), deli, ...
             string(cellineRatio), deli, string(bmDropRatio), ext));
           
% open a csv file to write the method
fid = fopen(csvfile, 'w');
% first write parameters in the first column
fprintf(fid, '%s, %f, %f, %f\n', method, threshold, cellineRatio, bmDropRatio);

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
genesToRem = findGenesToRemove(her2data, others, true);

% find rxns from genes
rxnsToRem = findRxnsWithGenesFromMModel(r1model, genesToRem);

for i = 1:length(rxnsToRem)
  mmodel = r1model;
  % modify a model by removing each reaction from the list
  rxnName = mmodel.rxns(rxnsToRem(i));
  mmodel = removeRxns(mmodel, rxnName);
  
  % get essential genes of a MODIFIED metabolic model
  essGenesOfModel = getESSGenesFromMModel(mmodel, method, bmDropRatio);
  % write a list of essential genes from default metabolic model into csv file
  toCSV = essGenesOfModel';
  fprintf(fid, '[%s],', char(rxnName));
  fprintf(fid, '%s,', toCSV{1:end-1});
  fprintf(fid, '%s\n', toCSV{end});

  % TODO: calculate the accuracy of the modified model (by Mr. Hong)
  
end

fclose(fid);

