% prefix direcotry for the datasets folder
datasets_dir = '/Users/jun/workspace/CancerModeling/datasets/';

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
threshold = -0.5;
cellineRatio = 0.9;
method = 'FBA';

% get HER2 model
[her2data, others] = getHER2fromCancerData(cdata);

% get essential genes from HER2 model
essGenes = findESSGenesFromCancerData(cdata, threshold, cellineRatio);

% check if the genes are essential in a metabolic model
%[r1essGeneIndexes, r1Sol, r1optSol] = checkESSinMModel(r1model, method, essGenes);
%[r2essGeneIndexes, r2Sol, r2optSol] = checkESSinMModel(r2model, method, essGenes);

% TODO: first calculate the base accuracy (by Mr. Hong)


% choose the lower expression genes from HER2
% 3rd param:  using median = true
%             using average = false
genesToRem = findGenesToRemove(her2data, others, true);

% find rxns from genes
rxnsToRem = findRxnsWithGenesFromMModel(r1model, genesToRem);

% TODO: modify a model by removing each reaction from the list

% TODO: calculate the accuracy of the modified model (by Mr. Hong)
