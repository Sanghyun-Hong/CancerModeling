% prefix direcotry for the datasets folder
datasets_dir = '';

% load ESSbrca.mat
load(strcat(datasets_dir, 'ESSbrca.mat'));
% load recon1.mat
load(strcat(datasets_dir, 'recon1.mat'));
% load recon2.mat
load(strcat(datasets_dir, 'recon2.mat'));

% initialize variables
cmodel = ESSbrca;
r1model = recon1;
r2model = recon2;
threshold = -0.5;
cellineRatio = 0.5;
method = 'FBA';

% get HER2 model
her2model = getHER2fromCModel(cmodel);

% get essential genes from HER2 model
essGenes = findESSGenesFromCModel(cmodel, threshold, cellineRatio);

% check if the genes are essential in a metabolic model
[r1essGeneIndexes, r1Sol, r1optSol] = checkESSinMModel(r1model, method, essGenes);
[r2essGeneIndexes, r2Sol, r2optSol] = checkESSinMModel(r2model, method, essGenes);
