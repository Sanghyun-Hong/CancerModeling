function [ her2model ] = getHER2fromCModel( cmodel )
% Get a cacner model (most likely breat cancer model
% and then extract the subtype; HER2; from the given model

% cmodel: cacner model of breast cancer
% her2model: output --- extracted from cmodel

% a list of celllines associated to the HER2 subtype of breast cancer
% For more info, look at http://doi.org/10.1016/j.cell.2015.11.062
her2cl = {'au565'; 'bt474'; 'efm192a'; 'hcc1419'; 'hcc1569'; 'hcc1954'; ...
  'hcc202'; 'hcc2218'; 'mdamb330'; 'mdamb361'; 'ocubm'; 'skbr3'; 'sum190'; ...
  'sum225'; 'uacc812'; 'uacc893'; 'zr7530'};

her2model = cmodel;
% Extract from cmodel a list of cellines related to HER2
relCellines = ismember(cmodel.cellines, her2cl);
% Create HER2 subtype model
her2model.celllines = cmodel.cellines(relCellines);
her2model.ESS = cmodel.ESS(:, relCellines);

end

