function [ her2data, others ] = getHER2fromCancerData( cdata )
% Get a cacner model (most likely breat cancer model
% and then extract the subtype; HER2; from the given model

% cdata: cacner data of breast cancer
% her2model: output --- extracted from cdata

% a list of celllines associated to the HER2 subtype of breast cancer
% For more info, look at http://doi.org/10.1016/j.cell.2015.11.062
her2cl = {'au565'; 'bt474'; 'efm192a'; 'hcc1419'; 'hcc1569'; 'hcc1954'; ...
  'hcc202'; 'hcc2218'; 'mdamb330'; 'mdamb361'; 'ocubm'; 'skbr3'; 'sum190'; ...
  'sum225'; 'uacc812'; 'uacc893'; 'zr7530'};

her2data = cdata;
others = cdata;
% Extract from cdata a list of cellines related to HER2
relCellines = ismember(cdata.cellines, her2cl);
% Extract HER2 cancer data
her2data.cellines = cdata.cellines(relCellines);
her2data.ESS = cdata.ESS(:, relCellines);
% Extract other cancer data
others.cellines = cdata.cellines(not(relCellines));
others.ESS = cdata.ESS(:, not(relCellines));

end

