clear;

clRatio = 0.9;
essThreshold = -0.5;
bmDropRatio = 0.1;
pValThreshold = 0.01;

runProjectScript(clRatio, essThreshold, bmDropRatio, pValThreshold);

% clRatio = [0.9, 0.8];
% dropRatio = [0.1, 0.05];
% pValThreshold = [0.05, 0.01, 0.005, 0.001];

% for i = 1:length(clRatio)
%   for j = 1:length(dropRatio)
%     for k = 1:length(pvalThreshold)
%       runProjectScript(clRatio(i), dropRatio(j), pValThreshold(k));
%     end
%   end
% end
% 
