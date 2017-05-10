clear;

clRatio = 0.9;
essThreshold = -0.5;
bmDropRatio = 0.1;
pValThreshold = 0.01;

runProjectScript(clRatio, essThreshold, bmDropRatio, pValThreshold);

% clRatio = [0.9];
% essThreshold = -0.5;
% dropRatio = [0.1];
% pValThreshold = [0.05];
% 
% for i = 1:length(clRatio)
%   for j = 1:length(bmDropRatio)
%     for k = 1:length(pvalThreshold)
%       runProjectScript(clRatio(i), essThreshold, bmDropRatio(j), pValThreshold(k));
%     end
%   end
% end

