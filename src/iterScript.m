clear;

clRatio = [0.7 0.6 0.5];
dropRatio = [0.05 0.1 0.3];
exprCompMethod = {'median', 'average'};

for i = 1:length(clRatio)
  for j = 1:length(exprCompMethod)
    for k = 1:length(dropRatio)      
      runProjectScript(clRatio(i), dropRatio(k), char(exprCompMethod(j)));
    end
  end
end

