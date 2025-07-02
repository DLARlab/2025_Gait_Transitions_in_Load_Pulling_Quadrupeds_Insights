function [avgForces,peakForces] = P1_LoadingForceStatistics(leash_force, strides_midflight)
% ComputeLeashForcePerStride computes mean leash force per stride
% Inputs:
%   leash_force: N x 1 vector of leash force values
%   strides_midflight: 2 x M matrix with stride start/end indices (per column)
% Output:
%   avgForces: 1 x M vector of average leash force for each stride

numStrides = size(strides_midflight, 2);
avgForces = zeros(1, numStrides);
peakForces = zeros(1, numStrides);

for i = 1:numStrides
    startIdx = strides_midflight(1, i);
    endIdx = strides_midflight(2, i);
    
    if startIdx >= 1 && endIdx <= length(leash_force) && endIdx >= startIdx
        segment = leash_force(startIdx:endIdx);
        avgForces(i) = mean(segment);
        peakForces(i) = max(segment);
    else
        avgForces(i) = NaN; % if indices are invalid
    end
end
end
