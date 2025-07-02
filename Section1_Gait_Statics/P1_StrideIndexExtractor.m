function [stridesIds, stridesIds_midflight, gaitTypes, gaitTypes_midflight] = P1_StrideIndexExtractor(stridesequences)
% P1_StrideIndexExtractor detects strides and gait types from 4-limb stance/flight data
% Input:  stridesequences (N x 4), each row is a time frame; columns are [LH, LF, RF, RH]
% Output:
%   stridesIds          (3 x M)     = [start of flight; end of flight; end of stride]
%   stridesIds_midflight (2 x M-1)  = [startIdx; endIdx] using midpoints between flights
%   gaitTypes           (1 x M)     = string array: gait classification for each stride
%   gaitTypes_midflight (1 x M-1)   = gait type labels corresponding to midflight strides

% Transpose input to 4 x N format for compatibility
stridesequences = stridesequences.';  

[numLegs, numSamples] = size(stridesequences);
stridesIds = [];
gaitTypes = [];

idx = 2;
inFlightPhase = false;
flightStartIdx = -1;
flightEndIdx = -1;
maxFlightLen = 0;
maxFlightStartIdx = -1;
maxFlightEndIdx = -1;

legCompleted = zeros(1, numLegs);
trackingStride = false;

while idx <= numSamples
    prevState = stridesequences(:, idx - 1);
    currState = stridesequences(:, idx);

    % Detect start of flight phase
    if any(prevState == 1) && all(currState == 0)
        inFlightPhase = true;
        flightStartIdx = idx;
    end

    % Detect end of flight phase
    if all(prevState == 0) && any(currState == 1) && inFlightPhase
        flightEndIdx = idx - 1;
        inFlightPhase = false;

        currFlightLen = flightEndIdx - flightStartIdx + 1;
        if currFlightLen > maxFlightLen
            maxFlightLen = currFlightLen;
            maxFlightStartIdx = flightStartIdx;
            maxFlightEndIdx = flightEndIdx;

            trackingStride = true;
            legCompleted = zeros(1, numLegs);
        end
    end

    % Track touchdown/liftoff for each leg
    if trackingStride && maxFlightStartIdx > 0
        for leg = 1:numLegs
            seg = stridesequences(leg, maxFlightStartIdx:idx);
            if any(diff(seg) == 1) && any(diff(seg) == -1)
                legCompleted(leg) = 1;
            end
        end

        if all(legCompleted)
            strideEnd = -1;
            tempIdx = idx;
            while tempIdx < numSamples
                if any(stridesequences(:, tempIdx - 1) == 1) && all(stridesequences(:, tempIdx) == 0)
                    strideEnd = tempIdx - 1;
                    break;
                end
                tempIdx = tempIdx + 1;
            end

            if strideEnd > maxFlightStartIdx
                stridesIds = [stridesIds, [maxFlightStartIdx; maxFlightEndIdx; strideEnd]];
                strideSegment = stridesequences(:, maxFlightStartIdx:strideEnd);
                gait = ClassifyGaitFromTouchdown(strideSegment);
                gaitTypes = [gaitTypes, gait];

                % Reset
                idx = strideEnd + 1;
                inFlightPhase = false;
                flightStartIdx = -1;
                flightEndIdx = -1;
                maxFlightLen = 0;
                maxFlightStartIdx = -1;
                maxFlightEndIdx = -1;
                legCompleted = zeros(1, numLegs);
                trackingStride = false;
                continue;
            else
                break;
            end
        end
    end

    idx = idx + 1;
end

gaitTypes = string(gaitTypes);  % Convert to string array
gaitTypes_midflight = gaitTypes(1:end-1);
stridesIds_midflight = RefineStrideBoundaries(stridesIds);
gaitTypes_midflight = string(gaitTypes_midflight);

end

%% Refine stride boundaries using flight midpoints
function midflightstrides = RefineStrideBoundaries(strides)
numStrides = size(strides, 2);
if numStrides < 2
    midflightstrides = [];
    return;
end

midflightstrides = zeros(2, numStrides - 1);
for i = 1:numStrides - 1
    startMid = round((strides(1, i) + strides(2, i)) / 2);
    endMid   = round((strides(1, i + 1) + strides(2, i + 1)) / 2) - 1;
    midflightstrides(:, i) = [startMid; endMid];
end
end

%% Gait classification from touchdown timing
function gaitType = ClassifyGaitFromTouchdown(strideSegment)
% Leg order: [LH; LF; RF; RH] = [1 2 3 4]
numFrames = size(strideSegment, 2);
touchdownFrame = NaN(1, 4);

for leg = 1:4
    for t = 2:numFrames
        if strideSegment(leg, t - 1) == 0 && strideSegment(leg, t) == 1
            touchdownFrame(leg) = t;
            break;
        end
    end
end

[~, order] = sort(touchdownFrame);

if isequal(order, [4 1 2 3])
    gaitType = "RGL";
    % gaitType = "RGL";
elseif isequal(order, [1 4 3 2])
    gaitType = "RGR";
    % gaitType = "RGN";
elseif isequal(order, [4 1 3 2])
    gaitType = "TGL";
elseif isequal(order, [1 4 2 3])
    gaitType = "TGR";
else
    gaitType = "UNKNOWN";
end
end
