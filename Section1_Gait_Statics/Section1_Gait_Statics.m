function Section1_Gait_Statics()

    %  --- Setup GUI Window --- 
    Screensize = get(0, 'ScreenSize');
    
    figHeight = Screensize(4) * 0.6;
    figWidth = figHeight * 4/3; % Two axes side by side with 4:3 aspect ratio + padding

    fig = uifigure('Name', 'Footfall Viewer', 'Position',...
        [0.5*Screensize(3) - 0.5*figWidth, 0.5*Screensize(4) - 0.5*figHeight, figWidth, figHeight],...
        'Resize', 'off');
    
    buttonWidth = 0.10*Screensize(4);
    buttonHeight = 0.04*Screensize(4);

    currentfolder = pwd;


    %  --- Dropdown and Folder Button --- 

    folderButton = uibutton(fig, 'Text', 'Select Folder', ...
        'Position', [0.02 * figWidth, 0.91 * figHeight, buttonWidth, buttonHeight], ...
        'ButtonPushedFcn', @(src, event) SelectFolder());

    datasetDropdown = uidropdown(fig, ...
        'Position', [0.02 * figWidth + 1.2*buttonWidth, 0.91 * figHeight, 2*buttonWidth, buttonHeight], ...
        'Items', GetCSVFiles(), ...
        'ValueChangedFcn', @(src, event) LoadSelectedFile());

    ax = uiaxes(fig, 'Position', [0.4*figWidth, 0.2*figHeight, 0.8* figHeight, 0.6*figHeight], ...
        'Box', 'on');

    %  --- Initialize dropdown and load file if valid --- 
    initialFiles = GetCSVFiles();
    datasetDropdown.Items = initialFiles;
    if ~strcmp(initialFiles{1}, '<none>')
        datasetDropdown.Value = initialFiles{1};
        LoadSelectedFile();
    end

    % --- Nested Functions ---

    % Function to get .csv files
    function files = GetCSVFiles()
        csvFiles = dir(fullfile(currentfolder, '*.csv'));
        validFiles = {};
        for i = 1:length(csvFiles)
            try
                data = readtable(fullfile(currentfolder, csvFiles(i).name));
                if all(ismember({'LF', 'RF', 'LH', 'RH'}, data.Properties.VariableNames))
                    validFiles{end+1} = csvFiles(i).name;
                end
            catch
                continue;
            end
        end
        if isempty(validFiles)
            files = {'<none>'};
        else
            files = validFiles;
        end
    end

    % Folder selection and dropdown update
    function SelectFolder()
        folder = uigetdir(currentfolder);  % Start from previous folder
        if folder ~= 0
            currentfolder = folder;
            datasetDropdown.Items = GetCSVFiles();
            if ~strcmp(datasetDropdown.Items{1}, '<none>')
                datasetDropdown.Value = datasetDropdown.Items{1};
                LoadSelectedFile();
            else
                % Clear/reset UI state if no valid files
                % Add additional UI resets as needed
            end
            figure(fig);
        end
    end

    % Load selected CSV and extract footfall sequences
    function LoadSelectedFile()
        file = datasetDropdown.Value;
        individual_name = string(file(1:10)) + ' ' + string(file(12));
        if strcmp(file, '<none>'), return; end
        raw_data = readtable(fullfile(currentfolder, file));
        if ~all(ismember({'LF', 'RF', 'LH', 'RH'}, raw_data.Properties.VariableNames))
            uialert(fig, 'Selected file does not contain required footfall columns.', 'Error');
            return;
        end
        stridesequences = raw_data(:, {'LH', 'LF', 'RF', 'RH'});
        strideMatrix = table2array(stridesequences);    
        [stridesIds_midflight, gaitTypes_midflight] = StrideIndexExtractor(strideMatrix);
        
        % --- Call Visualization Function ---
         [allGaits, gaitCounts, transitionPercent]  = PlotGaitTransition(ax, gaitTypes_midflight, individual_name);
        % --- Build the table
        BuildAndShowGaitStatsTable(fig, allGaits, gaitCounts, transitionPercent, individual_name)
    end
end


%% Function to Plot Gait Transition Statics
function  [allGaits, gaitCounts, transitionPercent]  = PlotGaitTransition(ax, gaitTypes_midflight, individual_name)
    cla(ax);
    hold(ax, 'on');
    axis(ax, 'equal');
    axis(ax, 'off');
    box(ax, 'on');

    allGaits = ["TL", "TR", "RL", "RR"];
    allGaitColors = [
        8 163 119;
        8 118 179;
        8 163 119;
        8 163 179] / 256;

    positions = [
        -1,  1;
         1,  1;
        -1, -1;
         1, -1];

    fill(ax, [-2 2 2 -2], [-2 -2 0 0], [0.9 0.9 0.9], 'EdgeColor', 'none');
    fill(ax, [-2 2 2 -2], [0 0 2 2], [1 1 1], 'EdgeColor', 'none');

    numGaits = numel(allGaits);

    transitionPercent = zeros(numGaits);
    for i = 1:numGaits
        from = allGaits(i);
        idx_from = find(gaitTypes_midflight(1:end-1) == from);
        total_from = numel(idx_from);
        if total_from == 0, continue; end
        for j = 1:numGaits
            to = allGaits(j);
            count = sum(gaitTypes_midflight(idx_from + 1) == to);
            transitionPercent(i, j) = 100 * count / total_from;
        end
    end

    gaitCounts = zeros(1, numGaits);
    for i = 1:numGaits
        gaitCounts(i) = sum(gaitTypes_midflight == allGaits(i));
    end

    validIdx = gaitCounts > 0;
    gaitPercent = gaitCounts(validIdx) / sum(gaitCounts);
    gaitPercentFull = zeros(1, numGaits);
    gaitPercentFull(validIdx) = gaitPercent;

    % Dynamic scaling based on axes size
    axPos = ax.Position;
    baseDotSize = 0.001 * min(axPos(3:4));
    dotSizesFull = zeros(1, numGaits);
    dotSizesFull(validIdx) = baseDotSize * gaitPercent;

    baseLineWidth = 0.003 * min(axPos(3:4));
    arrowLineWidths = zeros(numGaits);
    for i = 1:numGaits
        fromTotal = sum(transitionPercent(i, :));
        if fromTotal > 0
            scale = 5 * baseLineWidth;
            arrowLineWidths(i, :) = baseLineWidth + scale * gaitPercentFull(i) .* sqrt(transitionPercent(i, :) / 100);
        end
    end

    % Parameters
    theta = linspace(0, 2*pi, 100);  % Circle resolution
    
    for i = 1:numGaits
        pos = positions(i, :);
        dotRadius = dotSizesFull(i);  % Directly in axis units now
    
        if validIdx(i)
            % Create circle patch centered at pos with radius dotRadius
            xCircle = pos(1) + dotRadius * cos(theta);
            yCircle = pos(2) + dotRadius * sin(theta);
    
            patch(ax, xCircle, yCircle, allGaitColors(i,:), ...
                'EdgeColor', 'k', 'LineWidth', baseLineWidth, ...
                'FaceAlpha', 1, 'FaceColor', allGaitColors(i,:));
        else
            % Draw an 'x' or empty circle for invalid entries
            line(ax, [pos(1) - dotRadius, pos(1) + dotRadius], ...
                     [pos(2) - dotRadius, pos(2) + dotRadius], ...
                     'Color', [0.5 0.5 0.5], 'LineWidth', baseLineWidth);
            line(ax, [pos(1) - dotRadius, pos(1) + dotRadius], ...
                     [pos(2) + dotRadius, pos(2) - dotRadius], ...
                     'Color', [0.5 0.5 0.5], 'LineWidth', baseLineWidth);
        end
    
        % Label
        labelOffsetScale = 0.6;  % adjust spacing factor if needed
        arrowTheta = atan2(pos(2), pos(1));
        textOffset = labelOffsetScale * (dotRadius + baseLineWidth + 0.001*ax.Position(3)) * [cos(arrowTheta), sin(arrowTheta)];

        text(ax, pos(1) + textOffset(1), pos(2) + textOffset(2), allGaits(i), ...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
    end

    for i = 1:numGaits
        for j = 1:numGaits
            perc = transitionPercent(i, j);
            if perc == 0, continue; end
            lw = arrowLineWidths(i, j);
            fromPos = positions(i, :);
            toPos = positions(j, :);

            if i == j
                % Ensure arc stays outside the dot plus padding
                dotRadius =  dotSizesFull(i);  % approximate radius from area
                padding =  0.02*lw;                 % small buffer beyond dot + outline
                r = dotRadius + padding;


                arrowTheta = atan2(fromPos(2), fromPos(1));
                thetaArc = linspace(1.25*pi + arrowTheta , 2.75*pi + arrowTheta ,100);
                labelOffsetDir = [fromPos(1) 0];
                xc = fromPos(1) + r * cos(thetaArc);
                yc = fromPos(2) + r * sin(thetaArc);
                plot(ax, xc, yc, 'Color', allGaitColors(i,:), 'LineWidth', lw);

                arrowEnd = [xc(end), yc(end)];
                dirVec = [xc(end) - xc(end-1), yc(end) - yc(end-1)];
                dirVec = dirVec / norm(dirVec);
                headLength = 0.04 * lw;
                headWidth  = 0.02 * lw;
                sideVec = headWidth * [-dirVec(2), dirVec(1)];
                tip   = arrowEnd + headLength * dirVec;
                left  = arrowEnd + sideVec;
                right = arrowEnd - sideVec;
                fill(ax, [tip(1), left(1), right(1)], [tip(2), left(2), right(2)], allGaitColors(i,:), 'EdgeColor', 'none');

                midIdx = round(length(xc)/2);
                arcLabelOffset = 0.0003*ax.Position(3)+ 0.03*dotRadius + 0.03 *lw;  % safe buffer outside arc
                labelPos = [xc(midIdx), yc(midIdx)] + arcLabelOffset * [cos(arrowTheta), sin(arrowTheta)];

                text(ax, labelPos(1), labelPos(2), sprintf('%.2f%%', perc), 'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', allGaitColors(i,:));
            else
                vec = toPos - fromPos;
                unitVec = vec / norm(vec);
                perpVec = [-unitVec(2), unitVec(1)];
                shiftAmount = 0.05;
                offset = shiftAmount * perpVec;
                startPos = fromPos + 0.15 * vec + offset;
                endPos = toPos - 0.15 * vec + offset;
                dirVec = endPos - startPos;
                dirVec = dirVec / norm(dirVec);
                headLength = 0.03 * lw;
                headWidth  = 0.015 * lw;
                sideVec = headWidth * [-dirVec(2), dirVec(1)];
                tip   = endPos;
                left  = endPos - headLength * dirVec + sideVec;
                right = endPos - headLength * dirVec - sideVec;
                fill(ax, [tip(1), left(1), right(1)], [tip(2), left(2), right(2)], allGaitColors(i,:), 'EdgeColor', 'none');
                trimmedEnd = endPos - headLength * dirVec;
                plot(ax, [startPos(1), trimmedEnd(1)], [startPos(2), trimmedEnd(2)], 'Color', allGaitColors(i,:), 'LineWidth', lw);

                midPos = (startPos + endPos) / 2;
                arrowDir = (toPos - fromPos);
                arrowDir = arrowDir / norm(arrowDir);
                perpVec = [-arrowDir(2), arrowDir(1)];
                offsetPerp = 0.12 * perpVec;
                offsetSource = 0.40 * (fromPos - midPos) / norm(fromPos - midPos);
                labelPos = midPos + offsetPerp + offsetSource;
                angleDeg = atan2d(arrowDir(2), arrowDir(1));
                text(ax, labelPos(1), labelPos(2), sprintf('%.2f%%', perc), 'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', allGaitColors(i,:), 'Rotation', angleDeg);
            end
        end
    end

    title(ax, {'Gait Transition Percentages', 'For ' + string(individual_name)}, 'FontSize', 0.035*ax.Position(3));
end

%% Function to build the table
function BuildAndShowGaitStatsTable(fig, allGaits, gaitCounts, transitionPercent, individual_name)
    numGaits = numel(allGaits);
    rowLabels = cellstr(allGaits);
    numRows = numGaits + 1;
    numCols = numGaits + 2;

    % Build data matrix
    data = cell(numRows, numCols);
    validIdx = gaitCounts > 0;

    % Column names: → transitions, then count
    columnNames = [{'Gait Type'}, cellstr("→" + rowLabels), {'Count'}];

    for i = 1:numGaits
        data{i, 1} = rowLabels{i};  % Gait name
        for j = 1:numGaits
            if ~validIdx(i) || transitionPercent(i, j) == 0
                data{i, j + 1} = '×';
            else
                data{i, j + 1} = sprintf('%d', round(transitionPercent(i, j) * gaitCounts(i) / 100));
            end
        end
        data{i, end} = gaitCounts(i);  % Count column
    end

    % Final row: total
    data{end, 1} = 'Total';
    for j = 1:numGaits
        totalIncoming = 0;
        for i = 1:numGaits
            if validIdx(i)
                totalIncoming = totalIncoming + round(transitionPercent(i, j) * gaitCounts(i) / 100);
            end
        end
        data{end, j + 1} = num2str(totalIncoming);
    end
    data{end, end} = sum(gaitCounts);

    % Layout
    figWidth = fig.Position(3);
    figHeight = fig.Position(4);
    xPos = 0.02 * figWidth;
    yPos = 0.35 * figHeight;
    tableWidth = 0.6 * figHeight;
    tableHeight = 0.3 * figHeight;
   

    columnWidthPixels = tableWidth / numCols;
    colWidths = num2cell(repmat(columnWidthPixels, 1, numCols));

    % Draw table (no FontSize specified)
    delete(findall(fig, 'Type', 'uitable'));
    uitable(fig, ...
        'Data', data, ...
        'Position', [xPos, yPos, tableWidth, tableHeight], ...
        'ColumnName', columnNames, ...
        'RowName', {}, ...
        'ColumnWidth', colWidths, ...
        'RowStriping', 'off');

    % Add title label
    delete(findall(fig, 'Type', 'uilabel'));
    table_title = 'Gait Transition Statistics\n for' + individual_name;
    uilabel(fig, ...
        'Text', sprintf(table_title), ...
        'Position', [xPos, yPos + tableHeight + 0.01*figHeight, tableWidth, 0.15*tableWidth], ...
        'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold', ...
        'FontSize', 0.045*tableWidth, ...
        'Tag', 'TableTitle');
end


%% stride index extraction
function [stridesIds_midflight, gaitTypes_midflight] = StrideIndexExtractor(stridesequences)
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
    gaitType = "RL";
    % gaitType = "RGL";
elseif isequal(order, [1 4 3 2])
    gaitType = "RR";
    % gaitType = "RGN";
elseif isequal(order, [4 1 3 2])
    gaitType = "TL";
elseif isequal(order, [1 4 2 3])
    gaitType = "TR";
else
    gaitType = "UNKNOWN";
end
end

