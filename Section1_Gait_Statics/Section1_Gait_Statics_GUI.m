function Section1_Gait_Statics_GUI
    %% Setup GUI Window
    Screensize = get(0, 'ScreenSize');

    figHeight = Screensize(4) * 0.6;
    figWidth = figHeight * 4/3; % Two axes side by side with 4:3 aspect ratio + padding

    fig = uifigure('Name', 'Footfall Viewer', 'Position', [0.5*Screensize(3) - 0.5*figWidth, 0.5*Screensize(4) - 0.5*figHeight, figWidth, figHeight]);
    
    buttonWidth = 0.10*Screensize(4);
    buttonHeight = 0.04*Screensize(4);

    currentfolder = pwd;

    %% Dropdown and Folder Button

    folderButton = uibutton(fig, 'Text', 'Select Folder', ...
        'Position', [0.02 * figWidth, 0.91 * figHeight, buttonWidth, buttonHeight], ...
        'ButtonPushedFcn', @(src, event) SelectFolder());

    datasetDropdown = uidropdown(fig, ...
        'Position', [0.02 * figWidth + 1.2*buttonWidth, 0.91 * figHeight, 2*buttonWidth, buttonHeight], ...
        'Items', {});

    %% Axes (Blank for Now)

%     ax1 = uiaxes(fig, 'Position', [0.05*figWidth, 0.05*figHeight, 0.4* figHeight, 0.3* figHeight], ...
%         'Box', 'on');
    ax = uiaxes(fig, 'Position', [0.4*figWidth, 0.2*figHeight, 0.8* figHeight, 0.6*figHeight], ...
        'Box', 'on');

    %% Initialize dropdown and load file if valid
    initialFiles = GetCSVFiles();
    datasetDropdown.Items = initialFiles;
    if ~strcmp(initialFiles{1}, '<none>')
        datasetDropdown.Value = initialFiles{1};
        LoadSelectedFile();
    end

    %% --- Nested Functions ---

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
        if strcmp(file, '<none>'), return; end
        raw_data = readtable(fullfile(currentfolder, file));
        if ~all(ismember({'LF', 'RF', 'LH', 'RH'}, raw_data.Properties.VariableNames))
            uialert(fig, 'Selected file does not contain required footfall columns.', 'Error');
            return;
        end
        stridesequences = ExtractFootfall(raw_data);
        strideMatrix = table2array(stridesequences);    
        [stridesIds_midflight, gaitTypes_midflight] = S1_StrideIndexExtractor(strideMatrix);
        
        % --- Call Visualization Function ---
         [allGaits, gaitCounts, transitionPercent]  = PlotGaitTransition(ax, gaitTypes_midflight);
        % --- Build the table
        BuildAndShowGaitStatsTable(fig, ax.Position, allGaits, gaitCounts, transitionPercent)
    end
end

%% Function to build the table
function BuildAndShowGaitStatsTable(fig, axPos, allGaits, gaitCounts, transitionPercent)
    numGaits = numel(allGaits);
    rowLabels = cellstr(allGaits);
    columnNames = [{'Gait Type'}, strcat('-->', rowLabels)];

    validIdx = gaitCounts > 0;

    % Build table data
    data = cell(numGaits, numGaits + 1);  % 4 rows × (1 label + 4 transitions)
    for i = 1:numGaits
        data{i, 1} = rowLabels{i};
        if ~validIdx(i)
            data(i, 2:end) = {'×'};
        else
            for j = 1:numGaits
                val = transitionPercent(i, j);
                if gaitCounts(i) == 0 || val == 0
                    data{i, j + 1} = '×';
                else
                    data{i, j + 1} = sprintf('%d', round(val * gaitCounts(i) / 100));
                end
            end
        end
    end

    % Position table left of the axes
    tableWidth = 0.35 * axPos(3);  % relative to figure
    tableHeight = axPos(4);
    tableLeft = axPos(1) - tableWidth - 20;  % small gap
    tableBottom = axPos(2);

    % Remove previous tables (if any)
    delete(findall(fig, 'Type', 'uitable'));

    % Create new table
    uitable(fig, ...
        'Data', data, ...
        'Position', [tableLeft, tableBottom, tableWidth, tableHeight], ...
        'ColumnName', [{'Gait Type'}, strcat('-->', rowLabels).'], ...
        'RowName', {}, ...
        'FontSize', 12, ...
        'ColumnWidth', num2cell([80, repmat(70, 1, numGaits)]));
end

%% Function to extract footfall sequences 
function stridesequences = ExtractFootfall(raw_data)
    stridesequences = raw_data(:, {'LH', 'LF', 'RF', 'RH'});
    % Display or use footfall_data as needed
end

%% Function to Plot Gait Transition Statics
function  [allGaits, gaitCounts, transitionPercent]  = PlotGaitTransition(ax, gaitTypes_midflight)
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
        labelOffsetScale = 0.8;  % adjust spacing factor if needed
        arrowTheta = atan2(pos(2), pos(1));
        textOffset = labelOffsetScale * (dotRadius + baseLineWidth) * [cos(arrowTheta), sin(arrowTheta)];

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
                arcLabelOffset = 0.005*dotRadius + 0.1 *lw;  % safe buffer outside arc
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

    title(ax, {'Gait Transition Percentages', 'For Alice'}, 'FontSize', 14);
end

