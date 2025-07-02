%% Data Extraction

%% Extract the data from Excel
clear
clc
% select horizontal position first
[filename_data, pathname_data] = uigetfile('*.csv','Select a File');
raw_data = P1_ImportFile(filename_data);

% extract dog name
dog_name_index = strfind(filename_data,'_');
dog_name       = filename_data(1:dog_name_index-1);



% extract timestep
frequency      = str2num(filename_data(dog_name_index(3)+1:dog_name_index(3)+3)); 
dt             = 1/frequency; 
T              = 0:dt:dt*(size(raw_data,1)-1);


% extract footfall sequence
LH_touchdown_sequence = raw_data(:,find(strcmp(raw_data.Properties.VariableNames,'LH'))).LH;
LF_touchdown_sequence = raw_data(:,find(strcmp(raw_data.Properties.VariableNames,'LF'))).LF;
RF_touchdown_sequence = raw_data(:,find(strcmp(raw_data.Properties.VariableNames,'RF'))).RF;
RH_touchdown_sequence = raw_data(:,find(strcmp(raw_data.Properties.VariableNames,'RH'))).RH;

% % extract loading force
% loading_force           = raw_data(:,find(strcmp(raw_data.Properties.VariableNames,'force'))).force;

% extrac load with calibration
load     = raw_data(:,find(strcmp(raw_data.Properties.VariableNames,'load'))).load;
loading_force = 0.00259*(load + 6027.006); % for Alice
% loading_force = 0.00262*(load-86579.096); % for Chips

% extract horizontal velocity and pitching velocity

a_horizontal          = raw_data(:,find(strcmp(raw_data.Properties.VariableNames,'ay'))).ay;
a_vertical            = raw_data(:,find(strcmp(raw_data.Properties.VariableNames,'az'))).az;
o_pitching            = raw_data(:,find(strcmp(raw_data.Properties.VariableNames,'gx'))).gx;

% normalize data:       ay,az units: g          gx unit: degree/s

% normalize angular velocity data
o_pitching = deg2rad(o_pitching);
%% Data Normalization

m_Alice = 30; % Mass of Alice: 30kg
m_Chips = 25; % Mass of Chips: 25kg

ll_Alice = 0.5; % Leg length of Alice: 50cm
ll_Chips = 0.48; % Leg length of Chips: 48cm

g = 9.81; % gravity constant

%% Gait Index and Types Buffer
stridesequences = [LH_touchdown_sequence LF_touchdown_sequence RF_touchdown_sequence RH_touchdown_sequence];
[stridesIds, stridesIds_midflight, gaitTypes, gaitTypes_midflight] = P1_StrideIndexExtractor(stridesequences);
%% Plot Selected Strides

close all
clc

% index_start = 627;
% intercal = 47;

% stirde_index  = 5;
% index_start = strides(1,stirde_index);
% intercal = strides(3,stirde_index) - strides(1,stirde_index);

% stirde_index  = 5;
% index_start = strides_midflight(1,stirde_index);
% intercal = strides_midflight(2,stirde_index) - strides_midflight(1,stirde_index);

stirde_index  = 19;
index_start = stridesIds_midflight(1,stirde_index);
intercal = stridesIds_midflight(2,stirde_index+1) - stridesIds_midflight(1,stirde_index);


index = index_start:intercal+index_start;

[LH_TD_XData,LH_TD_YData] = P1_patch_indices_TD(LH_touchdown_sequence(index));
[LF_TD_XData,LF_TD_YData] = P1_patch_indices_TD(LF_touchdown_sequence(index));
[RF_TD_XData,RF_TD_YData] = P1_patch_indices_TD(RF_touchdown_sequence(index));
[RH_TD_XData,RH_TD_YData] = P1_patch_indices_TD(RH_touchdown_sequence(index));

%     td_indicator = clip(LH_touchdown_sequence(index)+LF_touchdown_sequence(index)+RF_touchdown_sequence(index)+RH_touchdown_sequence(index), 0, 1);
td_indicator = (LH_touchdown_sequence(index)+LF_touchdown_sequence(index)+RF_touchdown_sequence(index)+RH_touchdown_sequence(index))~=0;



figure
subplot(3,1,1)
hold on
pbaspect([10 1 1])

xlim([0 intercal])

ylim([0 4.5])
yticks([1 2 3 4])
yticklabels({'RH', 'RF', 'LF', 'LH'})
patch(LH_TD_XData,LH_TD_YData+3,'black')
patch(LF_TD_XData,LF_TD_YData+2,'black')
patch(RF_TD_XData,RF_TD_YData+1,'black')
patch(RH_TD_XData,RH_TD_YData,'black')
title('Footfall Pattern')

subplot(3,1,2)
plot(index,loading_force(index))
xlim([index_start intercal+index_start])
title('loading Force')

subplot(3,1,3)
plot(index,td_indicator)
xlim([index_start intercal+index_start])
ylim([0 2])

%% Plot Strides
intercal = 1000;

for i = 1:1

    index = intercal*(i-1)+1:intercal*i;
    [LH_TD_XData,LH_TD_YData] = P1_patch_indices_TD(LH_touchdown_sequence(index));
    [LF_TD_XData,LF_TD_YData] = P1_patch_indices_TD(LF_touchdown_sequence(index));
    [RF_TD_XData,RF_TD_YData] = P1_patch_indices_TD(RF_touchdown_sequence(index));
    [RH_TD_XData,RH_TD_YData] = P1_patch_indices_TD(RH_touchdown_sequence(index));

%     td_indicator = clip(LH_touchdown_sequence(index)+LF_touchdown_sequence(index)+RF_touchdown_sequence(index)+RH_touchdown_sequence(index), 0, 1);
    td_indicator = (LH_touchdown_sequence(index)+LF_touchdown_sequence(index)+RF_touchdown_sequence(index)+RH_touchdown_sequence(index))~=0;
    
 

    figure
    subplot(3,1,1)
    hold on
    pbaspect([10 1 1])

    xlim([0 length(index)])

    ylim([0 4.5])
    yticks([1 2 3 4])
    yticklabels({'RH', 'RF', 'LF', 'LH'})
    patch(LH_TD_XData,LH_TD_YData+3,'black')
    patch(LF_TD_XData,LF_TD_YData+2,'black')
    patch(RF_TD_XData,RF_TD_YData+1,'black')
    patch(RH_TD_XData,RH_TD_YData,'black')
    title('Footfall Pattern')

    subplot(3,1,2)
    plot(index,loading_force(index))
    xlim([intercal*(i-1)+1 intercal*i])
    title('loading Force')

    subplot(3,1,3)
    plot(index,td_indicator)
    xlim([intercal*(i-1)+1 intercal*i])
    ylim([0 2])



end
%% Gait Frequency Visualization
clc

% Define all gait labels and colors
% allGaitLabels = ["RGL", "RGR", "TGL", "TGR"];
% allGaitColors = [
%     0.20, 0.60, 0.80;   % RGL
%     0.90, 0.40, 0.40;   % RGR
%     0.40, 0.80, 0.40;   % TGL
%     0.90, 0.75, 0.20;   % TGR
% ];

allGaitLabels = ["TGL", "TGR", "RGL", "RGR", ];
allGaitColors = [
    0.40, 0.80, 0.40;   % TGL
    0.90, 0.75, 0.20;   % TGR
    0.20, 0.60, 0.80;   % RGL
    0.90, 0.40, 0.40;   % RGR
];

% Count gaits
gaitCounts = zeros(1, numel(allGaitLabels));
for i = 1:numel(allGaitLabels)
    gaitCounts(i) = sum(gaitTypes_midflight == allGaitLabels(i));
end

% Filter valid gaits (nonzero only)
validIdx = gaitCounts > 0;
gaitLabels = allGaitLabels(validIdx);
gaitColors = allGaitColors(validIdx, :);
gaitPercent = 100 * gaitCounts(validIdx) / sum(gaitCounts);

% Plot one bar at a time to control handles
figure; hold on;
barHandles = gobjects(1, numel(gaitLabels));
for i = 1:numel(gaitLabels)
    barHandles(i) = bar(i, gaitPercent(i), 'FaceColor', gaitColors(i, :), 'BarWidth', 0.6);
end

set(gca, 'XTick', 1:numel(gaitLabels), 'XTickLabel', gaitLabels, 'FontSize', 12);
ylabel('Percentage of Strides (%)');
xlabel('Gait Type');
title('Gait Type Distribution');
ylim([0 100]);
grid on;

% Only pass plotted bar handles to legend
legend(barHandles, gaitLabels, 'Location', 'northeastoutside');


%% Loading Force Statistics
clc

[avgForces,peakForces] = P1_LoadingForceStatistics(loading_force, stridesIds_midflight);
plot_average = 0;
plot_peak = 0;

if plot_average == 1
    % Setup right Y-axis for loading force
    yyaxis right
    ylabel('Average loading Force [N]')
    force_plot = avgForces;
else
    % Setup right Y-axis for loading force
    yyaxis right
    ylabel('Peak loading Force [N]')
    force_plot = peakForces;
end


% Compute min/max loading force for axis scaling
validForces = force_plot(~isnan(force_plot));
ymin = min(validForces);
ymax = max(validForces);
ylim([ymin, ymax]);

% Plot unscaled average loading force per gait type
numStrides = numel(force_plot);
scatterHandles = gobjects(1, numel(gaitLabels));

for i = 1:numel(gaitLabels)
    gait = gaitLabels(i);
    color = gaitColors(i, :);

    idx = find(gaitTypes_midflight == gait);
    idx = idx(idx <= numStrides);

    x = i + 0.2 * (rand(1, numel(idx)) - 0.5);  % jittered X
    y = force_plot(idx);                        % raw loading force

    scatterHandles(i) = scatter(x, y, 40, 'filled', ...
        'MarkerFaceColor', color, ...
        'MarkerEdgeColor', 'k', ...
        'MarkerFaceAlpha', 0.7, ...
        'DisplayName', gait + " (force)");
end


% Plot average force per gait as large darker dots
for i = 1:numel(gaitLabels)
    gait = gaitLabels(i);
    color = gaitColors(i, :);

    % Find valid strides of this gait
    idx = find(gaitTypes_midflight == gait);
    idx = idx(idx <= numStrides);

    if ~isempty(idx)
        meanForce = mean(force_plot(idx), 'omitnan');

        % Slightly darker color by reducing brightness
        darkColor = color * 0.6;

        % Plot larger scatter at center of gait bar (x = i)
        scatter(i, meanForce, 100, 'filled', ...
            'MarkerFaceColor', darkColor, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceAlpha', 1.0, ...
            'LineWidth', 1.2);
    end
end


% Update legend (bars + scatter)
yyaxis left  % make sure legend binds to left axis
legend([barHandles, scatterHandles], ...
       [gaitLabels + " (count)", gaitLabels + " (force)"], ...
       'Location', 'northeastoutside');
%% Gait Transition Statistics
clc
% Define gait labels and their colors
% allGaits = ["RGL", "RGR", "TGL", "TGR"];
% allGaitColors = [
%     0.20, 0.60, 0.80;   % RGL
%     0.90, 0.40, 0.40;   % RGR
%     0.40, 0.80, 0.40;   % TGL
%     0.90, 0.75, 0.20    % TGR
% ];

allGaits = ["TGL", "TGR", "RGL", "RGR", ];
allGaitColors = [
    0.40, 0.80, 0.40;   % TGL
    0.90, 0.75, 0.20;   % TGR
    0.20, 0.60, 0.80;   % RGL
    0.90, 0.40, 0.40;   % RGR
];

numGaits = numel(allGaits);

% Build all 16 transitions and count them
allTransitions = strings(1, numGaits^2);
transitionColors = zeros(numGaits^2, 3);  % Preallocate RGB triplets
k = 1;
for i = 1:numGaits
    for j = 1:numGaits
        allTransitions(k) = allGaits(i) + " → " + allGaits(j);
        transitionColors(k, :) = allGaitColors(i, :);  % Use FROM gait color
        k = k + 1;
    end
end

% Build actual transition list from gaitTypes
transitions = strings(1, numel(gaitTypes_midflight)-1);
for i = 1:numel(transitions)
    transitions(i) = gaitTypes_midflight(i) + " → " + gaitTypes_midflight(i + 1);
end

% Count occurrences
transitionCounts = zeros(1, numGaits^2);
for i = 1:numel(allTransitions)
    transitionCounts(i) = sum(transitions == allTransitions(i));
end

% Plot with bar color per source gait
figure;
b = bar(transitionCounts, 'FaceColor', 'flat');
for i = 1:numel(allTransitions)
    b.CData(i, :) = transitionColors(i, :);
end

% Axis labels
set(gca, 'XTick', 1:numel(allTransitions), ...
         'XTickLabel', allTransitions, ...
         'XTickLabelRotation', 45, ...
         'FontSize', 10);
ylabel('Number of Transitions');
xlabel('Gait Transition');
title('Gait Transition Statistics');
grid on;

% Optional: Add legend by FROM gait (one of each color only)
hold on;
legendHandles = gobjects(1, numGaits);
for i = 1:numGaits
    legendHandles(i) = plot(NaN, NaN, 's', 'MarkerFaceColor', allGaitColors(i, :), ...
        'MarkerEdgeColor', 'k', 'DisplayName', allGaits(i) + " as FROM");
end
legend(legendHandles, 'Location', 'northeastoutside');

%% Gait Transition Statistics 2
close all
clc

% % Define gait types and colors
% allGaits = ["RGL", "RGR", "TGL", "TGR"];
% allGaitColors = [
%     0.20, 0.60, 0.80;   % RGL
%     0.90, 0.40, 0.40;   % RGR
%     0.40, 0.80, 0.40;   % TGL
%     0.90, 0.75, 0.20    % TGR
% ];

allGaits = ["TGL", "TGR", "RGL", "RGR", ];
allGaitColors = [
    0.40, 0.80, 0.40;   % TGL
    0.90, 0.75, 0.20;   % TGR
    0.20, 0.60, 0.80;   % RGL
    0.90, 0.40, 0.40;   % RGR
];

numGaits = numel(allGaits);

% % Define square positions for each gait node
% positions = [
%     -1,  1;  % RGL (top-left)
%      1,  1;  % RGR (top-right)
%      1, -1;  % TGL (bottom-right)
%     -1, -1;  % TGR (bottom-left)
% ];

% Define square positions for each gait node
positions = [
    -1,  1;  % TGL (top-right)
     1,  1;  % TGR (top-left)
    -1, -1;  % RGR (bottom-right)
     1, -1;  % RGL (bottom-left)
];


% Build transition list
transitions = strings(1, numel(gaitTypes_midflight)-1);
for i = 1:numel(transitions)
    transitions(i) = gaitTypes_midflight(i) + " → " + gaitTypes_midflight(i + 1);
end

% Count and calculate per-origin percentages
transitionPercent = zeros(numGaits, numGaits);  % rows = FROM, cols = TO

for i = 1:numGaits
    from = allGaits(i);
    idx_from = find(gaitTypes_midflight(1:end-1) == from);
    total_from = numel(idx_from);

    if total_from == 0
        continue;
    end

    for j = 1:numGaits
        to = allGaits(j);
        count = sum(gaitTypes_midflight(idx_from + 1) == to);
        transitionPercent(i, j) = 100 * count / total_from;
    end
end

% --- Scaling factors ---

% Dot size scaling: total outgoing transitions relative to global total
% Count gaits
gaitCounts = zeros(1, numel(allGaits));
for i = 1:numel(allGaits)
    gaitCounts(i) = sum(gaitTypes_midflight == allGaits(i));
end

% Filter valid gaits (nonzero only)
validIdx = gaitCounts > 0;
gaitLabels = allGaits(validIdx);
gaitColors = allGaitColors(validIdx, :);
gaitPercent =  gaitCounts(validIdx) / sum(gaitCounts);
dotSizes = 50 * gaitPercent;  % size = baseline + scaled

% Fill gaitPercent for all gaits, 0 for invalid ones
gaitPercentFull = zeros(1, numGaits);
gaitPercentFull(validIdx) = gaitPercent;
dotSizesFull = zeros(1, numGaits);
dotSizesFull(validIdx) = dotSizes;



% Arrow width scaling per gait (row-wise normalization)
arrowLineWidths = zeros(numGaits, numGaits);
for i = 1:numGaits
    fromTotal = sum(transitionPercent(i, :));
    if fromTotal > 0
        % arrowLineWidths(i, :) = 0.05 + 12*gaitPercent(i) * (transitionPercent(i, :) / fromTotal);  % line width scaling

        % Nonlinear scaling: sqrt of transition %, weighted by gait %
        base = 1.5;    % minimum linewidth
        scale = 10;     % scaling factor
        
    arrowLineWidths(i, :) = base + scale * gaitPercentFull(i) .* sqrt(transitionPercent(i, :) / 100);


    end
end

% --- Plotting ---
figure; hold on; axis equal; axis off;

for i = 1:numGaits
    if validIdx(i)
        % Plot valid gait with color and size
        scatter(positions(i,1), positions(i,2), dotSizesFull(i)^2, ...
            'MarkerFaceColor', allGaitColors(i,:), ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.2);
    else
        % Plot invalid gait as gray hollow placeholder with X
        scatter(positions(i,1), positions(i,2), 80, ...
            'MarkerEdgeColor', [0.5 0.5 0.5], ...
            'MarkerFaceColor', 'none', ...
            'LineWidth', 1.2, 'Marker', 'x');
    end

    % Gait label
    arrowTheta = atan2(positions(i,2), positions(i,1));
    textOffset = 0.5 * [cos(arrowTheta), sin(arrowTheta)];

    text(positions(i,1) + textOffset(1), positions(i,2) + textOffset(2), ...
        allGaits(i), 'HorizontalAlignment', 'center', ...
        'FontSize', 12, 'FontWeight', 'bold');
end



% --- Arrows & Labels ---
for i = 1:numGaits
    for j = 1:numGaits
        fromPos = positions(i, :);
        toPos = positions(j, :);
        perc = transitionPercent(i, j);

        if perc == 0
            continue;  % Skip zero-percentage transitions
        end

        lw = arrowLineWidths(i, j);  % line width based on source gait

        if i == j
            % --- Self-transition: arc above or below dot ---
            % r = 0.3;
            r = dotSizesFull(i)/200 + 0.05;

            arrowTheta = atan2(positions(i,2), positions(i,1));
            thetaArc = linspace(1.25*pi + arrowTheta , 2.75*pi + arrowTheta ,100);
            labelOffsetDir = [positions(i,1) 0];
            % if fromPos(2) > 0  % top nodes: arc above
            %     thetaArc = linspace(0, pi, 100);  % CCW arc above
            %     labelOffsetDir = [0, 1];  % upward
            % else  % bottom nodes: arc below
            %     thetaArc = linspace(pi, 2*pi, 100);  % CCW arc below
            %     labelOffsetDir = [0, -1];  % downward
            % end
            
            xc = fromPos(1) + r * cos(thetaArc);
            yc = fromPos(2) + r * sin(thetaArc);
            plot(xc, yc, 'Color', allGaitColors(i,:), 'LineWidth', lw);
            
            % Arrowhead at arc end using filled triangle (scaled with lw)
            arrowEnd = [xc(end), yc(end)];
            dirVec = [xc(end) - xc(end-1), yc(end) - yc(end-1)];
            dirVec = dirVec / norm(dirVec);
            
            headLength = 0.02 * lw;  % length of arrowhead
            headWidth  = 0.01 * lw;  % half-width of base
            
            sideVec = headWidth * [-dirVec(2), dirVec(1)];  % perpendicular for width
            tip   = arrowEnd + headLength * dirVec;
            left  = arrowEnd + sideVec;
            right = arrowEnd - sideVec;
            
            fill([tip(1), left(1), right(1)], [tip(2), left(2), right(2)], ...
                 allGaitColors(i,:), 'EdgeColor', 'none');


            
            % Label near center of arc
            midIdx = round(length(xc)/2);
            labelPos = [xc(midIdx), yc(midIdx)] + 0.2 * labelOffsetDir;
            text(labelPos(1), labelPos(2), sprintf('%.2f%%', perc), ...
                'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', allGaitColors(i,:));

        else
            % --- Shifted straight arrows for i ≠ j ---
            vec = toPos - fromPos;
            unitVec = vec / norm(vec);
            perpVec = [-unitVec(2), unitVec(1)];  % 90° rotation

            shiftAmount = 0.05;
            offset = shiftAmount * perpVec;

            % Start and end points
            startPos = fromPos + 0.15 * vec + offset;
            endPos = toPos - 0.15 * vec + offset;
            
            % Draw triangle arrowhead at end
            dirVec = endPos - startPos;
            dirVec = dirVec / norm(dirVec);
            
            headLength = 0.03 * lw;
            headWidth  = 0.015 * lw;
            
            sideVec = headWidth * [-dirVec(2), dirVec(1)];
            tip   = endPos;
            left  = endPos - headLength * dirVec + sideVec;
            right = endPos - headLength * dirVec - sideVec;
            
            fill([tip(1), left(1), right(1)], [tip(2), left(2), right(2)], ...
                 allGaitColors(i,:), 'EdgeColor', 'none');
            
            trimmedEnd = endPos - headLength * dirVec;
            
            plot([startPos(1), trimmedEnd(1)], [startPos(2), trimmedEnd(2)], ...
                 'Color', allGaitColors(i,:), 'LineWidth', lw);

            % Label placement (avoid overlap)
            midPos = (startPos + endPos) / 2;
            arrowDir = (toPos - fromPos);
            arrowDir = arrowDir / norm(arrowDir);
            perpVec = [-arrowDir(2), arrowDir(1)];
            offsetPerp = 0.12 * perpVec;
            offsetSource = 0.40 * (fromPos - midPos) / norm(fromPos - midPos);
            labelPos = midPos + offsetPerp + offsetSource;

            angleDeg = atan2d(arrowDir(2), arrowDir(1));
            text(labelPos(1), labelPos(2), ...
                sprintf('%.2f%%', perc), ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 10, ...
                'Color', allGaitColors(i,:), ...
                'Rotation', angleDeg);

        end
    end
end

% Move title above the topmost point
title({'Gait Transition Percentages', 'For Alice'}, 'FontSize', 14);
% titleHandle.Units = 'data';
% titleHandle.Position(2) = max(positions(:,2)) + 0.3;  % adjust vertical position