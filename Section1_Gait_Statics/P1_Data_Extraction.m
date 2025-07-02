%% Data Extraction

%% Extract the data from Excel
clear
clc
% select horizontal position first
[filename_data, pathname_data] = uigetfile('*.csv','Select a File');
raw_data = P1_ImportFile(filename_data);

% extract dog name
dog_name_index = strfind(filename_data,'_');
dog_name       = filename_data(dog_name_index(1)+1:dog_name_index(2)-1);



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

if strcmp(dog_name,"Alice")

    loading_force = 0.00259*(load + 6027.006); % for Alice

elseif strcmp(dog_name,"Chips")

    loading_force = 0.00262*(load-86579.096); % for Chips
end

% extract horizontal velocity and pitching velocity

a_horizontal          = raw_data(:,find(strcmp(raw_data.Properties.VariableNames,'ay'))).ay;
a_vertical            = raw_data(:,find(strcmp(raw_data.Properties.VariableNames,'az'))).az;
o_pitching            = raw_data(:,find(strcmp(raw_data.Properties.VariableNames,'gx'))).gx;

% normalize data:       ay,az units: g          gx unit: degree/s

% normalize angular velocity data
o_pitching = deg2rad(o_pitching);


%% Gait Index and Types Buffer
stridesequences = [LH_touchdown_sequence LF_touchdown_sequence RF_touchdown_sequence RH_touchdown_sequence];
[stridesIds, stridesIds_midflight, gaitTypes, gaitTypes_midflight] = P1_StrideIndexExtractor(stridesequences);

%% Creating Template Gait Data (Averaged, Normalized & Envelopes including ±1σ)
close all; clc;

% Data Normalization parameters
if strcmp(dog_name, "Alice")
    m = 30;       % Mass of Alice: 30kg
    l_leg = 0.5;  % Leg length of Alice: 50cm
elseif strcmp(dog_name, "Chips")
    m = 25;       % Mass of Chips: 25kg
    l_leg = 0.48; % Leg length of Chips: 48cm
end

g = 9.81; % gravity constant

% Normalization factors
T_scale = sqrt(l_leg / g);
a_scale = m;
o_scale = sqrt(g / l_leg);
F_scale = m * g;

gaitLabels = ["RGL", "RGR", "TGL", "TGR"];
numGaits = numel(gaitLabels);
GaitDataTemplates_Averaged_Normalized = struct();

% Extract stride indices & gait types
stridesequences = [LH_touchdown_sequence LF_touchdown_sequence RF_touchdown_sequence RH_touchdown_sequence];
[stridesIds, stridesIds_midflight, gaitTypes, gaitTypes_midflight] = P1_StrideIndexExtractor(stridesequences);

for g = 1:numGaits
    gait    = gaitLabels(g);
    gaitIdx = find(gaitTypes_midflight == gait);
    if isempty(gaitIdx), continue; end

    % Collect per-stride signals & footfall timing
    L  = {}; AH = {}; AV = {}; OP = {}; FT = {};
    strideLengths = zeros(1, numel(gaitIdx));
    for k = 1:numel(gaitIdx)
        range    = stridesIds_midflight(:, gaitIdx(k));
        idxStart = range(1);
        idxEnd   = range(2);
        len      = idxEnd - idxStart + 1;
        strideLengths(k) = len;

        L{end+1}  = loading_force(idxStart:idxEnd);
        AH{end+1} = a_horizontal(idxStart:idxEnd);
        AV{end+1} = a_vertical(idxStart:idxEnd);
        OP{end+1} = o_pitching(idxStart:idxEnd);

        % Footfall timing (normalized indices)
        footfall1 = stridesequences(idxStart:idxEnd, :);
        ft1 = nan(1,8);
        for limb = 1:4
            stance = footfall1(:, limb);
            td_idx = find(diff([0; stance]) == 1, 1);
            lo_idx = find(diff([stance; 0]) == -1, 1);
            ft1(2*limb-1) = td_idx / len;
            ft1(2*limb)   = lo_idx / len;
        end
        FT{end+1} = ft1;
    end

    avgStrideLen = round(mean(strideLengths));

    % Resample each stride to common length
    Lr  = cellfun(@(x) DataLengthResampler(x, 1:avgStrideLen), L,  'UniformOutput', false);
    AHr = cellfun(@(x) DataLengthResampler(x, 1:avgStrideLen), AH, 'UniformOutput', false);
    AVr = cellfun(@(x) DataLengthResampler(x, 1:avgStrideLen), AV, 'UniformOutput', false);
    OPr = cellfun(@(x) DataLengthResampler(x, 1:avgStrideLen), OP, 'UniformOutput', false);

    Lmat  = cell2mat(cellfun(@(x) x(:), Lr,  'UniformOutput', false));
    AHmat = cell2mat(cellfun(@(x) x(:), AHr, 'UniformOutput', false));
    AVmat = cell2mat(cellfun(@(x) x(:), AVr, 'UniformOutput', false));
    OPmat = cell2mat(cellfun(@(x) x(:), OPr, 'UniformOutput', false));

    % Compute statistics and normalize
    L_mean  = mean(Lmat,  2, 'omitnan') / F_scale;
    AH_mean = mean(AHmat, 2, 'omitnan') / a_scale;
    AV_mean = mean(AVmat, 2, 'omitnan') / a_scale;
    OP_mean = mean(OPmat, 2, 'omitnan') / o_scale;

    L_std  = std(Lmat,  0, 2, 'omitnan') / F_scale;
    AH_std = std(AHmat, 0, 2, 'omitnan') / a_scale;
    AV_std = std(AVmat, 0, 2, 'omitnan') / a_scale;
    OP_std = std(OPmat, 0, 2, 'omitnan') / o_scale;

    L_min   = min(Lmat,  [], 2, 'omitnan') / F_scale;
    L_max   = max(Lmat,  [], 2, 'omitnan') / F_scale;
    AH_min  = min(AHmat, [], 2, 'omitnan') / a_scale;
    AH_max  = max(AHmat, [], 2, 'omitnan') / a_scale;
    AV_min  = min(AVmat, [], 2, 'omitnan') / a_scale;
    AV_max  = max(AVmat, [], 2, 'omitnan') / a_scale;
    OP_min  = min(OPmat, [], 2, 'omitnan') / o_scale;
    OP_max  = max(OPmat, [], 2, 'omitnan') / o_scale;

    % Footfall timing statistics
    ft_all     = cell2mat(FT');                % [numStrides × 8]
    ft_average = mean(ft_all, 1, 'omitnan');   % 1×8
    ft_std     = std(ft_all, 0, 1, 'omitnan'); % 1×8 ±1σ
    ft_min     = min(ft_all,  [], 1, 'omitnan');
    ft_max     = max(ft_all,  [], 1, 'omitnan');

    % Normalized time vector
    avg_time = ((0:avgStrideLen-1) * dt) / T_scale;

    % Store into struct
    S = struct();
    S.t_exp                   = avg_time;
    S.loading_force_mean      = L_mean;
    S.loading_force_std       = L_std;
    S.loading_force_min       = L_min;
    S.loading_force_max       = L_max;
    S.a_horizontal_mean       = AH_mean;
    S.a_horizontal_std        = AH_std;
    S.a_horizontal_min        = AH_min;
    S.a_horizontal_max        = AH_max;
    S.a_vertical_mean         = AV_mean;
    S.a_vertical_std          = AV_std;
    S.a_vertical_min          = AV_min;
    S.a_vertical_max          = AV_max;
    S.o_pitching_mean         = OP_mean;
    S.o_pitching_std          = OP_std;
    S.o_pitching_min          = OP_min;
    S.o_pitching_max          = OP_max;
    S.ft_exp                  = ft_average;
    S.ft_std                  = ft_std;
    S.ft_min                  = ft_min;
    S.ft_max                  = ft_max;

    GaitDataTemplates_Averaged_Normalized.(gait) = S;
end

% Visualization with Min/Max & ±1σ Envelopes
figure;
for g = 1:numGaits
    gait = gaitLabels(g);
    if ~isfield(GaitDataTemplates_Averaged_Normalized, gait)
        continue;
    end

    S = GaitDataTemplates_Averaged_Normalized.(gait);
    t = S.t_exp;
    cols = lines(4);

    subplot(2,2,g);
    hold on;

    % Loading force envelopes
    fill([t, fliplr(t)], [S.loading_force_min', fliplr(S.loading_force_max')], cols(1,:), 'FaceAlpha',0.2, 'EdgeColor','none');
    fill([t, fliplr(t)], [S.loading_force_mean'+S.loading_force_std', fliplr(S.loading_force_mean'-S.loading_force_std')], cols(1,:), 'FaceAlpha',0.3, 'EdgeColor','none');
    plot(t, S.loading_force_mean, 'LineWidth',1.5,'Color',cols(1,:));

    % Horizontal accel envelopes
    fill([t, fliplr(t)], [S.a_horizontal_min', fliplr(S.a_horizontal_max')], cols(2,:), 'FaceAlpha',0.2,'EdgeColor','none');
    fill([t, fliplr(t)], [S.a_horizontal_mean'+S.a_horizontal_std', fliplr(S.a_horizontal_mean'-S.a_horizontal_std')], cols(2,:), 'FaceAlpha',0.3,'EdgeColor','none');
    plot(t, S.a_horizontal_mean, 'LineWidth',1.5,'Color',cols(2,:));

    % Vertical accel envelopes
    fill([t, fliplr(t)], [S.a_vertical_min', fliplr(S.a_vertical_max')], cols(3,:), 'FaceAlpha',0.2,'EdgeColor','none');
    fill([t, fliplr(t)], [S.a_vertical_mean'+S.a_vertical_std', fliplr(S.a_vertical_mean'-S.a_vertical_std')], cols(3,:), 'FaceAlpha',0.3,'EdgeColor','none');
    plot(t, S.a_vertical_mean, 'LineWidth',1.5,'Color',cols(3,:));

    % Pitching rate envelopes
    fill([t, fliplr(t)], [S.o_pitching_min', fliplr(S.o_pitching_max')], cols(4,:), 'FaceAlpha',0.2,'EdgeColor','none');
    fill([t, fliplr(t)], [S.o_pitching_mean'+S.o_pitching_std', fliplr(S.o_pitching_mean'-S.o_pitching_std')], cols(4,:), 'FaceAlpha',0.3,'EdgeColor','none');
    plot(t, S.o_pitching_mean, 'LineWidth',1.5,'Color',cols(4,:));

    title(gait);
    xlabel('Normalized Time');
    legend({'min/max env','±1σ env','mean'}, 'Location','best');
    grid on;
end
sgtitle('Normalized Averaged Gait Signals with Min/Max & ±1σ Envelopes');





%% Save template gaits

save( sprintf('P1_GaitDataTemplates_%s.mat', dog_name), 'GaitDataTemplates_Averaged_Normalized' );

%% Simplified Normalization & Visualization with Struct Packing for Two Strides

stride_idx = 22;
% transition_type = "TGR-RGL";
% transition_type = "RGL-TGR";

transition_type = "RGL-TGL";

% Intervals for two strides
int1 = stridesIds_midflight(1, stride_idx):stridesIds_midflight(2, stride_idx);
int2 = stridesIds_midflight(1, stride_idx+1):stridesIds_midflight(2, stride_idx+1);

% Constants
% For Alice
m = 30;     % kg
l_leg = 0.5; % m

g = 9.81;   % m/s^2

% Normalization factors
T_scale = sqrt(l_leg/g);
F_scale = m*g;
a_scale = m;
o_scale = sqrt(g/l_leg);

% Helper function for normalization
norm_data = @(data, interval, scale) data(interval) / scale;

% Normalized stride data (2x1 cells)
TransitionTemplate_Normalized.loading_force_exp = {norm_data(loading_force, int1, F_scale); norm_data(loading_force, int2, F_scale)};
TransitionTemplate_Normalized.a_horizontal_exp  = {norm_data(a_horizontal, int1, a_scale); norm_data(a_horizontal, int2, a_scale)};
TransitionTemplate_Normalized.a_vertical_exp    = {norm_data(a_vertical, int1, a_scale); norm_data(a_vertical, int2, a_scale)};
TransitionTemplate_Normalized.o_pitching_exp    = {norm_data(o_pitching, int1, o_scale); norm_data(o_pitching, int2, o_scale)};

% Footfall timing extraction to match original format
stride_len = [length(int1), length(int2)];
TransitionTemplate_Normalized.ft_exp = zeros(2,8);
for k = 1:2
    interval = eval(sprintf('int%d', k));
    footfall = stridesequences(interval, :);
    for limb = 1:4
        stance = footfall(:, limb);
        td_idx = find(diff([0; stance]) == 1, 1);
        lo_idx = find(diff([stance; 0]) == -1, 1);
        if isempty(td_idx), td_idx = NaN; end
        if isempty(lo_idx), lo_idx = NaN; end
        TransitionTemplate_Normalized.ft_exp(k, 2*limb-1:2*limb) = [td_idx, lo_idx] / stride_len(k);
    end
end
TransitionTemplate_Normalized.ft_exp(2,:) = TransitionTemplate_Normalized.ft_exp(2,:) + 1;
PlotFootfallSequence_Hildebrand(TransitionTemplate_Normalized.ft_exp,"Experimental Footfall Sequences")

% Normalized time vectors (2x1 cells)
DT = dt / T_scale;
TransitionTemplate_Normalized.t_exp = {(0:stride_len(1)-1)*DT; (0:stride_len(2)-1)*DT};
TransitionTemplate_Normalized.trans_type = transition_type;

% Plotting
figure;
for k = 1:2
    subplot(2,1,k);
    hold on;
    plot(TransitionTemplate_Normalized.t_exp{k}, TransitionTemplate_Normalized.loading_force_exp{k}, 'LineWidth',1.5, 'DisplayName', 'Loading Force');
    plot(TransitionTemplate_Normalized.t_exp{k}, TransitionTemplate_Normalized.a_horizontal_exp{k}, 'LineWidth',1.5, 'DisplayName', 'Accel. Horizontal');
    plot(TransitionTemplate_Normalized.t_exp{k}, TransitionTemplate_Normalized.a_vertical_exp{k}, 'LineWidth',1.5, 'DisplayName', 'Accel. Vertical');
    plot(TransitionTemplate_Normalized.t_exp{k}, TransitionTemplate_Normalized.o_pitching_exp{k}, 'LineWidth',1.5, 'DisplayName', 'Pitching Rate');
    xlabel('Normalized Time');
    ylabel('Normalized Values');
    title(sprintf('Stride %d (%s)', stride_idx + (k-1), transition_type));
    legend;
    grid on;
end



%% Save template gaits

% save('P1_TransitionTemplate_Alice_TGR_RGL.mat',"TransitionTemplate_Normalized")
% save('P1_TransitionTemplate_Alice_RGL_TGR.mat',"TransitionTemplate_Normalized")
save('P1_TransitionTemplate_Alice_RGL_TGL.mat',"TransitionTemplate_Normalized")

