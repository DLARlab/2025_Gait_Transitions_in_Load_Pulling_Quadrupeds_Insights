function Section2_Single_Stride_Replication()
    clc;
    global X_accum term_weights gait_data  T Y P TF_Sim  R2  graphicsObjs currentfolder isSimulating;
    X_accum = [];
    term_weights = struct();
    gait_data = struct();
    T = []; Y = []; P = []; TF_Sim = [];
    R2 = struct();
    graphicsObjs = struct();
    currentfolder = pwd;
    isSimulating = false;

    %% GUI Setup
    screenSize = get(groot, 'ScreenSize');
    fig = uifigure('Name','SLIP Quadruped Viewer', ...
                   'Position', [screenSize(3)*0.5-screenSize(4)*0.6,  screenSize(4)*0.5 - screenSize(4)*0.4, ...
                                screenSize(4)*1.2, screenSize(4)*0.80]);

    % Define proportions
    buttonWidth = 0.08 * screenSize(4);
    buttonHeight = 0.04 * screenSize(4);
    topY = screenSize(4)*0.80 - 1.5*buttonHeight;

    % Folder button
    folderButton = uibutton(fig, 'Text', 'Select Folder', ...
        'Position', [0.10*screenSize(4), topY, buttonWidth, buttonHeight], ...
        'ButtonPushedFcn', @(src, event) SelectFolder());

    % Dropdown
    dropdownWidth = 0.12 * screenSize(4);
    matDropdown = uidropdown(fig, ...
        'Position', [0.10*screenSize(4) + buttonWidth + 0.05*buttonWidth, topY, dropdownWidth, buttonHeight], ...
        'Items', GetMatFiles(), ...
        'ValueChangedFcn', @(src, event) LoadSelectedFile());

    % Slider (no ticks or labels)
    sliderX = 0.08*screenSize(4) + buttonWidth + dropdownWidth + 3.5*buttonWidth;
    slider = uislider(fig, ...
        'Position', [sliderX, topY + buttonHeight/2, 3 * buttonWidth, 3], ...
        'Limits', [0, 1], ...
        'ValueChangingFcn', @(s, e) HandleSliderChange(e.Value), ...
        'MajorTicks', [], 'MinorTicks', []);

    % Tick labels manually
    tickLabels = {'0', '0.5', '1'};
    tickPositions = linspace(sliderX, sliderX + 3 * buttonWidth, 3);
    for i = 1:length(tickLabels)
        uilabel(fig, 'Text', tickLabels{i}, 'Position', [tickPositions(i) - 10, topY + buttonHeight/2 - 20, 30, 20]);
    end

    % Normalized time input
    normTimeInput = uieditfield(fig, 'numeric', ...
        'Position', [slider.Position(1) + slider.Position(3) + 0.1*buttonWidth, topY, 0.5*buttonWidth, buttonHeight], ...
        'Limits', [0, 1], ...
        'ValueDisplayFormat', '%.2f', ...
        'ValueChangedFcn', @(src, event) HandleNormTimeInput());

    uilabel(fig, 'Text', '/ Stride Cycle', ...
        'Position', [normTimeInput.Position(1) + normTimeInput.Position(3) + 0.05*buttonWidth, topY, 1.25*buttonWidth, buttonHeight]);

    % Run Simulation button
    runButton = uibutton(fig, 'Text', 'Run Simulation', ...
        'Position', [normTimeInput.Position(1) + normTimeInput.Position(3) + 1.25*buttonWidth + 0.05*buttonWidth, topY, 1.1*buttonWidth, buttonHeight], ...
        'ButtonPushedFcn', @(src, event) RunSimulation());

    % Stop Simulation button
    stopButton = uibutton(fig, 'Text', 'Stop Simulation', ...
        'Position', [runButton.Position(1) + runButton.Position(3) + 0.05*buttonWidth, topY, 1.1*buttonWidth, buttonHeight], ...
        'ButtonPushedFcn', @(src, event) StopSimulation());


    %% Axes
    ax = gobjects(1,6);
    figSize = fig.Position;
    PlotSize = [(3.5/10)*figSize(3),         (12/16)*(3.5/10)*figSize(3); ...
                (3.5/2.5)*figSize(3)/5,        (12/16)*(3.5/2.5)*(figSize(3)/5)];
    PlotPositions = [
        (1/2)*figSize(3)-(2/2)*PlotSize(1,1) - 1.0*buttonWidth,           (0.5/10)*figSize(4)+PlotSize(2,2) + buttonHeight,       1.5*PlotSize(1,1),               PlotSize(1,2);
        (1/2)*figSize(3)-(3/2)*PlotSize(2,1) - 1.5*buttonWidth,           (0.5/10)*figSize(4),                                    2.1*PlotSize(2,1),                 PlotSize(2,2);
        (1/2)*figSize(3)+(1/2)*PlotSize(2,1) - 0.0*buttonWidth,           (0.5/10)*figSize(4),                                    PlotSize(2,1),                   PlotSize(2,2)];

    for axidx = 1:size(PlotPositions,1)
        ax(axidx) = uiaxes(fig);
        ax(axidx).Position = PlotPositions(axidx, :);
    end

    % Readout boxes for R2 values
    baseX = ax(1).Position(1) + ax(1).Position(3) + 2*buttonWidth;
    boxWidth = 2*buttonWidth;
    boxHeight = buttonHeight;
    spacingY = 1.2*buttonHeight;
    baseY = ax(1).Position(2) + 0.5*ax(1).Position(4) - 0.5*boxHeight ;

    uilabel(fig, 'Text', '$R^2$. Footfall Timing', 'Interpreter', 'latex', 'Position', [baseX, baseY + 3*spacingY, boxWidth, boxHeight]);
    footfallTimingField = uieditfield(fig, 'numeric', 'Editable', 'off', 'Position', [baseX, baseY + 2*spacingY, boxWidth, boxHeight]);

    uilabel(fig, 'Text', '$R^2$. Loading Force', 'Interpreter', 'latex', 'Position', [baseX, baseY + 1*spacingY, boxWidth, boxHeight]);
    loadingForceField = uieditfield(fig, 'numeric', 'Editable', 'off', 'Position', [baseX, baseY + 0*spacingY, boxWidth, boxHeight]);

    uilabel(fig, 'Text', '$R^2$. Stride Duration', 'Interpreter', 'latex', 'Position', [baseX, baseY - spacingY, boxWidth, boxHeight]);
    strideDurationField = uieditfield(fig, 'numeric', 'Editable', 'off', 'Position', [baseX, baseY - 2*spacingY, boxWidth, boxHeight]);

    %% File dropdown
    initialFiles = GetMatFiles();
    matDropdown.Items = initialFiles;
    if ~strcmp(initialFiles{1}, '<none>')
        matDropdown.Value = initialFiles{1};
        LoadSelectedFile();
    end

    function files = GetMatFiles()
        matFiles = dir(fullfile(currentfolder, '*.mat'));
        validFiles = {};
        for i = 1:length(matFiles)
            s = load(fullfile(currentfolder, matFiles(i).name));
            if isfield(s, 'X_accum')
                validFiles{end+1} = matFiles(i).name;
            end
        end
        if isempty(validFiles)
            files = {'<none>'};
        else
            files = validFiles;
        end
    end

    function SelectFolder()
        folder = uigetdir(currentfolder);  % Start from previous folder
        if folder ~= 0
            currentfolder = folder;
            matDropdown.Items = GetMatFiles();
            if ~strcmp(matDropdown.Items{1}, '<none>')
                matDropdown.Value = matDropdown.Items{1};
                LoadSelectedFile();
            else
                X_accum = [];
            end
            figure(fig);
        end
    end

    function LoadSelectedFile()
        file = matDropdown.Value;
        if strcmp(file, '<none>'), return; end
        data = load(fullfile(currentfolder, file));
        if ~isfield(data, 'X_accum')
            uialert(fig, 'Selected file does not contain "X_accum".', 'Error');
            return;
        end
        X_accum = data.X_accum;
        gait_data = data.gait_data;
        term_weights = data.term_weights;
        InitializePlots();
    end

    function InitializePlots()
        % Load .mat file containing 'X_accum'
        if ~exist('X_accum', 'var')
            [file, path] = uigetfile('*.mat', 'Select a MAT-file containing X_accum');
            if isequal(file, 0)
                uialert(fig, 'File selection cancelled.', 'Warning');
                return;
            end
            loadedData = load(fullfile(path, file));
            if ~isfield(loadedData, 'X_accum')
                uialert(fig, 'Selected file does not contain "X_accum".', 'Error');
                return;
            end
            X_accum = loadedData.X_accum;
            gait_data = loadedData.gait_data;
            term_weights = loadedData.term_weights;
        end
    
        % Ensure X_accum is valid
        if isempty(X_accum) || ~isvector(X_accum)
            uialert(fig, 'Invalid solution format in X_accum.', 'Error');
            return;
        end
    
        clc;
    
        % Run simulation
        [~, T, Y, P, ~, TF_Sim, ~] = Quad_Load_ZeroFun_Transition_v2(X_accum);        
        ft_sim = P([1 2 3 4 7 8 5 6])./P(9);
        
        Number_of_Strides = size(TF_Sim,2);

        t_exp             = gait_data.t_exp;

        % Extract all fields from gait_data that contain 'loading_force'
        lf_fields = fieldnames(gait_data);
        lf_fields = lf_fields(contains(lf_fields, 'loading_force'));
        if numel(lf_fields) >1
            tf_exp = struct();
            for j = 1:numel(lf_fields)
                tf_exp.(lf_fields{j}) = gait_data.(lf_fields{j});
            end
        else
            tf_exp = gait_data.loading_force_exp;
        end
        
        % Extract all fields from gait_data that contain 'ft'
        ft_fields = fieldnames(gait_data);
        ft_fields = ft_fields(contains(ft_fields, 'ft'));
        if numel(ft_fields)>1
            ft_exp = struct();
            for k = 1:numel(ft_fields)
                ft_exp.(ft_fields{k}) = gait_data.(ft_fields{k});
            end
        else
            ft_exp = gait_data.ft_exp;
        end
        
        if isstruct(ft_exp)
            [~, ~, R2] = fms_NStridesObjectiveFcn_Quad_Load_v2(X_accum, t_exp, ft_exp.ft_exp, tf_exp.loading_force_mean, term_weights);
        else
            [~, ~, R2] = fms_NStridesObjectiveFcn_Quad_Load_v2(X_accum, t_exp, ft_exp, tf_exp, term_weights);
        end

        % Update readout boxes
        footfallTimingField.Value = R2.footfalltiming;
        loadingForceField.Value = R2.loadingforce;
        strideDurationField.Value = R2.strideduration;

        % Reset GUI inputs
        slider.Value = 0;
        normTimeInput.Value = 0;
    
    
        % Animation options
        options = struct('AnimationMode', 'Detailed');
        options = DefaultOptions(options);
    
        % Create graphics objects
        graphicsObjs.Animation = SLIP_Animation_Quad_Load(Y, P, ax(1).Position, ax(1), options);
        graphicsObjs.Footfall = SLIP_FootfallSequence_Quad( ft_sim, ft_exp, ax(2).Position, ax(2));
        graphicsObjs.TuglineForce = SLIP_TuglineForce( TF_Sim, tf_exp, Number_of_Strides, ax(3).Position, ax(3));

        % graphicsObjs.Trajectories = SLIP_Trajectories_Quad(T, Y, ax(3));
    end
    
    function RunSimulation()
        if isempty(T) || isempty(Y), return; end
        isSimulating = true;
        startNorm = slider.Value;
        totalFrames = round((1 - startNorm) * 120);
        frameNorms = linspace(startNorm, 1, totalFrames);
        for k = 1:totalFrames
            if ~isSimulating
                break;
            end
            t = frameNorms(k) * T(end);

            slider.Value = frameNorms(k);
            normTimeInput.Value = round(frameNorms(k), 2);

            UpdateFrame(t);
%             pause(1/30);  % ~30 fps
        end
    end
    
    function StopSimulation()
        isSimulating = false;
    end

    function UpdateFrame(t)
        if isempty(T)
            return; 
        end

        y = interp1(T' + linspace(0,1e-5,length(T)), Y, t, 'linear', 'extrap');

        T_ = T(T <= t);
        Y_ = Y(T <= t,:);

        if ~isempty(graphicsObjs.Animation)
            graphicsObjs.Animation.update(t, y, P); 
        end

    end

    function HandleSliderChange(val)
        val = min(max(val, 0), 1);
        val = round(val, 2);
        normTimeInput.Value = val;
        UpdateFrame(val * T(end));
    end

    function HandleNormTimeInput()
        val = normTimeInput.Value;
        val = min(max(val, 0), 1);
        val = round(val, 2);
        normTimeInput.Value = val;
        slider.Value = val;
        UpdateFrame(val * T(end));
    end

    function options = DefaultOptions(options)
        defaults = struct(...
            'ShowAnimation', 'On', 'ShowTrajectories', 'On', 'ShowPeriodicOrbit', 'On', ...
            'ShowGRF', 'On', 'ShowPhaseDiagram', 'Off', 'SaveAnimation', 'Off', ...
            'SaveTrajectories', 'Off', 'SavePeriodicOrbit', 'Off', 'SaveGRF', 'Off', ...
            'AnimationRepeatTime', '1');
        fields = fieldnames(defaults);
        for i = 1:length(fields)
            if ~isfield(options, fields{i}) || isempty(options.(fields{i}))
                options.(fields{i}) = defaults.(fields{i});
            end
        end
    end
end
