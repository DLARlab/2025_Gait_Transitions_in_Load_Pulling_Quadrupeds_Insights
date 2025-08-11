function Section3_Gait_Transition_Replication()
    clc;
    global X_accum term_weights gait_data sensitivity_data  T Y P TF_Sim  R2  graphicsObjs currentfolder isSimulating;
    X_accum = [];
    term_weights = struct();
    gait_data = struct();
    sensitivity_data = struct();
    T = []; Y = []; P = []; TF_Sim = [];
    R2 = struct();
    graphicsObjs = struct();
    currentfolder = pwd;
    isSimulating = false;

    %% GUI Setup
    screenSize = get(groot, 'ScreenSize');
    fig = uifigure('Name','SLIP Quadruped Viewer', ...
                   'Position', [screenSize(3)*0.5-screenSize(4)*0.6,  screenSize(4)*0.5 - screenSize(4)*0.45, ...
                                screenSize(4)*1.0, screenSize(4)*0.60]);

    % Define proportions
    figPos = fig.Position;
    buttonWidth = 0.06 * screenSize(4);
    buttonHeight = 0.03 * screenSize(4);
    topY = figPos(4) - 1.5*buttonHeight;

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
    for t = 1:length(tickLabels)
        uilabel(fig, 'Text', tickLabels{t}, 'Position', [tickPositions(t) - 10, topY + buttonHeight/2 - 20, 30, 20]);
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
    figW = figPos(3); figH = figPos(4);
    mX = 0.06*figW; mY = 0.06*figH; gX = 0.02*figW; gY = 0.02*figH; gX12 = gX;
    usableH = max(0, (topY - gY) - mY);
    leftW  = (figW - 2*mX - gX)/2; rightW = leftW;
    leftX  = mX; rightX = leftX + leftW + gX;

    % Keep left column width, update height for 3:1 ratio
    wL = leftW;
    hL = wL / 3;
    totalLeftHeight = 3*hL + 2*gY;
    yTop = mY + totalLeftHeight;

    y1 = mY; y2 = y1 + hL + gY; y3 = y2 + hL + gY;
    posL = [leftX, y3, 1.1*wL, hL; leftX, y2, wL, hL; leftX, y1, wL, hL];

    % Right column as before
    coef = (1 + 1/1.2);
    rhs_v = usableH - gY + (gX12/2)/1.2; bound_v = rhs_v / coef;
    rhs_l = (3/4)*leftW + gY + (gX12/2)/1.2; bound_l = rhs_l / coef;
    rowH4 = min([rightW/2, bound_v, bound_l]); rowH4 = max(rowH4, gX12/2 + eps);
    w4  = 2*rowH4; w56 = (w4 - gX12)/2; h56 = w56 / 1.2;
    x4 = rightX + (rightW - w4)/2; pos4 = [x4, yTop - rowH4, w4, rowH4];
    x5 = x4; x6 = x5 + w56 + gX12; pos5 = [x5, mY, w56, h56]; pos6 = [x6, mY, w56, h56];

    PlotPositions = [posL; pos4; pos5; pos6];
    for ax_id = 1:6
        ax(ax_id) = uiaxes(fig);
        ax(ax_id).Position = PlotPositions(ax_id,:);
    end



    % Readout boxes for R2 values
    baseX = x5;
    baseY = mY + h56;
    boxWidth = 2*buttonWidth;
    boxHeight = buttonHeight;



    uilabel(fig, 'Text', '$R^2$. Weighted', 'Interpreter', 'latex', 'Position', [baseX + 0.5*boxWidth, baseY + boxHeight, boxWidth, boxHeight]);
    R2WeightedField = uieditfield(fig, 'numeric', 'Editable', 'off', 'Position', [baseX + 1.5*boxWidth, baseY +  boxHeight, boxWidth, boxHeight]);


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
        for m = 1:length(matFiles)
            s = load(fullfile(currentfolder, matFiles(m).name));
            if isfield(s, 'X_accum')
                validFiles{end+1} = matFiles(m).name;
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
        gait_data = data.TransitionTemplate_Normalized;
        term_weights = data.term_weights;
        sensitivity_data = data.SensitivityStudyData;
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
            gait_data = loadedData.TransitionTemplate_Normalized;
            term_weights = loadedData.term_weights;
            sensitivity_data = data.SensitivityStudyData;
        end
    
        % Ensure X_accum is valid
        if isempty(X_accum) || ~isvector(X_accum)
            uialert(fig, 'Invalid solution format in X_accum.', 'Error');
            return;
        end
    
        clc;
    
        % Run simulation
%         [~, T, Y, P, ~, TF_Sim, ~] = Quad_Load_ZeroFun_Transition_v2(X_accum); 
        [~, ~, T, Y, ~, TF_Sim, P, ~, Number_of_Strides] = SimulateQuadLoadStrides(X_accum);
       
        ft_sim = P(:,[1 2 3 4 7 8 5 6])./P(:,9) + ([1:Number_of_Strides]'-1);


        
        t_exp             = gait_data.t_exp;

        % Extract all fields from gait_data that contain 'loading_force'
        lf_fields = fieldnames(gait_data);
        lf_fields = lf_fields(contains(lf_fields, 'loading_force'));
        if numel(lf_fields) >1
            lf_exp = struct();
            for j = 1:numel(lf_fields)
                lf_exp.(lf_fields{j}) = gait_data.(lf_fields{j});
            end
        else
            lf_exp = gait_data.loading_force_exp;
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
            [~, ~, R2] = fms_NStridesObjectiveFcn_Quad_Load_v2(X_accum, t_exp, ft_exp.ft_exp, lf_exp.loading_force_mean, term_weights);
        else
            [~, ~, R2] = fms_NStridesObjectiveFcn_Quad_Load_v2(X_accum, t_exp, ft_exp, lf_exp, term_weights);
        end

        % Update readout boxes
        R2WeightedField.Value = R2.weighted;


        % Reset GUI inputs
        slider.Limits = [0 Number_of_Strides];
        slider.Value = 0;
        normTimeInput.Value = 0;
    
    
        % Animation options
        options = struct('AnimationMode', 'Detailed');
        options = DefaultOptions(options);
    
        % Create graphics objects
        graphicsObjs.Footfall = SLIP_FootfallSequence_Quad( ft_sim, ft_exp, ax(1).Position, ax(1));
        graphicsObjs.LegTrajectories = SLIP_LegTrajectories_Quad(T, Y, Number_of_Strides, ax(2).Position, ax(2));
        graphicsObjs.TuglineForce = SLIP_TuglineForce( TF_Sim, lf_exp, Number_of_Strides, ax(3).Position, ax(3));
        pbaspect(graphicsObjs.TuglineForce.ax, [2*Number_of_Strides 1 1]);
        
        graphicsObjs.Animation = SLIP_Animation_Quad_Load(Y, P(1,:), ax(4).Position, ax(4), options);

        graphicsObjs.Sensitivity = SLIP_Sensitivity_Quad(sensitivity_data, [ax(5); ax(6)]);
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

        if t < P(1,9)
            P_curr = P(1,:);
        else
            P_curr = P(2,:);
            P_curr(1:9) = P_curr(1:9) + P(1,9);
        end

        if ~isempty(graphicsObjs.Animation)
            graphicsObjs.Animation.Para = P_curr;
            graphicsObjs.Animation.update(t, y, P_curr); 
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
        for f = 1:length(fields)
            if ~isfield(options, fields{f}) || isempty(options.(fields{f}))
                options.(fields{f}) = defaults.(fields{f});
            end
        end
    end
end
