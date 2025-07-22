function Section2_Single_Stride_Replication()
    clc;
    global results selectedIndex T Y P GRF graphicsObjs currentfolder isSimulating;
    results = [];
    selectedIndex = 1;
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
        'Position', [0.03*screenSize(4), topY, buttonWidth, buttonHeight], ...
        'ButtonPushedFcn', @(src, event) SelectFolder());

    % Dropdown
    dropdownWidth = 0.12 * screenSize(4);
    matDropdown = uidropdown(fig, ...
        'Position', [0.03*screenSize(4) + buttonWidth + 0.05*buttonWidth, topY, dropdownWidth, buttonHeight], ...
        'Items', GetMatFiles(), ...
        'ValueChangedFcn', @(src, event) LoadSelectedFile());

    % Slider (no ticks or labels)
    sliderX = 0.08*screenSize(4) + buttonWidth + dropdownWidth + 4.0*buttonWidth;
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
    halfHeight = PlotSize(2,2)/2;
    PlotPositions = [
        (1/2)*figSize(3)-(2/2 + 1/10)*PlotSize(1,1),            (0.5/10)*figSize(4)+PlotSize(2,2) + buttonHeight, PlotSize(1,1), PlotSize(1,2);
        (1/2)*figSize(3)-(0/2 - 1/10)*PlotSize(1,1),            (0.5/10)*figSize(4)+PlotSize(2,2) + buttonHeight, PlotSize(1,1), PlotSize(1,2);
        (1/2)*figSize(3)-(3/2)*PlotSize(2,1) - 0.5*buttonWidth, (0.5/10)*figSize(4), PlotSize(2,1),               PlotSize(2,2);
        (1/2)*figSize(3)-(1/2)*PlotSize(2,1),                   (0.5/10)*figSize(4)+halfHeight,                   PlotSize(2,1), halfHeight;
        (1/2)*figSize(3)-(1/2)*PlotSize(2,1),                   (0.5/10)*figSize(4),                              PlotSize(2,1), halfHeight;
        (1/2)*figSize(3)+(1/2)*PlotSize(2,1) + 0.5*buttonWidth, (0.5/10)*figSize(4),                              PlotSize(2,1), PlotSize(2,2)];
    for axidx = 1:6
        ax(axidx) = uiaxes(fig);
        ax(axidx).Position = PlotPositions(axidx, :);
    end

    %% File dropdown
    initialFiles = GetMatFiles();
    matDropdown.Items = initialFiles;
    if ~strcmp(initialFiles{1}, '<none>')
        matDropdown.Value = initialFiles{1};
        LoadSelectedFile();
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

    function files = GetMatFiles()
        matFiles = dir(fullfile(currentfolder, '*.mat'));
        validFiles = {};
        for i = 1:length(matFiles)
            s = load(fullfile(currentfolder, matFiles(i).name));
            if isfield(s, 'results')
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
                results = [];
                indexInput.Value = 1;
                numSolutionsLabel.Text = ' / Nah';
            end
            figure(fig);
        end
    end

    function LoadSelectedFile()
        file = matDropdown.Value;
        if strcmp(file, '<none>'), return; end
        data = load(fullfile(currentfolder, file));
        if ~isfield(data, 'results')
            uialert(fig, 'Selected file does not contain "results".', 'Error');
            return;
        end
        results = data.results;
        indexInput.Value = 1;
        numSolutionsLabel.Text = [' / ', num2str(size(results,2))];
        InitializePlots();
    end

    function InitializePlots()
        if isempty(results), return; end
        val = indexInput.Value;
        if isnan(val) || floor(val) ~= val
            uialert(fig, 'Unable to select solution: Integer input requested.', 'Error');
            return;
        end
        if val < 1 || val > size(results,2)
            uialert(fig, 'Unable to select solution: Index number out of range.', 'Error');
            return;
        end
        selectedIndex = val;

        clc

        if isempty(results), return; end
        X = results(1:22, selectedIndex);
        Para = results(23:end, selectedIndex);

        [~, T, Y, P, GRF, ~] = Quadrupedal_ZeroFun_v2(X, Para);
        slider.Value = 0;
        normTimeInput.Value = 0;

        [gait, abbr, color_plot, linetype] = Gait_Identification(results(:, selectedIndex));

        am = 'Detailed';
        options = struct('AnimationMode', am);
        options = DefaultOptions(options);

        graphicsObjs.Animation = SLIP_Animation_Quad(P, ax(1).Position, ax(1), options);
        graphicsObjs.Orbit = SLIP_PeriodicOrbit_Quad(Y, ax(2).Position, ax(2), color_plot);
        graphicsObjs.Trajectories = SLIP_Trajectories_Quad(T, Y, [ax(3); ax(4); ax(5)]);
        graphicsObjs.GRF = SLIP_GRF_Quad(T, GRF, ax(6));
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
        GRF_ = GRF(T <= t, :);

        if ~isempty(graphicsObjs.Animation)
            graphicsObjs.Animation.update(t, y, P); 
        end
        % if ~isempty(graphicsObjs.Orbit)
        %     graphicsObjs.Orbit.update(y); 
        % end
        % if ~isempty(graphicsObjs.Trajectories)
        %     graphicsObjs.Trajectories.update(T_, Y_); 
        % end
        % if ~isempty(graphicsObjs.GRF)
        %     graphicsObjs.GRF.update(T_, GRF_); 
        % end
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
