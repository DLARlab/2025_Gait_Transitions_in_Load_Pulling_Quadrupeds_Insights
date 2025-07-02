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

    ax1 = uiaxes(fig, 'Position', [0.05*figWidth, 0.05*figHeight, 0.4* figHeight, 0.3* figHeight], ...
        'Box', 'on');
    ax2 = uiaxes(fig, 'Position', [0.55*figWidth, 0.05*figHeight, 0.4* figHeight, 0.3*figHeight], ...
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
        ExtractFootfall(raw_data);
        % Optional: call plot or update routines like InitializePlots();
    end

    %% Function to extract footfall sequences (call after loading data)
    function ExtractFootfall(raw_data)
        footfall_data = raw_data(:, {'LF', 'LH', 'RF', 'RH'});
        % Display or use footfall_data as needed
    end
end
