% Add current script folder + all subfolders to MATLAB path

thisFile = matlab.desktop.editor.getActiveFilename();
if isempty(thisFile)
    error('No active script detected. Open this script in the MATLAB Editor and run it.');
end

thisDir = fileparts(thisFile);
addpath(genpath(thisDir));

