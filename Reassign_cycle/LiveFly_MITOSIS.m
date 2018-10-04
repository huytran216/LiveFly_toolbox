mov_folder='\\zserver.curie.fr\umr3664\equipe_dostatni\c_perez-romero\RAW\170405_CP1\table_summary';               %indicate the full path movies file
filename='th45_20.txt';

if strcmp(questdlg('Choose a configuration file?'),'Yes')
    [ConfigName,PathName,FilterIndex] = uigetfile(fullfile(mov_folder,'*.m'),'Select the configuration file');
    run(fullfile(PathName,ConfigName));
else
    dt=13.5;                % Time step in seconds that you read from the metafile
    A_pole=0;               % Anterior pole position 1=left
    channel=0;              % Channel for Nuclei image (0 or 1)
    x_resolution=0.1968;    % Pixel size (micron/pixels) on XY axis (optional)
    z_resolution=0.5;       % Pixel size (micron/pixels) on Z axis
end

nuc=fullfile(mov_folder,filename);
%% Fixing parameters
    % Cell cycle to correct:
    cycle_range=[10 11 12 13 14];
    % Threshold to detect fast/slow drifting speed
    drift_thresh=[0.07 0.07];
    % Maximum patching time before and after a trace (in seconds):
    % Default: [50 30 70]
    patch_before=100;   % Roughly equal first spot appearance (in second)
    patch_after=100;    % Roughly equal last spot till mitosis (in second)
    min_percent=70;     % Minimum trace length/interphase duration to consider for correction (in %). Must be greater than 50.
%% Run the analysis
reassign_cycle(nuc,dt,x_resolution,drift_thresh,[patch_before patch_after min_percent],[10:14]);