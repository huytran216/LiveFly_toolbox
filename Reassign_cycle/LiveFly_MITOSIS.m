function []=LiveFly_MITOSIS(mov_folder)
    if nargin==0
        mov_folder='\\isiserver.curie.net\umr3664\equipe_dostatni\g_fernandes\RAW\190424_GF2\table_summary\';               %indicate the full path movies file
    end

    if strcmp(questdlg('Choose a configuration file?'),'Yes')
        % If there is only one config file
        listing = dir(mov_folder);
        ConfigName='';
        PathName=mov_folder;
        for i=1:numel(listing)
            if strfind(listing(i).name,'config.m')
                if numel(ConfigName)
                    ConfigName='';
                    break;
                else
                    ConfigName=listing(i).name;                    
                end
            end
        end
        % Reuse config file if needed
        if ~numel(ConfigName)
            [ConfigName,PathName,~] = uigetfile(fullfile(mov_folder,'*.m'),'Select the configuration file');
        end
        ConfigName
        run(fullfile(PathName,ConfigName));
        filename=[ConfigName(1:end-9) '.txt'];
    else
        filename='th26_15.txt'; % Change file name here if params are set manually
        dt=13.5;                % Time step in seconds that you read from the metafile
        A_pole=0;               % Anterior pole position 1=left
        channel=0;              % Channel for Nuclei image (0 or 1)
        x_resolution=0.1968;    % Pixel size (micron/pixels) on XY axis (optional)
        z_resolution=0.5;       % Pixel size (micron/pixels) on Z axis
    end

    nuc=fullfile(mov_folder,filename);
    %% Fixing parameters
        % Cell cycle to correct:
        cycle_range=[10 11 12 13];
        % Threshold to detect fast/slow drifting speed
        drift_thresh=[0.07 0.07];
        % Maximum patching time before and after a trace (in seconds):
        % Default: [50 30 70]
        patch_before=150;   % Roughly equal first spot appearance (in second)
        patch_after=150;    % Roughly equal last spot till mitosis (in second)
        min_percent=70;     % Minimum trace length/interphase duration to consider for correction (in %). Must be greater than 50.
    %% Run the analysis
    reassign_cycle(nuc,dt,x_resolution,drift_thresh,[patch_before patch_after min_percent],cycle_range);