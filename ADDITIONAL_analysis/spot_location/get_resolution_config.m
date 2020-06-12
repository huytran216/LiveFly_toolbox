function [x_resolution,z_resolution]=get_resolution_config(filename)

    % Get initiated text
    A_ = cell(1, 4);
    A_{1} = 'Config file created on 19-Jul-2017';
    A_{2} = 'THIS IS THE CONFIGURATION FILE FOR SPOT DETECTION PROGRAM: ';
    A_{3} = 'If you need to reanalyze the data, change the params here';
    A_{4} = 'Then rerun the Spot detection with this file loaded';

    A_pole = 1;

    averaging_radius = 3;

    channel = 0;

    dt = 15.4;

    fact_r = 1.2;

    main_mov = 'B6H6HisMCPnoNLS1.lsm';

    nuclei_mov = 'RED_B6H6HisMCPnoNLS1_MAX';

    shift_left = 407;

    shift_right = 805;

    th = [25 15];

    th1 = 25;

    th2 = 15;

    voxels_max = 60;

    voxels_min = 10;

    x_resolution = 0.1968;

    z_resolution = 0.5;
    
    % Get folder name
    if strfind(filename,'_.')
        filename = filename(1:end-5);
    else
        filename = filename(1:end-4);
    end
    filename = filename(find(filename=='_',1,'first')+1:end);
    foldername = ['Z:\equipe_dostatni\c_perez-romero\RAW\' filename '\table_summary'];
    listing = dir(foldername);
    % find the config file
    for i=1:numel(listing)
        if strfind(listing(i).name,'_config.m')
            run(fullfile(foldername,listing(i).name));
        end
    end    
end
    