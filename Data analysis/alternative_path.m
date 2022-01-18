function newpath_ = alternative_path(fullpath)
% Set the path to find here:
% This is helpful when you send the general data files around

% The program will try to find the movies' data in these alternative folder,
% If not found, move on the next one.
% Make sure the global raw data folder's name is 'RAW'
    force_change = true;    % Force to change the folder
    
    if exist(fullpath,'dir')
        if ~force_change
            newpath_ = fullpath;
            return;
        end
    end
    newpath = ...
        {'D:\Users\HuyT\Dropox_curie\Dropbox (UMR3664)\Data\RAW',...    % Huy home computer
        'E:\Dropbox Curie\Dropbox (UMR3664)\Data\RAW',...               % Huy work computer
        'C:\Users\bvhutr\Dropbox (UMR3664)\Data\RAW',...                % Huy TUNI computer
        'C:\Users\bvhutr\Dropbox (UMR3664)\backup\Memory_analysis\Mounia data\RAW', ... Mounia data
        };
    newpath_ = fullpath;
    for i=1:numel(newpath)
        tmppos = findstr(fullpath,'RAW');
        tmpname = fullpath(tmppos+4:end);
        newpath_ = fullfile(newpath{i},tmpname());
        if exist(newpath_,'dir')
            return;
        end
    end