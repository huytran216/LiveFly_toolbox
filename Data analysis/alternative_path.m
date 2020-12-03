function newpath_ = alternative_path(fullpath)
% Set the path to find here:
% This is helpful when you send the general data files around

% The program will try to find the movies' data in these alternative folder,
% If not found, move on the next one.
% Make sure the global raw data folder's name is 'RAW'
    newpath = ...
        {'D:\Users\HuyT\Dropox_curie\Dropbox (UMR3664)\Data\RAW',...
        '',...
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