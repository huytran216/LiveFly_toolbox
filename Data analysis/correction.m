% Delimiter for the result file (default=' ')
delimiter=' ';

% 3D Image size: X is WIDTH; Y IS HEIGHT
xlen=0;
ylen=0;

% Cropping the frame: nuclei outside this frame will not be considered.
xlim_left=0;
xlim_right=0;
ylim_up=0;
ylim_down=0;

% Intensity column (default = 24 for thresholding method, 28 for gaussian method)
Icolumn=28;

% Frames to remove due to protein aggregate. 
% All spots detected in this frame are ignored
rm={};

%% Nuclei to ignore: enter the ID of the nuclei to ignore
nclist={66};


% In case we use a hybrid threshold, we have different result file
% specify range1 and range2 the frame range for differnt result file
% In case we only have one result file, set range1 to 1:20000 (all movie)
range1=1:10001;filename1='Result_file1';
range2=10002:30000;filename2='Result_file2';

%Cycle 13 maybe screwed up with agregates. We can manually add border position
%and see which nuclei can be kept.
border13=10000;
keepid13=[];

% Is the movie right or left censored, right censored or not?
% These variable are set to zero if the movie begin or end with mitosis.
% They are set to 1 if otherwise. 
% Default: 0
leftcensored=0;
rightcensored=0;


% Ignore frames that are right before/mitosis (from cycle 9, 10, 11, 12, 13, 14)
% This is helpful as protein aggregate usually appears at the end of cc13.
% Ignore this otherwise.
startcutcycle=[0 0 0 0 0 0];
endcutcycle=[10000 10000 10000 10000 10000 10000];