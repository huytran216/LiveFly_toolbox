%% Re-extract all features in all strains before running the script
% THe data should be in Data_analysis/
warning off;
%% Params
fld='../../Data analysis/';
load([fld 'feature_label.mat']);

% Number of samples - set to 0 when studying t0
Nsample_min=5;      % Minimum total nuclei per position
Nsample_indi=5; % Minimum nuclei per embryo per position

plot_embryo_error=1;    % plot error based on embryo diversity (1) or nuclei error (merged, 0)
shaded_error_bar = 1;   % Plot shaded errorbar or normal errorbar

kymo_intensity = 1;     % Plot kymograph of loci intensity (1) or Pspot (0) or ON (2)
compare_1x2x = false;

smooth_curve = 1;
%% Feature to plot, plot settings
fea_range=[21];
istrimed_range = [1]; % Applied for trimmed traces? (size = fea_range);

nc_range=[13];

AP_limit = [-30 20];
is_compare_experiments = true;   % Compare between experiments
is_compare_features = false;     % Compare between experiments, fea_range is vector, nc_range and compare_list usually scalar
is_make_kymo = false;             % Make kymograph
    fit_boundary= true;         % Find boundary or not?
    fit_sigmoid = true;         % Fit with a sigmoid curve or find half position (if fit_boundary=true)
    plot_vertically = false;     % Plot reporter vertically?
    plot_boundary = false;
%% Set up data list

dtset = struct('filename','','label','');
dtset(1).filename = 'hb-vk33';  dtset(1).label = 'hb-P2';

dtset(2).filename = 'B6-near';  dtset(2).label = 'B6';
dtset(3).filename = 'B9-near';  dtset(3).label = 'B9';
dtset(4).filename = 'B12-near'; dtset(4).label = 'B12';
dtset(5).filename = 'B6-far';   dtset(5).label = 'B6-far';
dtset(6).filename = 'B9-far';   dtset(6).label = 'B9-far';
dtset(7).filename = 'B12-far';  dtset(7).label = 'B12-far';

dtset(8).filename = 'hb-II';  dtset(8).label = 'rand. II';
dtset(9).filename = 'hb-III-Lucas2018';  dtset(9).label = 'rand. III';

dtset(10).filename = 'H6B6-near';   dtset(10).label = 'H6B6';

dtset(11).filename = 'Z6';  dtset(11).label = 'Z6';
dtset(12).filename = 'Z2B6-near';  dtset(12).label = 'Z2B6';
dtset(13).filename = 'Z7B6-near';  dtset(13).label = 'Z7B6';

dtset(14).filename = 'B6near+hbprom';  dtset(14).label = 'B6-hb-P2';

%compare_list = [1 2 5 3 7];isBcd1X = [0 0 0 0 0]; % For B6-B9-B12 comparison
compare_list = [1 2 3 4]; isBcd1X = [0 0 0 0]; % For hb-B6-H6B6 comparison
%compare_list = [1 2 3 7 1 2 3 7];isBcd1X=[0 0 0 0 1 1 1 1]; % For hb-B6-H6B6 comparison, 1x2x
%compare_list = [1 8 9]; isBcd1X =[0 0 0 ];% For vk33 vs random insertion
%compare_list = [1 2 3 4];isBcd1X=[0 0 0 0];
%compare_list = [2 4 10 12]; isBcd1X = compare_list*0;


avr = [600 750 1100];                % Mean nc duration
avr_cut = [450 500 1000];           % Cut time window (until mitosis)   
%% Set folder containing mean data (contain dash)
folder={};
folder{1}='tmp/';
folder{2}='tmp_trimmed/';
%% Cook label_list
DatasetLabel = {dtset(compare_list).label};
DatasetFile = {dtset(compare_list).filename};
for i=1:numel(compare_list)
    if isBcd1X(i)==1
        DatasetLabel{i}=[DatasetLabel{i} '-Bcd1X'];
        DatasetFile{i}=[DatasetFile{i} '-Bcd1X'];
    end
    if isBcd1X(i)==2
        DatasetLabel{i}=[DatasetLabel{i} '-dBcd'];
        DatasetFile{i}=[DatasetFile{i} '-dBcd'];
    end
end
%% Plots
if is_compare_experiments
    compare_experiments;
end
if is_compare_features
    compare_features;
end
if is_make_kymo
    show_kymo;
end