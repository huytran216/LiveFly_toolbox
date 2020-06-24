function LiveFly_VISUALIZER
debug=1;
%% Author
% Huy Tran
%% ideas
% Interface to manage all parameters
%% Parameters for positions and colors of figures and some handles
% The position of the game figure will be [figX1,figX2,figY1,figY2]
figX1 = 0;
figX2 = 1200;
figY1 = 0;
figY2 = 700;

feature_label={};feature_unit={};
Nfea=17;
load('feature_label.mat');
%% sets initial variables
% Dataset path
    AP=[-20 20];
    DatasetName = '';       % Dataset name
    DatasetPath = 'data';       % Dataset save path
    Nmov=0;                 % Number of movie
% Storage movie list
    emptystruct=orderfields(struct('tscnt',[],'Path','','Name','','correction',false, ...
        'nc9',false,'nc10',false,'nc11',false,'nc12',false,'nc13',false,'nc14',false, ...
        'shift_EL',0,'Nf',0,'dt',0,'Lx',0,'Ly',0,'BG',1,'BG_man',1,...
        'ShiftL',0,'ShiftR',0,'APole',1 ...
        ));
    DatasetList = emptystruct;
% Storage ensemble dataset
    datamat=struct('id',[],'Intensity',{},'Feature',[],'x',[],'y',[],'border',[],'mother',[],'daughter1',[],'daughter2',[],'tscnt',[],'cycle',[],'tlen',[],'start',[]);

% Storage for extracted features
    DatasetFeature=struct('idrec',[],'tsrec',[],'xrec',[],'yrec',[],'fearec',{},'xborder_rec',[],'hborder_rec',[],'wborder_rec',[],'vborder_rec',[],'xaxis_all',[],'yaxis_all',[],'fearec_all',{});
    FitResult = cell(1,2);
    
% Record for the heatmap
    heatmapI=struct;   % Heat map with absolute time (unscaled interphase)
    
% Record the mean curves for merged embryos {feature x cycle} x [position]
    mf_rec={};  % Mean curve
    sf_rec={};  % Standard deviation
    nf_rec={};  % Number of nuclei
    ef_rec={};  % Number of embryo
    
% Record the mean curves for individual embryos {feature x cycle x embryo} x [position]
    mf_indi={};  % Mean curve
    sf_indi={};  % Standard deviation
    nf_indi={};  % Number of nuclei

% Cell cycle range:
    nc_range=[9:14];                                % Valid nuclear cycle
    fea_ref=0;                                      % Reference feature for embryo alignment
    binwidth=3;                                     % Bin width
    tinterphase_min = [10 10 50 100 200 300];       % minimum interphase duration
    tinterphase_max = [1e4 1e4 1e4 1e4 1e4 1e4];    % maximum interphase duration
    fea_int = [8 9];                                % Which feature is intensity
    trimmed = false;                                % Trimmed traces or not (default = false)
    

Nsample_indi_min = 2;   % Number of samples required per bin for individual embryo
Nsample_all_min = 2;   % Number of samples required per bin for merged embryos
%% Create and hide the GUI figure as it is being constructed.
segfigure = figure('Visible','on','Tag','segfigure','Position',[figX1,figY1,figX2,figY2]);
set ( gcf, 'Color', [0 0 0] )

%% Create menu
set(gcf,'Menubar','None');
% Image or mask file loading, saving:
mFile = uimenu('Label','File_');

    hcreatedataset = uimenu(mFile,'Label','Create new dataset',...
    'Callback',@hcreatedataset_Callback); % load file

    hloadbutton = uimenu(mFile,'Label','Load dataset',...
    'Callback',@hloaddataset_Callback);
    
    hloadbutton = uimenu(mFile,'Label','Quick process...',...
    'Callback',@hquickprocess_Callback);

    hsavebutton = uimenu(mFile,'Label','Save dataset',...
    'Callback',@hsavedataset_Callback);

    hsavebutton = uimenu(mFile,'Label','Save as...',...
    'Callback',@hsaveasdataset_Callback);

    hexportxls = uimenu(mFile,'Label','Export xls...',...
        'Callback',@hexportxls_Callback);
    
    hquitbutton = uimenu(mFile,'Label','Quit',...
    'Callback',@hquit_Callback); % Quit seg

mHelp = uimenu('Label','Help_');

    hshow_feature_meaning = uimenu(mHelp,'Label','Show feature description',...
        'Callback',@show_feature_meaning);
%% Panel for movie info
hp = uipanel('Parent',gcf,'Title','Movie list','FontSize',12,...
             'BackgroundColor','white','Units','pixel',...
             'Position',[20 100 700 500]...
             );
         htable =   uitable('Parent',hp,'Position',[0 0 700 460],...
             'CellSelectionCallback',@htableselection_Callback,...
             'CellEditCallback',@htableedit_Callback);
         htable.ColumnName = {'Color','Name','shift_EL',...
             'nc9','nc10','nc11','nc12','nc13','nc14','BG','BG_man'};
         set(htable,'ColumnEditable',logical([0 0 1 1 1 1 1 1 1 0 1]));
         set(htable,'ColumnWidth',{60 100 60 40 40 40 40 40 40 60 60});
         % Currently selected row
        crrrow=0;
%% Create buttons and others

% Message: telling which dataset is loaded
    hmessages = uicontrol('Style','text','String','Create or load a dataset',...
        'Position',[50,600,510,90],'FontSize',14, 'HorizontalAlignment','left', ...
        'ForegroundColor','white','BackgroundColor','black'); % Display instructions

% Manipulate movie list
    haddmovie = uicontrol('Style','pushbutton','String', 'Add...',...
        'Position',[20,20,90,20],...
        'Callback',@haddmovie_Callback);

    hremovemovie = uicontrol('Style','pushbutton','String', 'Remove',...
        'Position',[120,20,90,20],...
        'Callback',@hremovemovie_Callback);

    hmoveup = uicontrol('Style','pushbutton','String', 'Move Up',...
        'Position',[220,20,90,20],...
        'Callback',@hupmovie_Callback);

    hmovedown = uicontrol('Style','pushbutton','String', 'Move Down',...
        'Position',[320,20,90,20],...
        'Callback',@hdownmovie_Callback);

    hbrowse = uicontrol('Style','pushbutton','String', 'Browse',...
        'Position',[420,20,90,20],...
        'Callback',@hbrowse_Callback);
    
    hbrowse = uicontrol('Style','pushbutton','String', 'Copy path',...
        'Position',[520,20,90,20],...
        'Callback',@hcopypath_Callback);

% Set AP axis:
    hmessages_AP = uicontrol('Style','text','String','AP axis range:',...
        'Position',[730,650,200,30],'FontSize',10, 'HorizontalAlignment','left',...
        'ForegroundColor','white','BackgroundColor','black');
    hmessages_APfrom = uicontrol('Style','text','String','A(%):',...
        'Position',[730,640,50,20],'FontSize',10, 'HorizontalAlignment','left',...
        'ForegroundColor','white','BackgroundColor','black');
    hmessages_APto = uicontrol('Style','text','String','P(%):',...
        'Position',[830,640,50,20],'FontSize',10, 'HorizontalAlignment','left',...
        'ForegroundColor','white','BackgroundColor','black');

    hAPfrom = uicontrol('Style','edit','String',num2str(AP(1)),...
        'Position',[770,640,50,20],'FontSize',10,'Callback',@hsetAP_Callback);
    hAPto = uicontrol('Style','edit','String',num2str(AP(2)),...
        'Position',[870,640,50,20],'FontSize',10,'Callback',@hsetAP_Callback);

% Movie alignment:
    hmessages_fea_ref = uicontrol('Style','text','String','Align by Feature:',...
        'Position',[730,600,200,20],'FontSize',10, 'HorizontalAlignment','left',...
        'ForegroundColor','white','BackgroundColor','black');
    hfea_ref = uicontrol('Style','popup','String',{'None','Manual'},...
        'Position',[850,600,70,20],'FontSize',10, 'HorizontalAlignment','left',...
        'Callback',@hset_fea_ref);

% Set bin width
    hmessages_binwdith = uicontrol('Style','text','String','Bin width (%EL):',...
        'Position',[730,570,200,20],'FontSize',10, 'HorizontalAlignment','left',...
        'ForegroundColor','white','BackgroundColor','black');

    hbinwidth = uicontrol('Style','edit','String',num2str(binwidth),...
        'Position',[850,570,70,20],'FontSize',10, 'HorizontalAlignment','left',...
        'Callback',@hset_binwidth);

% Fit option: Option when fitting with a Hill/Sigmoid function
    fitoption_Hill=logical([0 ...                           % 1: if sharing maximum level
        1 ...                                               % 1: if sharing Hill coefficient
        0 ...                                               % 1: if sharing border position
        0 ...                                               % 1: if fit based on mean curves
        ]);
    fitoption_Hill_text={'Similar maximum level', ...
        'Similar Hill coeff','Similar border','Fit mean curve'};
    hmessage_fitoption = uicontrol('Style','text','String','Hill fitting option:',...
        'Position',[950,650,200,30],'FontSize',10, 'HorizontalAlignment','left',...
        'ForegroundColor','white','BackgroundColor','black');
    hfitoption = uitable('Position',[950 570 210 90],...
        'CellEditCallback',@hfitoptionedit_Callback);
    hfitoption.ColumnName = {'Option','Yes?'};
    hfitoption.Data=[fitoption_Hill_text' num2cell(fitoption_Hill')];
    set(hfitoption,'ColumnEditable',logical([0 1]),'ColumnWidth',{130,40});

% Manipulate movie data
    hmessages_dataset = uicontrol('Style','text','String','Dataset manipulation',...
        'Position',[730,520,150,30],'FontSize',10, 'HorizontalAlignment','left',...
        'ForegroundColor','white','BackgroundColor','black');

    hloadmovie = uicontrol('Style','pushbutton','String', 'Load movies',...
        'Position',[730,500,100,30],...
        'Callback',@hloadmovie_Callback);

    hsetinterphase = uicontrol('Style','pushbutton','String', 'Refine interphase',...
        'Position',[850,500,100,30],...
        'Callback',@hsetinterphase_Callback);

    hshowfeature = uicontrol('Style','pushbutton','String', 'Extract feature',...
        'Position',[970,500,100,30],...
        'Callback',@hextractfeature_Callback);

% Manipulate individual movie
    hmessages_individual = uicontrol('Style','text','String','Show single movie feature',...
        'Position',[730,430,200,30],'FontSize',10, 'HorizontalAlignment','left',...
        'ForegroundColor','white','BackgroundColor','black');

    hindi_feature = uicontrol('Style','pushbutton','String', 'Show Feature',...
        'Position',[730,410,100,30],'FontSize',8,...
        'Callback',@hindi_feature_Callback);

    hindi_traces = uicontrol('Style','pushbutton','String', 'Show Traces',...
        'Position',[850,410,100,30],'FontSize',8,...
        'Callback',@hindi_traces_Callback);
    
    hindi_info = uicontrol('Style','pushbutton','String', 'Show Info',...
        'Position',[970,410,100,30],'FontSize',8,...
        'Callback',@hindi_info_Callback);
% Manipulate all movies
    hmessages_all = uicontrol('Style','text','String','Show all movie feature',...
        'Position',[730,330,200,30],'FontSize',10, 'HorizontalAlignment','left',...
        'ForegroundColor','white','BackgroundColor','black');

    hall_feature = uicontrol('Style','pushbutton','String', 'Summary',...
        'Position',[730,310,100,30],...
        'Callback',@hall_feature_Callback);
    
    hall_pattern = uicontrol('Style','pushbutton','String', 'View pattern',...
        'Position',[850,310,100,30],...
        'Callback',@hall_pattern_Callback);
    
    hall_kymo = uicontrol('Style','pushbutton','String', 'View kymograph',...
        'Position',[970,310,100,30],...
        'Callback',@hall_kymo_Callback);
    
    hall_table = uicontrol('Style','pushbutton','String', 'View fit table',...
        'Position',[1090,310,100,30],...
        'Callback',@hall_table_summary_Callback);
    
    hall_slide = uicontrol('Style','pushbutton','String', 'TimeMagnifier...',...
        'Position',[730,270,100,30],...
        'Callback',@hall_TimeMaginifier_Callback);
    
% Interested feature:
    hmessages_feature = uicontrol('Style','text','String','Interested feature:',...
        'Position',[730,230,200,30],'FontSize',10, 'HorizontalAlignment','left',...
        'ForegroundColor','white','BackgroundColor','black');
    hfeatable = uitable('Position',[730 10 180 220],...
        'CellEditCallback',@hfeatableedit_Callback);
    hfeatable.ColumnName = {'Feature','Show?'};
    fea_plot={};
    for sthtmp=1:Nfea
        fea_plot{sthtmp}=false;
    end
    set(hfeatable,'ColumnEditable',logical([0 1]),'ColumnWidth',{80,40});
    hfeatable.Data=[feature_label(:),fea_plot(:)];
    fea_plot=[];
% Way to show data:
    hmessages_showoption = uicontrol('Style','text','String','Options:',...
        'Position',[930,230,200,30],'FontSize',10, 'HorizontalAlignment','left',...
        'ForegroundColor','white','BackgroundColor','black');
    hshowoption = uitable('Position',[930 10 230 220],...
        'CellEditCallback',@hshowoptionedit_Callback);
    hshowoption.ColumnName = {'Option','Yes?'};
    showoption_list={'Save kymograph movies',...
        'Show fitted Hill curve',...
        'Separate embryos by color',...
        'Show embryo legend',...
        'Show merged mean curve',...
        'Show indi mean curve',...
        'Normalize fitted feature',...
        'Normalize by intensity',...
        'SError instead of SDev',...
        };
    showoption_plot={};
    for sthtmp=1:numel(showoption_list)
        showoption_plot{sthtmp}=false;
    end
    set(hshowoption,'ColumnEditable',logical([0 1]),'ColumnWidth',{150,40});
    hshowoption.Data=[showoption_list(:),showoption_plot(:)];
    showoption_plot=[];
    for sthtmp=1:numel(showoption_list)
        showoption_plot(sthtmp)=false;
    end
%% Final settings
% Assign the GUI a name to appear in the window title.
set(segfigure,'Name','GETDATA (Press F1 for help)') 
% Move the GUI to the center of the screen.
movegui(segfigure,'center')
% Make the GUI visible.
set(segfigure,'Visible','on');
%% Manage dataset
    function hquit_Callback(~,~)
        close(segfigure);
    end
    
    % Create new dataset
    function hcreatedataset_Callback(~,~)
        try 
            % Ask for input file name
            DatasetName_ = inputdlg('Enter a valid new dataset name');
            if isempty(regexp(DatasetName_{1}, '[/\*:?"<>|]', 'once'))
                DatasetName = DatasetName_{1};
                DatasetPath = 'data/';  % Default save location
                Nmov=0;

                % Storage movie list
                DatasetList = emptystruct;
                % Storage ensemble dataset
                datamat=struct('id',[],'Intensity',{},'Feature',[],'x',[],'y',[],'border',[],'mother',[],'daughter1',[],'daughter2',[],'tscnt',[],'cycle',[],'tlen',[],'start',[]);
                
                % Refresh the table
                Update_List();
             else
                 msgbox('Error/Invalid file name');
             end
        catch exception
            if debug
                throw(exception);
            else
                msgbox('Error/Invalid file name');
            end        
         end
    end

    % Load/Save dataset
    function hloaddataset_Callback(~,~)
        [FileName_,PathName_,~] = uigetfile('*.mat','Select the  file');
        try
            if FileName_
                DatasetName=FileName_;
                DatasetPath=PathName_;
                f=waitbar(0,'Loading');
                load(fullfile(PathName_,FileName_),'DatasetList','datamat','DatasetFeature','Nmov','AP','fea_ref','binwidth','tinterphase_min','tinterphase_max','fitoption_Hill','heatmapI','mf_rec','sf_rec','nf_rec','ef_rec','mf_rec','sf_rec','nf_rec','ef_rec','mf_indi','sf_indi','nf_indi');
                fitoption_Hill(end+1:4)=0;
                close(f);
                % Patch DatasetList if some new field is missing
                if ~isfield(DatasetList,'shift_EL')   % Reference feature
                    for i=1:numel(DatasetList)
                        DatasetList(i).shift_EL=0;
                    end
                end
                if ~isfield(DatasetList,'BG')
                    for i=1:numel(DatasetList)      % Spot Background intensity
                        DatasetList(i).BG=1;
                    end
                end
                if ~isfield(DatasetList,'BG_man')   % Background intensity (manually added)
                    for i=1:numel(DatasetList)
                        DatasetList(i).BG_man=1;
                    end
                end
                % Remove unnecessary field
                tmp1=fieldnames(DatasetList);
                for i=1:numel(tmp1)
                    if ~isfield(emptystruct,tmp1{i})
                        DatasetList=rmfield(DatasetList,tmp1{i});
                    end
                end
                % Order DatasetList
                DatasetList=orderfields(DatasetList);
                Update_List();
            end
        catch exception
            if debug
                throw(exception);
            else
                msgbox('Error loading Dataset');
            end            
        end
    end

    % Load/Save dataset
    function hquickprocess_Callback(~,~)
        % Imput for bulk process
        Outtext=cell(1,numel(nc_range));
        Deftext=cell(1,numel(nc_range));
        if numel(datamat)
            for i=1:numel(nc_range)
                Outtext{i}=['Trim: nc' num2str(nc_range(i)) ' (from to (in second)) ' ];
                if nc_range(i)==13
                    Deftext{i}=['600 750'];
                else
                    Deftext{i}=['0 10000'];
                end
            end
            dlg=inputdlg(Outtext,'Set the parameters for quick processing',[1 50],Deftext);
            for i=1:numel(nc_range)
                tmp=str2num(dlg{i});
                tcut1(i)=tmp(1);
                tcut2(i)=tmp(2);
                if tcut1>=tcut2
                    msgbox('Invalid input');
                    return;
                end
            end
        else
            msgbox('Load movie first');
            return;
        end
        % QUICK PROCESSING
            % Analyze untrimmed trace
            trimmed=false;
                hextractfeature_Callback();
                hall_kymo_Callback();
            % Open magnifier
                h=figure;
                %try
                    [tlower,tupper,cycle_range,posborder]=Magnifier(h,heatmapI,binwidth,nc_range(2:end),tcut1(2:end),tcut2(2:end));
                    for i=1:numel(cycle_range)
                        heatmapI(cycle_range(i)-8).posborder=posborder(1,i);
                    end
                    Re_Extract_feature(cycle_range,tlower,tupper);
                    trimmed=true;
                %catch
                %end
            % Analyze trimmed trace
                hextractfeature_Callback();
    end

    function hsavedataset_Callback(~,~)
        if exist(fullfile(DatasetPath,DatasetName),'file')
            if strcmp(questdlg('Overwrite?','Dataset file exist','Yes','No','No'),'Yes')
                okwrite=true;
            else
                okwrite=false;
            end
        else
            okwrite=true;
        end
        if okwrite
            f=waitbar(0,'Saving');
            save(fullfile(DatasetPath,DatasetName),'DatasetList','datamat','DatasetFeature','Nmov','AP','fea_ref','binwidth','tinterphase_min','tinterphase_max','fitoption_Hill','heatmapI','mf_rec','sf_rec','nf_rec','ef_rec','mf_indi','sf_indi','nf_indi');
            close(f);
        end
    end

    function hsaveasdataset_Callback(~,~)
        [FileName,PathName] = uiputfile('*.mat');
        okwrite=false;
        if FileName
            okwrite=true;
        end
        if okwrite
            DatasetPath=PathName;
            DatasetName=FileName;
            f=waitbar(0,'Saving');
            save(fullfile(DatasetPath,DatasetName),'DatasetList','datamat','DatasetFeature','Nmov','AP','fea_ref','binwidth','tinterphase_min','tinterphase_max','fitoption_Hill','heatmapI','mf_rec','sf_rec','nf_rec','ef_rec','mf_indi','sf_indi','nf_indi');
            close(f);
        end
        Update_List;
    end
    
    function hexportxls_Callback(~,~)
        [FileName,PathName] = uiputfile('DatasetPath/*.xls');
        okwrite=false;
        if FileName
            okwrite=true;
        end
        if okwrite
            f=waitbar(0,'Exporting');
            % Save datamat
                % Compiling outtab:
                outtab=zeros(numel(datamat),24);
                header={};
                outtab(:,1)=[datamat(:).tscnt]';header{1}='TS';
                outtab(:,2)=[datamat(:).id]';header{2}='id';
                outtab(:,3)=[datamat(:).oldid]';header{3}='oldid';
                outtab(:,4)=[datamat(:).cycle]';header{4}='cycle';
                outtab(:,5)=[datamat(:).x]';header{5}='x';
                outtab(:,6)=[datamat(:).y]';header{6}='y';
                outtab(:,7)=[datamat(:).mother]';header{7}='M';
                outtab(:,8)=[datamat(:).daughter1]';header{8}='D1';
                outtab(:,9)=[datamat(:).daughter2]';header{9}='D2';
                waitbar(0.2,f);
                for i=1:numel(feature_label)
                    tmp=arrayfun(@(x) subindex(datamat(x).Feature,i),1:numel(datamat));
                    outtab(:,i+9)=tmp(:);
                end
                header(10:9+numel(feature_label))=feature_label;
                waitbar(0.4,f);                
                outxls=mat2cell(outtab,ones(1,size(outtab,1)),ones(1,size(outtab,2)));
                %csvwrite(fullfile(FilePath,FileName),outtab);
                xlswrite(fullfile(PathName,FileName),[header;outxls],1);
            % Save border feature
                header={};
                header{1,1}='FEATURE:';
                header{2,1}='cycle';
                header{2,2}='TS';
                header{2,3}='xborder';
                header{2,5}='Hill';
                header{2,7}='Width';
                header{2,9}='MaxVal';
                header{2,4}='s_xborder';
                header{2,6}='s_Hill';
                header{2,8}='s_Width';
                header{2,10}='s_MaxVal';
                for feaidx=1:numel(feature_label)
                    header{1,2}=feature_label{feaidx};
                    outtab=zeros(Nmov*numel(nc_range),10)-1;
                    cnt=0;
                    for cycleidx=1:numel(nc_range)
                        if numel(DatasetFeature(cycleidx).xborder_rec)
                            xrec=[];
                            for movidx=1:Nmov
                                if size(DatasetFeature(cycleidx).hborder_rec,2)>=movidx
                                    if DatasetFeature(cycleidx).hborder_rec(feaidx,movidx)
                                        cnt=cnt+1;
                                        outline=[nc_range(cycleidx) movidx DatasetFeature(cycleidx).xborder_rec(feaidx,movidx) 0 DatasetFeature(cycleidx).hborder_rec(feaidx,movidx) 0 DatasetFeature(cycleidx).wborder_rec(feaidx,movidx) 0 DatasetFeature(cycleidx).vborder_rec(feaidx,movidx) 0];
                                        outtab(cnt,:)=outline;
                                        xrec=[xrec;outline([3 5 7 9])];
                                    end
                                end
                            end
                            cnt=cnt+1;
                            outline=[nc_range(cycleidx) 0 mean(xrec(:,1)) sqrt(var(xrec(:,1))) mean(xrec(:,2)) sqrt(var(xrec(:,2))) mean(xrec(:,3)) sqrt(var(xrec(:,3))) mean(xrec(:,4)) sqrt(var(xrec(:,4)))];
                            outtab(cnt,:)=outline;
                        end
                    end
                    outtab=outtab(1:cnt,:);
                    outxls=num2cell(outtab);
                    xlswrite(fullfile(PathName,FileName),[header;outxls],1+feaidx);
                    waitbar(0.4+0.6*feaidx/numel(feature_label),f);
                end
            close(f);
        end
    end
%% Manipulate movie list
    function haddmovie_Callback(~,~)
        if isempty(DatasetName)
            msgbox('Create/Load a dataset first');
        else
            % Load a movie
%            try
                [FileName_,PathName_,~] = uigetfile('correction.m','Select the correction.m file');
                if strcmp(FileName_,'correction.m')
                    Nmov=sum([DatasetList(:).tscnt]>0)+1;
                    % Add new empty field
                    DatasetList(Nmov) =  emptystruct;
                    DatasetList(Nmov).Path=PathName_;
                    DatasetList(Nmov).tscnt=Nmov;
                    % Get initial information
                    Update_DatasetList(Nmov);
                    % Update list
                    Update_List();
                end
%            catch
%                msgbox('No movie loaded');
%            end
        end
    end

    function hupmovie_Callback(~,~)
        % Move current movie up:
        if (crrrow<=Nmov)&&(crrrow>1)
            tmp=DatasetList(crrrow);
            DatasetList(crrrow)=DatasetList(crrrow-1);
            DatasetList(crrrow-1)=tmp;
        end
        Update_List;
    end
    
    function hdownmovie_Callback(~,~)
        % Move current movie down:
        if (crrrow<Nmov)&&(crrrow>0)
            tmp=DatasetList(crrrow);
            DatasetList(crrrow)=DatasetList(crrrow+1);
            DatasetList(crrrow+1)=tmp;
        end
        Update_List;
    end

    function hremovemovie_Callback(~,~)
        % Move current movie down:
        if crrrow
            DatasetList=DatasetList([1:crrrow-1 crrrow+1:end]);
            Nmov=Nmov-1;
            Update_List;
        end
    end
 
    function hbrowse_Callback(~,~)
        if crrrow
            explorer(DatasetList(crrrow).Path);
        end
    end

    function hcopypath_Callback(~,~)
        if crrrow
            clipboard('copy',DatasetList(crrrow).Path);
        end
    end
    
    function htableedit_Callback(~,~)
        datatmp=htable.Data;
        for i=1:Nmov
            DatasetList(i).shift_EL=datatmp{i,3};
            DatasetList(i).nc9=datatmp{i,4};
            DatasetList(i).nc10=datatmp{i,5};
            DatasetList(i).nc11=datatmp{i,6};
            DatasetList(i).nc12=datatmp{i,7};
            DatasetList(i).nc13=datatmp{i,8};
            DatasetList(i).nc14=datatmp{i,9};
            DatasetList(i).BG_man=datatmp{i,11};
        end
    end

    function htableselection_Callback(~,event)
        if numel(event.Indices)
            crrrow=event.Indices(1);
        end
    end

%% Data extraction
    function hloadmovie_Callback(~,~)
        if Nmov
            for i=1:Nmov
                Update_DatasetList(Nmov);
            end
            [datamat,DatasetList]=local_extract(DatasetList);
        else
            msgbox('No movie loaded');
        end
        Update_List;
    end

    function hsetinterphase_Callback(~,~)
        Outtext=cell(1,numel(nc_range));
        if numel(datamat)
            for i=1:numel(nc_range)
                Outtext{i}=['nc' num2str(nc_range(i)) ' (min max) ' ];
                Deftext{i}=[num2str(tinterphase_min(i)) ' ' num2str(tinterphase_max(i))];
            end
            Outtext{i+1}='Automatic (1 for yes, 0 for no)';
            Deftext{i+1}='1';
            dlg=inputdlg(Outtext,'Set the range for interphase duration (s)',[1 50],Deftext);
            for i=1:numel(nc_range)
                tmp=str2num(dlg{i});
                tinterphase_min(i)=tmp(1);
                tinterphase_max(i)=tmp(2);
            end
            tinterphase_min(i+1)=str2num(dlg{i+1});
            trimmed=false;
            refine_tinterphase(tinterphase_min,tinterphase_max);
        else
            msgbox('Load movie first');
        end
    end
    
    function refine_tinterphase(tinterphase_min,tinterphase_max)
        Re_Extract_feature();
        % Restore the feature
        for i=1:numel(datamat)
            datamat(i).Feature=datamat(i).Feature_store;
        end
        % Clean interphase duration as usual from upper lower bounds
            cnt=0;
            for cycleno=nc_range
                cnt=cnt+1;
                idselect=find([datamat(:).cycle]==cycleno);
                tmp1=idselect(arrayfun(@(x) subindex(datamat(x).Feature,10),idselect)<tinterphase_min(cnt));
                tmp2=idselect(arrayfun(@(x) subindex(datamat(x).Feature,10),idselect)>tinterphase_max(cnt));
                for i=[tmp1 tmp2]
                    datamat(i).Feature=-ones(1,numel(feature_label));
                end
            end
        if tinterphase_min(end)  
            % Refine interphase automatically
            cnt=0;
            for cycleno=nc_range
                cnt=cnt+1;
                ts_spec=find(arrayfun(@(x) getfield(DatasetList,{x},['nc' num2str(cycleno)]),1:Nmov));
                for tsidx=ts_spec                    
                    idselect=find(([datamat(:).cycle]==cycleno)&([datamat(:).tscnt]==tsidx));
                    tmp=arrayfun(@(x) subindex(datamat(x).Feature,10),idselect);
                    tinterphase=mode(tmp(tmp>0));
                    tmp1=idselect(arrayfun(@(x) subindex(datamat(x).Feature,10),idselect)<tinterphase*0.70);
                    tmp2=idselect(arrayfun(@(x) subindex(datamat(x).Feature,10),idselect)>tinterphase*1.30);
                    for i=[tmp1 tmp2]
                        datamat(i).Feature=-ones(1,numel(feature_label));
                    end
                end
            end
        end
    end

    function hextractfeature_Callback(~,~)
        f=waitbar(0,'Fitting features');
        normalize_feature=0;
        if numel(datamat)
            cnt=0;
            for cycleno=nc_range
                cnt=cnt+1;
                ts_spec=arrayfun(@(x) getfield(DatasetList,{x},['nc' num2str(cycleno)]),1:Nmov);
                [DatasetFeature(cnt).idrec,DatasetFeature(cnt).tsrec,DatasetFeature(cnt).xrec,DatasetFeature(cnt).yrec,DatasetFeature(cnt).fearec,...
                    DatasetFeature(cnt).xborder_rec,DatasetFeature(cnt).hborder_rec,DatasetFeature(cnt).wborder_rec,DatasetFeature(cnt).vborder_rec,...
                    DatasetFeature(cnt).xborder_CI,DatasetFeature(cnt).hborder_CI,DatasetFeature(cnt).wborder_CI,DatasetFeature(cnt).vborder_CI,...
                    DatasetFeature(cnt).xaxis_all,DatasetFeature(cnt).yaxis_all,DatasetFeature(cnt).fearec_all]=local_extract_feature(datamat,mf_indi,cycleno,find(ts_spec),normalize_feature,feature_label,AP,fitoption_Hill);
            end
        end
        close(f);
        if trimmed
            FitResult{2} = DatasetFeature;
        else
            FitResult{1} = DatasetFeature;
        end
        Update_Mean_Curves;
    end
    
%% Show stuffs from individual movies
    function hindi_feature_Callback(~,~)
        % Summarize selected movies
%         try
            idx=crrrow;
            if idx
                new_nc_range=find(arrayfun(@(x) getfield(DatasetList(crrrow),['nc' num2str(x)]),nc_range));
                if numel(new_nc_range)
                    for fea=fea_plot(:)'
                        figure('Name',['Movie: ' DatasetList(crrrow).Name ', Feature: ' feature_label{fea}]);
                        ymax=0;ymin=0;
                        for cnt=1:numel(new_nc_range)
                            % Show the figure
                            subplot(numel(new_nc_range),2,2*cnt-1);
                            allx=DatasetFeature(new_nc_range(cnt)).xaxis_all{fea,crrrow};
                            ally=DatasetFeature(new_nc_range(cnt)).yaxis_all{fea,crrrow};
                            allf=DatasetFeature(new_nc_range(cnt)).fearec_all{fea,crrrow};
                            h=scatter(allx(:),ally(:),0*allf(:)+20,allf(:),'filled');
                            map = ones(100,3);
                            map(:,1)=map(:,1)-linspace(0,1,100)';
                            map(:,3)=map(:,3)-linspace(0,1,100)';
                            colormap(map);
                            set(h,'MarkerEdgeColor','k');
                            set(gca,'YTick',[]);
                            
                            ylim([1.1*min(ally(:))-0.1*max(ally(:)),1.1*max(ally(:))-0.1*min(ally(:))]);
                            ylabel(['nc' num2str(nc_range(new_nc_range(cnt)))]);
                            xlim(AP);
                            if cnt==1
                                title([DatasetList(crrrow).Name]);
                            end
                            if cnt<numel(new_nc_range)
                                set(gca,'XTick',[]);
                            end
                            
                            set(gca,'color','k');
                            
                            subplot(numel(new_nc_range),2,2*cnt);
                            plot(allx(:),allf(:),'x');
                            xlim(AP);
                            if cnt==1
                                title(feature_label{fea});
                            end
                            if cnt<numel(new_nc_range)
                                set(gca,'XTick',[]);
                            end
                            if showoption_plot(2)
                                hold on;
                                plot([AP(1):AP(2)],2*DatasetFeature(new_nc_range(cnt)).vborder_rec(fea,crrrow)*sigmf([AP(1):AP(2)],[ -DatasetFeature(new_nc_range(cnt)).hborder_rec(fea,crrrow)*0.04 DatasetFeature(new_nc_range(cnt)).xborder_rec(fea,crrrow)]),'--k');
                            end
                            % Adjust the axis
                            ytmp=get(gca,'ylim');
                            ymax=max(ymax,ytmp(2));
                            ymin=min(ymin,ytmp(1));
                        end
                        for cnt=1:numel(new_nc_range)
                            subplot(numel(new_nc_range),2,2*cnt);
                            ylim([ymin-0.2*(ymax-ymin) ymax+0.2*(ymax-ymin)]);
                        end
                        tightfig;
                    end
                end
            end
%         catch exception
%             if debug
%                 throw(exception)
%             else
%                 msgbox('Error. Reanalyze the data');
%             end
%         end
    end
    
    function hindi_traces_Callback(~,~)
        Outtext={'ID from:','ID to','Pos from (%):','Pos to (%)','Only ON trace?:','Maximum Intensity'};
        Deftext={'0','10000','-50','50','0','30'};
        dlg=inputdlg(Outtext,'Show trace value (s)',[1 50],Deftext);
        
        new_nc_range=find(arrayfun(@(x) getfield(DatasetList(crrrow),['nc' num2str(x)]),nc_range));
        ncol=5;
        for cycleno=nc_range(new_nc_range)
            idselect=find(([datamat(:).cycle]==cycleno)&([datamat(:).tscnt]==crrrow)...
                &([datamat(:).oldid]>=str2num(dlg{1})) ...
                &([datamat(:).oldid]<=str2num(dlg{2})) ...
                &([datamat(:).x]*100-50>=str2num(dlg{3})) ...
                &([datamat(:).x]*100-50<=str2num(dlg{4})) ...
                );
            if str2num(dlg{5})
                idselect=idselect(arrayfun(@(x) subindex(datamat(x).Feature,4)>0.1,idselect));
            end
            if numel(idselect)>80
                idselect=idselect(1:80);
                msgbox(['Too more than 80 traces in cc' num2str(cycleno) '. Please narrow the criteria']);
            end
            nrow=ceil(numel(idselect)/ncol);
            figure('Name',['nc' num2str(cycleno)]);
            xmax=0;xmin=0;  % Modifier for time limit (xlim)
            for idx=1:numel(idselect)
                subplot(nrow,ncol,idx);
                ax=([datamat(idselect(idx)).Adjustedtime]-min(datamat(idselect(idx)).Adjustedtime)).*[datamat(idselect(idx)).dt];
                plot(ax,[datamat(idselect(idx)).AdjustedIntensity]);
                set(gca,'YTick',[]);
                set(gca,'XTick',[]);
                ylim([0 str2num(dlg{6})]);
                text(0,2000,{num2str(datamat(idselect(idx)).x*100-50,'%.1f')});
                xmax=max(xmax,max(ax));
                xmin=min(xmin,min(ax));
            end
            for idx=1:numel(idselect)
                subplot(nrow,ncol,idx);
                xlim([xmin xmax]);
            end
        end
    end

    function hindi_info_Callback(~,~)
        % Show info of the movies:
        header={};
        cnt=1;
        header{1,1}='FEATURE:';        
        header{1,2}='VALUE';
        cnt=cnt+1;
            header{cnt,1}='Name';
            header{cnt,2}=DatasetList(crrrow).Name;
        cnt=cnt+1;
            header{cnt,1}='dt';
            header{cnt,2}=DatasetList(crrrow).dt;
        cnt=cnt+1;
            header{cnt,1}='BG intensity';
            header{cnt,2}=DatasetList(crrrow).BG;
        cnt=cnt+1;
            header{cnt,1}='Intensity coeff';
            header{cnt,2}=DatasetList(crrrow).BG_man;
        cnt=cnt+1;
            header{cnt,1}='xlen';
            header{cnt,2}=DatasetList(crrrow).Lx;    
        cnt=cnt+1;
            header{cnt,1}='ylen';
            header{cnt,2}=DatasetList(crrrow).Ly;
        cnt=cnt+1;
            header{cnt,1}='ShiftL';
            header{cnt,2}=DatasetList(crrrow).ShiftL;
        cnt=cnt+1;
            header{cnt,1}='ShiftR';
            header{cnt,2}=DatasetList(crrrow).ShiftR;
        cnt=cnt+1;
            header{cnt,1}='APole';
            header{cnt,2}=DatasetList(crrrow).APole;
        cnt=cnt+1;
            header{cnt,1}='tphase nc9 (s)';
            try
                header{cnt,2}=DatasetFeature(1).vborder_rec(10,crrrow);
            catch
                header{cnt,2}='';
            end
        cnt=cnt+1;
            header{cnt,1}='tphase nc10 (s)';
            try
                header{cnt,2}=DatasetFeature(2).vborder_rec(10,crrrow)*2;
            catch
                header{cnt,2}='';
            end
        cnt=cnt+1;
            header{cnt,1}='tphase nc11 (s)';
            try
                header{cnt,2}=DatasetFeature(3).vborder_rec(10,crrrow)*2;
            catch
                header{cnt,2}='';
            end
        cnt=cnt+1;
            header{cnt,1}='tphase nc12 (s)';
            try
                header{cnt,2}=DatasetFeature(4).vborder_rec(10,crrrow)*2;    
            catch
                header{cnt,2}='';
            end
        cnt=cnt+1;
            header{cnt,1}='tphase nc13 (s)';
            try
                header{cnt,2}=DatasetFeature(5).vborder_rec(10,crrrow)*2;
            catch
                header{cnt,2}='';
            end
        cnt=cnt+1;
            header{cnt,1}='tphase nc14 (s)';
            try
                header{cnt,2}=DatasetFeature(6).vborder_rec(10,crrrow)*2;
            catch
                header{cnt,2}='';
            end
            
            % Put out new figure;
            htmp=figure('Name',['Movie : ' num2str(crrrow) ': ' DatasetList(crrrow).Name ],'Position',[100 100 500 500]);
            htmptable =   uitable('Parent',htmp,'Position',[0 0 500 500]);
            htmptable.Data = header;
    end
%% Show stuffs from all movies        
    function hall_feature_Callback(~,~)
        if numel(datamat)&numel(fea_plot)
            cnt0=0;
            for cycleno=nc_range
                cnt0=cnt0+1;
                ts_spec=find(arrayfun(@(x) getfield(DatasetList,{x},['nc' num2str(cycleno)]),1:Nmov));
                if numel(ts_spec)
                    figure('Name',['Dataset: ' DatasetName '. nc' num2str(cycleno)]);
                    cnt=0;
                    yrange={};
                    ymin=zeros(1,15);ymax=zeros(1,15);
                    for tsidx=ts_spec
                        for feaidx=fea_plot(:)'
                            if numel(DatasetFeature(cnt0).fearec_all{feaidx,tsidx})
                                cnt=cnt+1;
                                xborder=0;
                                % Begin plotting
                                subplot(numel(ts_spec),numel(fea_plot),cnt)
                                plot(DatasetFeature(cnt0).xaxis_all{feaidx,tsidx}-xborder,DatasetFeature(cnt0).fearec_all{feaidx,tsidx},'.b');hold on;
                                tmpx=[AP(1):AP(2)];
                                % Show fitted Hill curve if needed
                                if showoption_plot(2)   % Show fitted Hill coefficient
                                    plot(tmpx,2*DatasetFeature(cnt0).vborder_rec(feaidx,tsidx)*sigmf(tmpx,[ -DatasetFeature(cnt0).hborder_rec(feaidx,tsidx)*0.04 DatasetFeature(cnt0).xborder_rec(feaidx,tsidx)-xborder]),'--k');
                                end
                                % Set axis limit
                                xlim(AP);
                                yrange{feaidx}=[];
                                if tsidx==ts_spec(1)
                                    title(feature_label{feaidx});
                                end
                                if tsidx~=ts_spec(end)
                                    set(gca,'XTick',[]);
                                end
                                if feaidx==fea_plot(1)
                                    ylabel(['embryo ' num2str(tsidx)]);
                                end
                                % Adjust the axis
                                ytmp=get(gca,'ylim');
                                ymax(feaidx)=max(ymax(feaidx),ytmp(2));
                                ymin(feaidx)=min(ymin(feaidx),ytmp(1));
                            end
                        end
                    end
                    cnt=0;
                    for tsidx=ts_spec
                        for feaidx=fea_plot(:)'
                            if numel(DatasetFeature(cnt0).fearec_all{feaidx,tsidx})
                                cnt=cnt+1;
                                subplot(numel(ts_spec),numel(fea_plot),cnt);
                                ylim([ymin(feaidx)-0.1*(ymax(feaidx)-ymin(feaidx)) ymax(feaidx)+0.1*(ymax(feaidx)-ymin(feaidx))]);
                            end
                        end
                    end
                end
            end
        end
    end

    function hall_pattern_Callback(~,~)
        if numel(datamat)   % If data is loaded
            cnt=0;
            pos_range=[AP(1):AP(2)];
            pos_range=[-45:45];
            for cycleno=nc_range
                cnt=cnt+1;
                ts_spec=find(arrayfun(@(x) getfield(DatasetList,{x},['nc' num2str(cycleno)]),1:Nmov));
                if numel(ts_spec)
                    cnt1=0;
                    figure('Name',['Dataset: ' DatasetName '. nc' num2str(cycleno) ]);
                    for fea=fea_plot(:)'
                        cnt1=cnt1+1;
                        subplot(numel(fea_plot),1,cnt1);
                        % Put some alignment here
                        allx=[];allf=[];    % prealignment
                        allx_={};allf_={};  % posalignment
                        cnt2=0;
                        for tsidx=ts_spec
                            cnt2=cnt2+1;
                            xborder=0;
                            allx=[allx DatasetFeature(cnt).xaxis_all{fea,tsidx}-xborder];
                            allx_{cnt2}=DatasetFeature(cnt).xaxis_all{fea,tsidx};
                            allf=[allf DatasetFeature(cnt).fearec_all{fea,tsidx}];
                            allf_{cnt2}=DatasetFeature(cnt).fearec_all{fea,tsidx};
                        end
                        cnt2=0;leg={};
                        if showoption_plot(3)   % Seperate embryo by color
                            for tsidx=ts_spec
                                cnt2=cnt2+1;
                                plot(allx_{cnt2},allf_{cnt2}, ...
                                        'Marker','.','color',corder(tsidx),'LineStyle','none'); hold on;
                                leg{cnt2}=['embryo ' num2str(tsidx)];
                            end
                        else
                            cnt2=cnt2+1;
                            plot(allx,allf,...
                                'Marker','.','LineStyle','none');
                            leg{cnt2}='data';
                        end
                        if showoption_plot(5)   % Show single mean curve
                            hold on;
                            if showoption_plot(9)   % Show standard error instead of standard deviation
                                errorbar(pos_range,mf_rec{fea,cycleno},sf_rec{fea,cycleno}./sqrt(nf_rec{fea,cycleno}-1),'k');
                            else
                                errorbar(pos_range,mf_rec{fea,cycleno},sf_rec{fea,cycleno},'k');
                            end
                            leg{cnt2+1}='Mean';
                        end
                        
                        if showoption_plot(6)   % Show individual mean curves
                            hold on;
                            for tsidx=ts_spec
                                plot(pos_range,mf_indi{fea,cycleno,tsidx},'color',corder(tsidx));
                            end
                        end
                        if showoption_plot(4)   % Show legend
                            legend(leg);
                        end
                        ylabel(feature_label{fea});
                        xlim(AP);
                    end
                end
            end
        end
    end

    function hall_kymo_Callback(~,~)
        if numel(datamat)
            % Calculate the number of run
            cnt_=0;
            for cycleno=nc_range
                ts_spec=find(arrayfun(@(x) getfield(DatasetList,{x},['nc' num2str(cycleno)]),1:Nmov));
                cnt_=cnt_+numel(ts_spec);
            end
            f=waitbar(0,'Making Kymograph...');
            cnt=0;
            for cycleno=nc_range
                ts_spec=find(arrayfun(@(x) getfield(DatasetList,{x},['nc' num2str(cycleno)]),1:Nmov));
                if numel(ts_spec)                    
                    cnt=cnt+numel(ts_spec);
                    % Update progress bar
                    waitbar((cnt-1)/cnt_,f);
                    xborder=zeros(1,Nmov);
                    % Get interphase duration:
                    tphase=DatasetFeature(cycleno-8).vborder_rec(10,:)*2;
                    time_align=0;
                    time_normalize=0;
                    [Imap,Imap_std,Imap_indi,pos_range,time_axis,Cellcount]=local_draw_kymo(datamat,ts_spec,xborder,tphase,cycleno,binwidth,AP,time_normalize,time_align);
                    heatmapI(cycleno-8).Abs_map=Imap;
                    heatmapI(cycleno-8).Abs_map_std=Imap_std;
                    heatmapI(cycleno-8).Abs_map_indi=Imap_indi;
                    heatmapI(cycleno-8).Abs_time=time_axis;
                    heatmapI(cycleno-8).Abs_count=Cellcount;
                    time_align=0;
                    time_normalize=1;
                    if showoption_plot(1)
                        tmpName = DatasetName;
                        Ymax = max(DatasetFeature(cycleno-8).vborder_rec(7,:)*2);
                    else
                        tmpName = '';
                        Ymax = 0;
                    end
                    [Imap,Imap_std,Imap_indi,pos_range,time_axis,mtphase,Cellcount]=local_draw_kymo(datamat,ts_spec,xborder,tphase,cycleno,binwidth,AP,time_normalize,time_align,tmpName,Ymax);
                    heatmapI(cycleno-8).Rel_map=Imap;
                    heatmapI(cycleno-8).Rel_map_std=Imap_std;
                    heatmapI(cycleno-8).Rel_map_indi=Imap_indi;
                    heatmapI(cycleno-8).Rel_time=time_axis;
                    heatmapI(cycleno-8).Rel_count=Cellcount;
                    heatmapI(cycleno-8).mtphase=mtphase;
                    heatmapI(cycleno-8).tphase=tphase;
                    heatmapI(cycleno-8).pos_range=pos_range;
                    % Close all figures
                    all_figs = findobj(0, 'type', 'figure');
                    delete(setdiff(all_figs(2:end), segfigure));
                end
            end
            close(f);
        end
        %msgbox('Done. Press TimeMagnifier to view kymograph');
        trimmed=false;
        mkdir('tmp');
        [~,tmp]=fileparts(DatasetName);
        save(['tmp/' tmp],'pos_range','heatmapI','-append');
    end

    function hall_table_summary_Callback(~,~)
        header={};
        header{1,1}='FEATURE:';
        header{2,1}='cycle';
        header{2,2}='TS';
        header{2,3}='xborder';
        header{2,5}='Hill';
        header{2,7}='Width';
        header{2,9}='HalfMax';
        header{2,4}='s_xborder';
        header{2,6}='s_Hill';
        header{2,8}='s_Width';
        header{2,10}='s_HalfMax';
        if trimmed
            trim_range = [1 2]; % Only show trimmed traces' features
        else
            trim_range = [1];   % Only show whole traces' features
        end
        for istrim = trim_range
            for feaidx=fea_plot(:)'
                header{1,2}=feature_label{feaidx};
                outtab=zeros(Nmov*numel(nc_range),10)-1;
                cnt=0;
                for cycleidx=1:numel(nc_range)
                    if numel(FitResult{istrim}(cycleidx).xborder_rec)
                        xrec=[];
                        for movidx=1:Nmov
                            if size(FitResult{istrim}(cycleidx).hborder_rec,2)>=movidx
                                if FitResult{istrim}(cycleidx).hborder_rec(feaidx,movidx)
                                    cnt=cnt+1;
                                    outline=[nc_range(cycleidx) movidx FitResult{istrim}(cycleidx).xborder_rec(feaidx,movidx) 0 FitResult{istrim}(cycleidx).hborder_rec(feaidx,movidx) 0 FitResult{istrim}(cycleidx).wborder_rec(feaidx,movidx) 0 FitResult{istrim}(cycleidx).vborder_rec(feaidx,movidx) 0];
                                    outtab(cnt,:)=outline;
                                    xrec=[xrec;outline([3 5 7 9])];
                                end
                            end
                        end
                        % Mean data
                        cnt=cnt+1;
                        outline=[nc_range(cycleidx) 0 mean(xrec(:,1)) sqrt(var(xrec(:,1))) mean(xrec(:,2)) sqrt(var(xrec(:,2))) mean(xrec(:,3)) sqrt(var(xrec(:,3))) mean(xrec(:,4)) sqrt(var(xrec(:,4)))];
                        outtab(cnt,:)=outline;
                        % Confidence interval
                        cnt=cnt+1;
                        outline=[nc_range(cycleidx) -1 ...
                            FitResult{istrim}(cycleidx).xborder_CI(feaidx,1) FitResult{istrim}(cycleidx).xborder_CI(feaidx,2) ...
                            FitResult{istrim}(cycleidx).hborder_CI(feaidx,1) FitResult{istrim}(cycleidx).hborder_CI(feaidx,2) ...
                            FitResult{istrim}(cycleidx).wborder_CI(feaidx,1) FitResult{istrim}(cycleidx).wborder_CI(feaidx,2) ...
                            FitResult{istrim}(cycleidx).vborder_CI(feaidx,1) FitResult{istrim}(cycleidx).vborder_CI(feaidx,2) ...
                            ];
                        outtab(cnt,:)=outline;
                    end
                end
                outtab=outtab(1:cnt,:);
                outxls=mat2cell(outtab,ones(1,size(outtab,1)),ones(1,size(outtab,2)));

                % Put out new figure;
                if istrim==2
                    htmp=figure('Name',[feature_label{feaidx}  'for trimmed traces'],'Position',[100 100 800 500]);
                else
                    htmp=figure('Name',[feature_label{feaidx}  'for untrimmed traces'],'Position',[100 100 800 500]);
                end
                htmptable =   uitable('Parent',htmp,'Position',[0 0 800 500]);
                htmptable.Data = [header;outxls];
                set(htmptable,'ColumnEditable',true(1,10))
            end
            hall_plot_summary_Callback(istrim);
        end
    end

    function hall_plot_summary_Callback(istrim)
        for feaidx=fea_plot(:)'
            % Put out new figure for each feature;
            if istrim==2
                htmp=figure('Name',[feature_label{feaidx}  'for trimmed traces'],'Position',[100 100 800 500]);
            else
                htmp=figure('Name',[feature_label{feaidx}  'for untrimmed traces'],'Position',[100 100 800 500]);
            end                
            % Get number of nuclear cycle to plot
                nplot = 0;
                for cycleidx=1:numel(nc_range)
                    if numel(FitResult{istrim}(cycleidx).xborder_rec)
                        nplot = nplot+1;
                    end
                end
            % Get legend
                xlb = {};
                if fitoption_Hill(3)
                    xlb{1} = 'Border pos*';
                else
                    xlb{1} = 'Border pos';
                end
                if fitoption_Hill(2)
                    xlb{2} = 'Hill coeff*';
                else
                    xlb{2} = 'Hill coeff';
                end
                if fitoption_Hill(1)
                    xlb{3} = 'Plateau val*';
                else
                    xlb{3} = 'Plateau val';
                end
            % Plot statistics for each nuclear cycle
                mfea = [];
                ubfea=[];
                lbfea=[];
                cyclecnt=0;
                cyclelb={};     % Label for nc
                for cycleidx=1:numel(nc_range)
                    if numel(FitResult{istrim}(cycleidx).xborder_rec)
                        cyclecnt=cyclecnt+1;
                        % record the feature
                        xrec=[];
                        for movidx=1:Nmov
                            if size(FitResult{istrim}(cycleidx).hborder_rec,2)>=movidx
                                if FitResult{istrim}(cycleidx).hborder_rec(feaidx,movidx)
                                    xrec=[xrec;...
                                        [FitResult{istrim}(cycleidx).xborder_rec(feaidx,movidx) ...
                                        FitResult{istrim}(cycleidx).hborder_rec(feaidx,movidx) ...
                                        FitResult{istrim}(cycleidx).vborder_rec(feaidx,movidx)]];
                                end
                            end
                        end
                        % Mean data
                            mfea = [mfea; mean(xrec,1)];
                            fea_std = sqrt(var(xrec,[],1));
                        % Confidence interval
                        fea_bound = [FitResult{istrim}(cycleidx).xborder_CI(feaidx,1:2)' ...
                            FitResult{istrim}(cycleidx).hborder_CI(feaidx,1:2)' ...
                            FitResult{istrim}(cycleidx).vborder_CI(feaidx,1:2)' ];
                        % Show standard deviation or bounds
                            if fitoption_Hill(3)
                                x_lb(1) = fea_bound(1,1);
                                x_ub(1) = fea_bound(2,1);
                            else
                                x_lb(1) = fea_mean - fea_std(1);
                                x_ub(1) = fea_mean - fea_std(1);
                            end
                            if fitoption_Hill(2)
                                x_lb(2) = fea_bound(1,2);
                                x_ub(2) = fea_bound(2,2);
                            else
                                x_lb(2) = fea_mean - fea_std(2);
                                x_ub(2) = fea_mean - fea_std(2);
                            end
                            if fitoption_Hill(1)
                                x_lb(3) = fea_bound(1,3);
                                x_ub(3) = fea_bound(2,3);
                            else
                                x_lb(3) = fea_mean - fea_std(3);
                                x_ub(3) = fea_mean - fea_std(3);
                            end
                            ubfea = [ubfea;x_ub];
                            lbfea = [lbfea;x_lb];
                            cyclelb{cyclecnt} = ['nc' num2str(nc_range(cycleidx))];
                    end
                end
            % Plot everything 
                for j=1:3
                    subplot(3,1,j);
                    barplot_bound(mfea(:,j),lbfea(:,j),ubfea(:,j),cyclelb);
                    ylabel(xlb{j})
                end
        end
    end

    function hall_TimeMaginifier_Callback(~,~)
        h=figure;
        try
            [tlower,tupper,cycle_range,posborder]=Magnifier(h,heatmapI,binwidth);
            for i=1:numel(cycle_range)
                heatmapI(cycle_range(i)-8).posborder=posborder(1,i);
            end
            Re_Extract_feature(cycle_range,tlower,tupper);
            trimmed=true;
        catch
        end
    end
%% Auxiliary function
    function Re_Extract_feature(cycle_range,tlower,tupper)
        if nargin==0
            % Simple refresh
            cycle_range = nc_range;
            tlower = zeros(Nmov,numel(nc_range));
            tupper = zeros(Nmov,numel(nc_range))+10000;
        end
        cycle_range_=zeros(1,14);
        cycle_range_(cycle_range)=1;
        tlower_ = zeros(Nmov,14);tlower_(1:size(tlower,1),cycle_range)=tlower;
        tupper_ = zeros(Nmov,14)+10000;tupper_(1:size(tlower,1),cycle_range)=tupper;
        
        % Begin trimming:
        for i=1:numel(datamat)
            if cycle_range_(datamat(i).cycle)
                if (datamat(i).Feature(1)>=0)&&(tupper_(datamat(i).tscnt,datamat(i).cycle))
                    datamat(i).Feature=extract_feature(datamat(i).Intensity,datamat(i).time,...
                        0,datamat(i).dt,datamat(i).Imax,0,[tlower_(datamat(i).tscnt,datamat(i).cycle),tupper_(datamat(i).tscnt,datamat(i).cycle)]);
                end
            end
        end
    end
    
    function Update_Mean_Curves(~)
        pos_range=[AP(1):AP(2)];
        pos_range = [-45:45];
        % Get the mean curve for specific nuclear cycle
        if numel(datamat)   % If data is loaded
            % RESET the holders
            % Record the mean curves for merged embryos {feature x cycle} x [position]
            mf_rec={};  % Mean curve
            sf_rec={};  % Standard deviation
            nf_rec={};  % Number of nuclei
            ef_rec={};  % Number of embryo

            % Record the mean curves for individual embryos {feature x cycle x embryo} x [position]
            mf_indi={};  % Mean curve
            sf_indi={};  % Standard deviation
            nf_indi={};  % Number of nuclei
            cnt=0;
            for cycleno=nc_range
                cnt=cnt+1;
                ts_spec=find(arrayfun(@(x) getfield(DatasetList,{x},['nc' num2str(cycleno)]),1:Nmov));
                if numel(ts_spec)
                    cnt1=0;
                    for fea=1:numel(feature_label)
                        cnt1=cnt1+1;
                        % Put some alignment here
                        allx=[];allf=[];    % Merged feature
                        allx_={};allf_={};  % Individual feature
                        cnt2=0;
                        for tsidx=ts_spec
                            cnt2=cnt2+1;
                            xborder=0;
                            allx=[allx DatasetFeature(cnt).xaxis_all{fea,tsidx}-xborder];
                            allx_{cnt2}=DatasetFeature(cnt).xaxis_all{fea,tsidx}-xborder;
                            allf=[allf DatasetFeature(cnt).fearec_all{fea,tsidx}];
                            allf_{cnt2}=DatasetFeature(cnt).fearec_all{fea,tsidx};
                        end
                        % Calculate merged mean curve:
                            cnt3=0;
                            mf=[];  % mean
                            sf=[];  % standard deviation
                            nf=[];  % Number of samples
                            for pos=pos_range
                                cnt3=cnt3+1;
                                tmp=(allx-pos+binwidth/2).*(allx-pos-binwidth/2)<=0;
                                if sum(tmp)>Nsample_all_min
                                    mf(cnt3)=mean(allf(tmp));
                                    sf(cnt3)=sqrt(var(allf(tmp)));
                                    nf(cnt3)=sum(tmp);
                                else
                                    mf(cnt3)=NaN;
                                    sf(cnt3)=NaN;
                                    nf(cnt3)=0;
                                 end
                            end
                            mf_rec{fea,cycleno}=mf;
                            sf_rec{fea,cycleno}=sf;
                            nf_rec{fea,cycleno}=nf;
                            ef_rec{fea,cycleno}=pos_range*0;
                            
                        % Calculate individual mean curves                        
                            cnt2=0;
                            for tsidx=ts_spec
                                cnt2=cnt2+1;
                                cnt3=0;
                                mf=[];  % mean
                                sf=[];  % standard deviation
                                nf=[];  % Number of samples
                                for pos=pos_range
                                    cnt3=cnt3+1;
                                    tmp=(allx_{cnt2}-pos+binwidth/2).*(allx_{cnt2}-pos-binwidth/2)<=0;
                                    if sum(tmp)>Nsample_indi_min
                                        mf(cnt3)=mean(allf_{cnt2}(tmp));
                                        sf(cnt3)=sqrt(var(allf_{cnt2}(tmp)));
                                        nf(cnt3)=numel(allf_{cnt2}(tmp));
                                        ef_rec{fea,cycleno}=ef_rec{fea,cycleno}+1;
                                    else
                                        mf(cnt3)=NaN;
                                        sf(cnt3)=NaN;
                                        nf(cnt3)=0;
                                    end
                                end
                                mf_indi{fea,cycleno,tsidx}=mf;
                                sf_indi{fea,cycleno,tsidx}=sf;
                                nf_indi{fea,cycleno,tsidx}=nf;
                            end
                    end
                end
            end
            [~,tmp]=fileparts(DatasetName);
            if trimmed
                mkdir('tmp_trimmed');
                FitRes = FitResult{2};
                save(['tmp_trimmed/' tmp],'mf_rec','sf_rec','nf_rec','pos_range','heatmapI','mf_indi','sf_indi','nf_indi','FitRes');
            else
                mkdir('tmp');
                FitRes = FitResult{1};
                save(['tmp/' tmp],'mf_rec','sf_rec','nf_rec','pos_range','heatmapI','mf_indi','sf_indi','nf_indi','FitRes');
            end
        end
    end

    function hsetAP_Callback(~,~)
        AP(1)=str2num(get(hAPfrom,'String'));
        AP(2)=str2num(get(hAPto,'String'));
        if (AP(1)>-100)&&(AP(2)>-100)
            %Update_Mean_Curves;
        end
    end

    function hset_fea_ref(~,~)
        fea_ref=get(hfea_ref,'Value')-1;
    end

    function hset_binwidth(~,~)
        binwidth=str2num(get(hbinwidth,'String'));
    end

    function Update_DatasetList(idx)
        filename=fullfile(DatasetList(idx).Path,'correction.m');
        C1 = strsplit(DatasetList(idx).Path,'/');
        C2 = strsplit(DatasetList(idx).Path,'\');
        if numel(C2)>numel(C1)
            C1=C2;
        end
        if regexp(C1{end-1},'\d')
            DatasetList(idx).Name=C1{end-1};
        else
            if regexp(C1{end-2},'\d')
            DatasetList(idx).Name=C1{end-2};
            else
                DatasetList(idx).Name=C1{end-3};
            end
        end
        if exist(filename,'file')
            run(filename);
            DatasetList(idx).Lx=xlen;
            DatasetList(idx).Ly=ylen;
        else
            msgbox('File not found');
        end
    end

    function Update_List()
        set(hmessages,'String',fullfile(DatasetPath,DatasetName));
        set(hfea_ref,'Value',fea_ref+1);
        set(hbinwidth,'String',num2str(binwidth));
        datatmp=cell(Nmov,9);
        for i=1:Nmov
            datatmp{i,1}=['<html><body bgcolor=' rgb2hex(round(corder(i)*255)) '>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp</body></html>'];
            datatmp{i,2}=DatasetList(i).Name;
            datatmp{i,3}=DatasetList(i).shift_EL;
            datatmp{i,4}=DatasetList(i).nc9;
            datatmp{i,5}=DatasetList(i).nc10;
            datatmp{i,6}=DatasetList(i).nc11;
            datatmp{i,7}=DatasetList(i).nc12;
            datatmp{i,8}=DatasetList(i).nc13;
            datatmp{i,9}=DatasetList(i).nc14;
            datatmp{i,10}=DatasetList(i).BG;
            datatmp{i,11}=DatasetList(i).BG_man;
        end
        set(hAPfrom,'String',num2str(AP(1)));
        set(hAPto,'String',num2str(AP(2)));
        tmp=hfitoption.Data;
        tmp(:,2)=num2cell(fitoption_Hill');
        hfitoption.Data=tmp;
        htable.Data=datatmp;
    end

    function hfitoptionedit_Callback(~,~)
        fitoption_Hill = cell2mat(hfitoption.Data(:,2));
    end

    function hfeatableedit_Callback(~,~)
        fea_plot = find(cell2mat(hfeatable.Data(:,2)));
    end

    function hshowoptionedit_Callback(~,~)
        showoption_plot = cell2mat(hshowoption.Data(:,2));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELP HELP HELP HELP HELP HELP HELP HELP HELP HELP HELP HELP HELP %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function show_feature_meaning(~,~)
        Msg={'FEATURE DESCRIPTION:',...
          '1. ON/OFF', ...
          '2. Relative First Activation time ', ...
          '3. Relative Last Activation time ', ...
          '4. Relative Total active duration ', ...
          '5. Relative Total spot duration of existence ', ...
          '6. Relative Time to max intensity ', ...
          '7. Brightest spot intensity ', ...
          '8. Mean RNA produced ', ...
          '9. Total RNA produced ', ...
          '10. Total interphase duration ', ...
          '11. Absolute First Activation time ', ...
          '12. Absolute Last Activation time ', ...
          '13. Absolute Total active duration ', ...
          '14. Absolute Total spot duration of existence ', ...
          '15. Absolute Time to max intensity ', ...
            };
        msgbox(Msg,'Feature description','help')
    end
    %% Junk correction file:
    % Actually not needed, but needed anyway for the stupid static
    % workspace problem
    delimiter=' ';
    xlen=0;
    ylen=0;
    xlim_left=0;
    xlim_right=0;
    ylim_up=0;
    ylim_down=0;
    Icolumn=24;
    rm={};
    nclist={};
    range1=1:10001;filename1='Result_file1';
    range2=10002:30000;filename2='Result_file2';
    border13=10000;
    keepid13=[];
    leftcensored=0;
    rightcensored=0;
    startcutcycle=[0 0 0 0 0 0];
    endcutcycle=[10000 10000 10000 10000 10000 10000];
end
