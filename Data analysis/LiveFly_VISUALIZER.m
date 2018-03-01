function GETDATA
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

feature_label={};Nfea=15;
feature_unit={};
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
        'nc_ref',0,'Nf',0,'dt',0,'Lx',0,'Ly',0,'BG',1,'BG_man',1 ...
        ));
    DatasetList = emptystruct;
% Storage ensemble dataset
    datamat=struct('id',[],'Intensity',{},'Feature',[],'x',[],'y',[],'border',[],'mother',[],'daughter1',[],'daughter2',[],'tscnt',[],'cycle',[],'tlen',[],'start',[]);

% Storage for extracted features
    DatasetFeature=struct('idrec',[],'tsrec',[],'xrec',[],'yrec',[],'fearec',{},'xborder_rec',[],'hborder_rec',[],'wborder_rec',[],'vborder_rec',[],'xaxis_all',[],'yaxis_all',[],'fearec_all',{});
% Cell cycle range:
    nc_range=[9:14];                                % Valid nuclear cycle
    fea_ref=0;                                      % Reference feature for embryo alignment
    binwidth=3;                                     % Bin width
    tinterphase_min = [10 10 50 100 200 300];       % minimum interphase duration
    tinterphase_max = [1e4 1e4 1e4 1e4 1e4 1e4];    % maximum interphase duration
    fea_int = [8 9];                                % Which feature is intensity

% Default Color Order
    defcolor = ...
    [     0   0.4470    0.7410;...
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560;...
    0.4660    0.6740    0.1880;...
    0.3010    0.7450    0.9330;...
    0.6350    0.0780    0.1840;...
    1         0         0;...
    0         1         0;...
    0         0         1;...
         0   0.4470    0.7410;...
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560;...
    0.4660    0.6740    0.1880;...
    0.3010    0.7450    0.9330;...
    0.6350    0.0780    0.1840;...
    1         0         0;...
    0         1         0;...
    0         0         1;...
    ];
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
         htable.ColumnName = {'Color','Name','Ref cycle',...
             'nc9','nc10','nc11','nc12','nc13','nc14','BG','BG_man'};
         set(htable,'ColumnEditable',logical([0 0 1 1 1 1 1 1 1 0 1]));
         set(htable,'ColumnWidth',{80 120 60 40 40 40 40 40 40 60 60});
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
    hfea_ref = uicontrol('Style','popup','String',{'None','ON'},...
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
        1 ...                                       % 1: if sharing Hill coefficient
        ]);
    fitoption_Hill_text={'Similar maximum level', ...                      
        'Similar Hill coeff'};                                        
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
        'Callback',@hall_table_Callback);

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
    showoption_list={'Align',...
        'Show fitted Hill curve',...
        'Separate embryos by color',...
        'Show embryo legend',...
        'Show merged mean curve',...
        'Show indi mean curve',...
        'Normalize feature',...
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
                load(fullfile(PathName_,FileName_),'DatasetList','datamat','DatasetFeature','Nmov','AP','fea_ref','binwidth','tinterphase_min','tinterphase_max','fitoption_Hill');
                close(f);
                % Patch DatasetList if some new field is missing
                if ~isfield(DatasetList,'nc_ref')   % Reference feature
                    for i=1:numel(DatasetList)
                        DatasetList(i).nc_ref=0;
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
        catch
            msgbox('Error loading Dataset');
        end
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
            save(fullfile(DatasetPath,DatasetName),'DatasetList','datamat','DatasetFeature','Nmov','AP','fea_ref','binwidth','tinterphase_min','tinterphase_max','fitoption_Hill');
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
            save(fullfile(DatasetPath,DatasetName),'DatasetList','datamat','DatasetFeature','Nmov','AP','fea_ref','binwidth','tinterphase_min','tinterphase_max','fitoption_Hill');
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
                header(10:24)=feature_label;
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
    
    function htableedit_Callback(~,~)
        datatmp=htable.Data;
        for i=1:Nmov
            DatasetList(i).nc_ref=datatmp{i,3};
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
        f=waitbar(0,'Loading');
        if Nmov
            for i=1:Nmov
                waitbar((i-1)/Nmov,f);
                Update_DatasetList(Nmov);
            end
            [datamat,DatasetList]=local_extract(DatasetList);
        else
            msgbox('No movie loaded');
        end
        Update_List;
        close(f);
    end

    function hsetinterphase_Callback(~,~)
        Outtext=cell(1,numel(nc_range));
        if numel(datamat)
            for i=1:numel(nc_range)
                Outtext{i}=['nc' num2str(nc_range(i)) ' (min max) ' ];
                Deftext{i}=[num2str(tinterphase_min(i)) ' ' num2str(tinterphase_max(i))];
            end
            dlg=inputdlg(Outtext,'Set the range for interphase duration (s)',[1 50],Deftext);
            for i=1:numel(nc_range)
                tmp=str2num(dlg{i});
                tinterphase_min(i)=tmp(1);
                tinterphase_max(i)=tmp(2);
            end
            refine_tinterphase(tinterphase_min,tinterphase_max);
        else
            msgbox('Load movie first');
        end
    end
    
    function refine_tinterphase(tinterphase_min,tinterphase_max)
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
    end

    function hextractfeature_Callback(~,~)
        f=waitbar(0,'Loading');
        if numel(datamat)
            cnt=0;
            for cycleno=nc_range
                cnt=cnt+1;
                ts_spec=arrayfun(@(x) getfield(DatasetList,{x},['nc' num2str(cycleno)]),1:Nmov);
                [DatasetFeature(cnt).idrec,DatasetFeature(cnt).tsrec,DatasetFeature(cnt).xrec,DatasetFeature(cnt).yrec,DatasetFeature(cnt).fearec,...
                    DatasetFeature(cnt).xborder_rec,DatasetFeature(cnt).hborder_rec,DatasetFeature(cnt).wborder_rec,DatasetFeature(cnt).vborder_rec,...
                    DatasetFeature(cnt).xaxis_all,DatasetFeature(cnt).yaxis_all,DatasetFeature(cnt).fearec_all]=local_extract_feature(datamat,cycleno,find(ts_spec),0,feature_label,AP,fitoption_Hill);
            end
        end
        close(f);
    end
    
%% Show stuffs from individual movies
    function hindi_feature_Callback(~,~)
        % Summarize selected movies
        try
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
                            if strcmp(feature_unit{fea},'%')
                                allf=allf*100;
                                ratio=100;
                            else
                                ratio=1;
                            end
                            h=scatter(allx(:),ally(:),0*allf(:)+20,allf(:),'filled');
                            colorbar
                            set(h,'MarkerEdgeColor','k');
                            ylabel(['nc' num2str(nc_range(new_nc_range(cnt)))]);
                            xlim(AP);
                            if cnt==1
                                title(DatasetList(crrrow).Name,'interpreter','none');
                            end
                            if cnt<numel(new_nc_range)
                                set(gca,'XTick',[]);
                            else
                                xlabel('AP axis (%)');
                            end
                            subplot(numel(new_nc_range),2,2*cnt);
                            plot(allx,allf,'x');
                            xlim(AP);
                            if cnt==1
                                title([feature_label{fea} ' (' feature_unit{fea} ')']);
                            end
                            if cnt<numel(new_nc_range)
                                set(gca,'XTick',[]);
                            else
                                xlabel('AP axis (%)');
                            end
                            if showoption_plot(2)
                                hold on;
                                plot([AP(1):AP(2)],ratio*2*DatasetFeature(new_nc_range(cnt)).vborder_rec(fea,crrrow)*sigmf([AP(1):AP(2)],[ -DatasetFeature(new_nc_range(cnt)).hborder_rec(fea,crrrow)*0.04 DatasetFeature(new_nc_range(cnt)).xborder_rec(fea,crrrow)]),'--k');
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
        catch
            msgbox('Error. Reanalyze the data');
        end
    end
    
    function hindi_traces_Callback(~,~)
        Outtext={'ID from:','ID to','Pos from (%):','Pos to (%)','Only ON trace?:','Maximum Intensity'};
        Deftext={'0','10000','-50','50','0','3000'};
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
            if numel(idselect)<=80
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
            else
                msgbox(['Too more than 80 traces in cc' num2str(cycleno) '. Please narrow the criteria']);
            end
        end
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
                            if strcmp(feature_unit{feaidx},'%')
                                    ratio=100;
                                else
                                    ratio=1;
                                end
                            if numel(DatasetFeature(cnt0).fearec_all{feaidx,tsidx})
                                cnt=cnt+1;
                                % Check for border position if needed
                                if fea_ref & DatasetList(tsidx).nc_ref & showoption_plot(1)
                                    xborder=DatasetFeature(DatasetList(tsidx).nc_ref-8).xborder_rec(fea_ref,tsidx);
                                else
                                    xborder=0;
                                end
                                % Begin plotting
                                subplot(numel(ts_spec),numel(fea_plot),cnt)
                                plot(DatasetFeature(cnt0).xaxis_all{feaidx,tsidx}-xborder,ratio*DatasetFeature(cnt0).fearec_all{feaidx,tsidx},'.b');hold on;
                                tmpx=[AP(1):AP(2)];
                                % Show fitted Hill curve if needed
                                if showoption_plot(2)   % Show fitted Hill coefficient
                                    plot(tmpx,ratio*2*DatasetFeature(cnt0).vborder_rec(feaidx,tsidx)*sigmf(tmpx,[ -DatasetFeature(cnt0).hborder_rec(feaidx,tsidx)*0.04 DatasetFeature(cnt0).xborder_rec(feaidx,tsidx)-xborder]),'--k');
                                end
                                % Set axis limit
                                xlim(AP);
                                yrange{feaidx}=[];
                                if tsidx==ts_spec(1)
                                    title([feature_label{feaidx} ' (' feature_unit{feaidx} ')']);
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
            for cycleno=nc_range
                cnt=cnt+1;
                ts_spec=find(arrayfun(@(x) getfield(DatasetList,{x},['nc' num2str(cycleno)]),1:Nmov));
                               
                if numel(ts_spec)
                    cnt1=0;
                    figure('Name',['Dataset: ' DatasetName '. nc' num2str(cycleno) ]);
                    
                    for fea=fea_plot(:)'
                        if strcmp(feature_unit{fea},'%')
                            ratio=100;
                        else
                            ratio=1;
                        end
                        cnt1=cnt1+1;
                        subplot(numel(fea_plot),1,cnt1);
                        % Put some alignment here
                        allx=[];allf=[];
                        allx_={};allf_={};
                        cnt2=0;
                        for tsidx=ts_spec
                            cnt2=cnt2+1;
                            if fea_ref & DatasetList(tsidx).nc_ref & showoption_plot(1) % Align embryo
                                xborder=DatasetFeature(DatasetList(tsidx).nc_ref-8).xborder_rec(fea_ref,tsidx);
                            else
                                xborder=0;
                            end
                            allx=[allx DatasetFeature(cnt).xaxis_all{fea,tsidx}-xborder];
                            allx_{cnt2}=DatasetFeature(cnt).xaxis_all{fea,tsidx};
                            allf=[allf DatasetFeature(cnt).fearec_all{fea,tsidx}];
                            allf_{cnt2}=DatasetFeature(cnt).fearec_all{fea,tsidx};
                        end
                        cnt2=0;leg={};
                        if showoption_plot(3)   % Seperate embryo by color
                            for tsidx=ts_spec
                                cnt2=cnt2+1;
                                plot(allx_{cnt2},ratio*allf_{cnt2}, ...
                                        'Marker','.','color',defcolor(tsidx,:),'LineStyle','none'); hold on;
                                leg{cnt2}=['embryo ' num2str(tsidx)];
                            end
                        else
                            cnt2=cnt2+1;
                            plot(allx,ratio*allf,...
                                'Marker','.','LineStyle','none');
                            leg{cnt2}='data';
                        end
                        if showoption_plot(5)   % Show single mean curve
                            hold on;
                            pos_range=[AP(1):AP(2)];
                            cnt3=0;
                            mf=[];
                            sf=[];
                            for pos=pos_range
                                cnt3=cnt3+1;
                                tmp=(allx-pos+binwidth/2).*(allx-pos-binwidth/2)<=0;
                                mf(cnt3)=mean(allf(tmp));
                                sf(cnt3)=sqrt(var(allf(tmp)));
                            end
                            errorbar(pos_range,ratio*mf,ratio*sf,'k');
                            leg{cnt2+1}='Mean';
                        end
                        if showoption_plot(6)   % Show individual mean curves
                            hold on;
                            pos_range=[AP(1):AP(2)];
                            cnt2=0;
                            for tsidx=ts_spec
                                cnt2=cnt2+1;
                                cnt3=0;
                                mf=[];
                                sf=[];
                                for pos=pos_range
                                    cnt3=cnt3+1;
                                    tmp=(allx_{cnt2}-pos+binwidth/2).*(allx_{cnt2}-pos-binwidth/2)<=0;
                                    if numel(tmp)>3
                                        mf(cnt3)=mean(allf_{cnt2}(tmp));
                                        sf(cnt3)=sqrt(var(allf_{cnt2}(tmp)));
                                    else
                                        mf(cnt3)=NaN;
                                        sf(cnt3)=NaN;
                                    end
                                end
                                plot(pos_range,ratio*mf,'color',defcolor(tsidx,:));
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
            cnt=0;
            for cycleno=nc_range
                cnt=cnt+1;
                ts_spec=find(arrayfun(@(x) getfield(DatasetList,{x},['nc' num2str(cycleno)]),1:Nmov));
                if numel(ts_spec)
                    xborder=zeros(1,Nmov);
                    for tsidx=ts_spec
                        if fea_ref & DatasetList(tsidx).nc_ref & showoption_plot(1) % Align embryo
                            xborder(i)=DatasetFeature(DatasetList(tsidx).nc_ref-8).xborder_rec(fea_ref,tsidx);
                        end
                    end
                    time_align=0;time_normalize=0;
                    local_draw_kymo(datamat,ts_spec,xborder,cycleno,binwidth,AP,time_normalize,time_align);
                end
            end
        end
    end

    function hall_table_Callback(~,~)
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
        for feaidx=fea_plot(:)'
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
            outxls=mat2cell(outtab,ones(1,size(outtab,1)),ones(1,size(outtab,2)));
            
            % Put out new figure;
            htmp=figure('Name',feature_label{feaidx},'Position',[100 100 800 500]);
            htmptable =   uitable('Parent',htmp,'Position',[0 0 800 500]);
            htmptable.Data = [header;outxls];            
        end
    end
%% Auxiliary function
    function hsetAP_Callback(~,~)
        AP(1)=str2num(get(hAPfrom,'String'));
        AP(2)=str2num(get(hAPto,'String'));
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
            datatmp{i,1}=['<html><body bgcolor=' rgb2hex(round(defcolor(i,:)*255)) '>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp</body></html>'];
            datatmp{i,2}=DatasetList(i).Name;
            datatmp{i,3}=DatasetList(i).nc_ref;
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
    range1=1:10001;filename1='Result_file1';
    range2=10002:30000;filename2='Result_file2';
    border13=10000;
    keepid13=[];
    leftcensored=0;
    rightcensored=0;
    startcutcycle=[0 0 0 0 0 0];
    endcutcycle=[10000 10000 10000 10000 10000 10000];
end
