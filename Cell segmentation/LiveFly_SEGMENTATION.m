function LiveFly_SEGMENTATION
%% Author
% Mathieu Coppey
% Huy Tran
% Alec Douglas
%% ideas
% http://fr.mathworks.com/help/images/ref/fspecial.html?refresh=true
% http://www.mathworks.com/matlabcentral/answers/57077-how-do-you-perform-a-difference-of-gaussian-filter-on-an-image
% http://www.activovision.com/octavi/doku.php?id=integral_images
% http://fr.mathworks.com/help/images/examples/marker-controlled-watershed-segmentation.html
% http://fr.mathworks.com/help/images/ref/watershed.html
%% Parameters for positions and colors of figures and some handles
% The position of the game figure will be [figX1,figX2,figY1,figY2]
figX1 = 0;
figX2 = 1100;
figY1 = 0;
figY2 = 800;

%% sets initial variables

FileName = '';
PathName = '';
FilterIndex = '';
frame = 1;
frameinit = 1;
framefin = 150;
Nf = 1e5;           % Update Nf if only wish to mask till that frame
showBW = 0;
Img = [];
BW = [];BW=uint16(BW);
BWL = [];BWL=uint16(BWL);
Nnuc = [];
nucActive = {};
openflt = 1;
seg = 0;            % Propagation method
Cycle0=0;           % Cell cycle at frame 0
frametop=1;         % Top frame
framebottom=1000;   % Bottom frame
nuc = struct;       % Record of nuclei
zoommode=1;         % Mode of zooming in
brightness=1;       % Level of brightness in the image axis (not affecting segmentation)

tform_rec={};       % Storage for alignment
alignratio=1;       % Resize value before alignment
%% Create and hide the GUI figure as it is being constructed.
segfigure = figure('Visible','on','Tag','segfigure','Position',[figX1,figY1,figX2,figY2]);
set ( gcf, 'Color', [0 0 0] )

% To capture mouse and keyboard
set(segfigure,'WindowKeyPressFcn',@keypress);
set(segfigure,'WindowScrollWheelFcn',@scroll);
set(segfigure,'WindowButtonUpFcn',@click);
%% Create menu
set(gcf,'Menubar','None');
% Image or mask file loading, saving:
mFile = uimenu('Label','File_');

    hloadbutton = uimenu(mFile,'Label','Load image',...
    'Callback',@hloadbutton_Callback); % load file

    % Buttons for SAVE/LOAD MASK and TRACK files
    hsavebw = uimenu(mFile,'Label','Save Mask',...
    'Callback',@hsavebw_Callback);

    hloadbw = uimenu(mFile,'Label','Load Mask',...
    'Callback',@hloadbw_Callback);

    hquitbutton = uimenu(mFile,'Label','Quit',...
    'Callback',@hquitbutton_Callback); % Quit seg

% Edit mask:
mView = uimenu('Label','View_');
    
    hzoom_ = uimenu(mView,'Label','Zoom in/out',...
    'Callback',@hzoom_Callback);

    hshow = uimenu(mView,'Label','Cycle View (X)',...
    'Callback',@hshow_Callback);

    hcellcount = uimenu(mView,'Label','View cell count',...
        'Callback',@hcellcount_Callback);
    
    hbrightness = uimenu(mView,'Label','Adjust intensity:');
        
        hbrightup = uimenu(hbrightness,'Label','Increase brightness (1.5x)',...
            'Callback',@hbrightup_Callback);
        hbrightdown = uimenu(hbrightness,'Label','Decrease brightness (1.5x)',...
            'Callback',@hbrightdown_Callback);
        hbrightreset = uimenu(hbrightness,'Label','Reset brightness',...
            'Callback',@hbrightreset_Callback);
        
% Edit mask: segment, adjust mask size, add, remove
mEdit = uimenu('Label','Edit_');

    hremove = uimenu(mEdit,'Label','Remove Mask (R)',...
    'Callback',@hremove_Callback);

    hadd_ = uimenu(mEdit,'Label','Add new mask:');
        
        haddfree = uimenu(hadd_,'Label','Free hand (A)',...
            'Callback',@haddfree_Callback);
        
        haddelip = uimenu(hadd_,'Label','Elipse (Shift+A)',...
            'Callback',@haddelip_Callback);
    
        hsegment_ = uimenu(mEdit,'Label','Segment (S)',...
            'Callback',@hsegment_Callback);
        
        hpropplus = uimenu(mEdit,'Label','Propagate Next (Space)',...
            'Callback',@hpropplus_Callback);
        
        hpropminus = uimenu(mEdit,'Label','Propagate Prev (Shift+Space)',...
            'Callback',@hpropminus_Callback);
        
        hrealign = uimenu(mEdit,'Label','Realign',...
            'Callback',@hrealign_Callback);
% Cell tracking
mTrack = uimenu('Label','Track_');
    htrack = uimenu(mTrack,'Label','Track',...
        'Callback',@htrack_Callback);

    hchecktrack = uimenu(mTrack,'Label', 'Check Track',...
        'Callback',@hchecktrack_Callback);

    hvisualizetract = uimenu(mTrack,'Label','Plot lineage tree',...
        'Callback',@hvisualizetrack_Callback);

    hsetcycle = uimenu(mTrack,'Label', 'Set cycles',...
        'Callback',@hsetcycle_Callback);
    
% Help
mHelp = uimenu('Label','Help_');
%% Create buttons and others
    
hmessages = uicontrol('Style','text','String','Load a movie or press F1 for Help',...
    'Position',[50,730,400,50],'FontSize',14, ...
    'ForegroundColor','white','BackgroundColor','black'); % Display instructions

hpath = uicontrol('Style','text','String','Path',...
    'Position',[610,730,400,50],'FontSize',10);

hframe = uicontrol('Style','slider',...
    'Min',1,'Max',Nf, 'value',frame,...
    'Position',[100,10,900,20],...
    'SliderStep',[1/150 0.2],'Callback',@hframe_Callback);

hzoom = uicontrol('Style','pushbutton','String', 'Zoom',...
    'Position',[0,10,90,20],...
    'Callback',@hzoom_Callback);

hviewmode = uicontrol('Style','text','String','View mode: 0',...
    'Position',[1010,10,90,20],'FontSize',8);
% Segmentation parameters
htextchoice = uicontrol('Style','text','String','Seg. params',...
    'Position',[20,750,80,30],'FontSize',10, ...
    'ForegroundColor','white','BackgroundColor','black');

hfill = uicontrol('Style','checkbox','String', 'Fill',...
    'Position',[50,730,50,20]);

hpropchoice = uicontrol('Style', 'popup',...
           'String', 'Don''t use prev/next|Use prev|Use next',...
           'Position', [120 720 130 30],...
           'Callback', @hpropchoice_Callback);
htextfilter = uicontrol('Style','text','String','Nuclei size',...
    'Position',[40,690,80,30],'FontSize',10, ...
    'ForegroundColor','white','BackgroundColor','black');       
hfilter = uicontrol('Style','slider',...
    'Min',1,'Max',40, 'value',openflt,...
    'Position',[120,700,130,20],...
    'SliderStep',[1/40 0.2],'Callback',@hfilter_Callback);
htextthresh = uicontrol('Style','text','String','Threshold',...
    'Position',[40,660,80,30],'FontSize',10, ...
    'ForegroundColor','white','BackgroundColor','black');
hthresh = uicontrol('Style','slider',...
    'Min',0,'Max',255, 'value',openflt,...
    'Position',[120,670,130,20],...
    'SliderStep',[1/256 0.1],'Callback',@hfilter_Callback);

% Segmentation buttons
hsegment = uicontrol('Style','pushbutton','String', 'Segment (S)',...
    'Position',[260,700,60,20],...
    'Callback',@hsegment_Callback);

hsegchoice = uicontrol('Style', 'popup',...
           'String', 'Thresh1|Thresh2|Circle',...
           'Position', [260 670 60 20]);
       
hbwplus = uicontrol('Style','pushbutton','String', 'BW+',...
    'Position',[260,730,30,20],...
    'Callback',@hbwplus_Callback);

hbwminus = uicontrol('Style','pushbutton','String', 'BW-',...
    'Position',[290,730,30,20],...
    'Callback',@hbwminus_Callback);

hcheckecc = uicontrol('Style','checkbox','String', 'Show ecc',...
    'Position',[340,730,80,20]);

% Auto adjust size:
hautosize = uicontrol('Style','pushbutton','String', 'Auto-size:',...
    'Position',[450,700,40,20],...
    'Enable', 'off');

% Auto adjust size value:
hautosizef = uicontrol('Style','edit','String', '0',...
    'Position',[490,700,30,20]);

% Top frame - remove nuclei at the top border
htopframe = uicontrol('Style','pushbutton','String', 'Top. fr',...
    'Position',[340,700,40,20],...
    'Enable', 'off');

% Top frame - remove nuclei at the bottom border
hbottomframe = uicontrol('Style','pushbutton','String', 'Bot. fr',...
    'Position',[340,670,40,20],...
    'Enable', 'off');

% Value for top frame
htopframef = uicontrol('Style','edit','String', num2str(framebottom),...
    'Position',[390,670,30,20],...
    'Callback',@htopframef_Callback);

% Value for bottom frame
hbottomframef = uicontrol('Style','edit','String', num2str(frametop),...
    'Position',[390,700,30,20],...
    'Callback',@hbottomframef_Callback);


%% Create axes and display image
ImgOp = zeros(100,100);
sI = size(ImgOp);
panelI_X1 = 40;
panelI_X2 = panelI_X1+sI(2);
panelI_Y1 = 40;
panelI_Y2 = panelI_Y1+sI(1);

ha = axes('Parent',segfigure);
imshow(ImgOp);

hbw = axes('Parent',segfigure);
BWop = zeros(100,100);
imshow(BWop);
% A mouse cursor here
curradi=10;
hcursor = viscircles(ha,[0 0],curradi);

% Change units to normalized so components resize automatically.
set([segfigure,ha,hmessages...
    hbw,hpath,...
    htextfilter,htextthresh,htextchoice,hframe,hfilter,hthresh,...
    hsegment,hsegchoice,hbwplus,hbwminus,hfill,hpropchoice...
    hcheckecc,hzoom,hviewmode,...
    htopframe,hbottomframe,htopframef,hbottomframef,...
    hautosize,hautosizef,...
    ],...
    'Units','normalized');
%% Final settings
% Assign the GUI a name to appear in the window title.
set(segfigure,'Name','Segment NUP V2 (Press F1 for help)') 
% Move the GUI to the center of the screen.
movegui(segfigure,'center')
% Make the GUI visible.
set(segfigure,'Visible','on')
%% Callbacks
    function click(~,~)
        % Not used...
    end
    function hquitbutton_Callback(~,~)
        close(segfigure);
    end

    function hloadbutton_Callback(~,~)
        [FileName_,PathName_,FilterIndex_] = uigetfile('*.tif','Select the tif file');
        if FileName_
            % Prompt for file information
                FileName=FileName_;
                PathName=PathName_;
                FilterIndex=FilterIndex_;            
                info = imfinfo([PathName,FileName]);
                inpNf=inputdlg({'Enter max frame number (Type 0 if you want analyze the whole movie)',...
                    'Cell cycle at frame 1'});
                if str2num(inpNf{1})
                    Nf= str2num(inpNf{1});
                else
                    Nf = numel(info);
                end
                if str2num(inpNf{2})
                    Cycle0= str2num(inpNf{2});
                else
                    Cycle0 = 0;
                end
                Nf = min([Nf,numel(info)]);
                Wd = info(1).Width;
                Hg = info(1).Height;
                sI = [Hg,Wd];
                Img = [];
                BW = zeros(Hg,Wd,Nf,'uint16');
                BWL = zeros(Hg,Wd,Nf,'uint16');
                Nnuc = zeros(Nf,1);
                nucActive=cell(1,Nf);
            % Loading the image
                h = waitbar(0,'Loading image');
                for i=1:Nf
                    waitbar(i/Nf,h)
                    Img = cat(3,Img,imread([PathName,FileName],i));
                end
                close(h);
            % Setup the axis
                distavai=get(hsegchoice,'Position');
                sI = size(Img(:,:,1));
                sI_=(distavai(2)-0.07)*sI/sI(1)/2;
                figaspect=get(segfigure,'PaperPosition');
                figaspect=figaspect(4)/figaspect(3);
                panelI_X1 = 0.05;
                panelI_X2 = sI_(2)*figaspect;
                panelI_Y1 = 0.05;
                panelI_Y2 = sI_(1);
                zoommode=0;
                set(ha,'Units','Normalized');
                set(ha,'Position',[panelI_X1,panelI_Y1,panelI_X2,panelI_Y2]);
                axis normal
                set(hbw,'Units','Normalized');
                set(hbw,'Position',[panelI_X1,panelI_Y1+sI_(1),panelI_X2,panelI_Y2]);
                BW(:,:,1) = 0;
                axis normal
                showI(Img(:,:,1),BW(:,:,1));
                set(hframe,'max',Nf,'SliderStep',[1/(Nf-1) 0.2],'Value',1);
                frametop=1;
                framebottom=sI(1);
                set(htopframef,'String',num2str(1));
                set(hbottomframef,'String',num2str(sI(1)));
                set(hpath,'String',[PathName,FileName,'. Nframe = ' num2str(Nf)]);
                set(hmessages,'String',['Loaded. Frame ',num2str(frame)]);
            % Load alignment information if any
                k = strfind(FileName,'.');
                if exist([PathName,FileName(1:k(end)-1) '_align.mat'],'file')
                    load([PathName,FileName(1:k(end)-1) '_align.mat'],'tform_rec','alignratio');
                else
                    tform_rec=cell(1,Nf);   % Transformation info
                    alignratio=1;           % Resize
                end
        end
    end

    function hframe_Callback(~,~)
        frame = get(hframe,'Value');
        frame = round(frame);
        %set(hframe,'Value',frame);
        showI;        
    end

    function scroll(~,evnt)
        deltaframe = evnt.VerticalScrollCount;
        frame = frame + deltaframe;
        if frame<1
            frame=1;
        end
        if frame>Nf
            frame = Nf;
        end
        set(hframe,'Value',frame);        
        showI;
    end

    function hpropchoice_Callback(~,~)
        choice = get(hpropchoice,'value');
        switch choice 
            case 1
                seg = 0;% Segment as usual. New cells can emerge. This should be used during e.g. first frame, mitosis
            case 2
                seg = 1;% Use previous frame's mask for the segmentation. No new cell will be added
            case 3
                seg = 2;% Use next frame's mask for the segmentation. No new cell will be added
        end
    end

    % Add filter to the current background images
    function hfilter_Callback(~,~)
        img=Img(:,:,frame);
        % Filter the image
        openflt = round(get(hfilter,'Value'));
        set(hfilter,'Value',openflt);
        ratiostrel=1.5;
        f1 = fspecial('Gaussian', openflt, openflt/3);
        % Extract the background from image:
        background = imopen(img,strel('disk',round(openflt*ratiostrel)));
        img=imfilter(img-background,f1,'replicate');
        % Perform thresholding
        openthresh = round(get(hthresh,'Value'))/255;
        if openthresh==0
            openthresh=graythresh(img);
        end
        % Increasing contrast only:
        %BW(:,:,frame) = im2bw(img,openthresh);
        BWtmp=(im2bw(img,openthresh));
        BW(:,:,frame) = double(BWtmp);
        if get(hfill,'Value')
            BW(:,:,frame) = imfill(BW(:,:,frame),'holes');
        end
        showop(BW(:,:,frame));
        set(hmessages,'String',['Frame ',num2str(frame)]);
        segmode=0;
    end
    
    % Show button - add mask over the image in Panel 2.
    function hshow_Callback(~,~)
        showBW=mod(showBW+1,4);
        if (showBW==3)&&(frame==1)
            showBW=0;
        end
        set(hviewmode,'String',['View mode: ' num2str(showBW)]);
        showI;
    end

    % Remove a mask - Modification on Panel 1.
    function hremove_Callback(~,~)
        if zoommode==1
            [x,y] = getpts(ha);
        else
            [x,y] = getpts(hbw);
        end
        bwi=BW(:,:,frame);
        bwi(bwi==BW(round(y),round(x),frame))=0;
        BW(:,:,frame) = bwi;
        showI;
    end
    
    % Add a mask - either freehand or elipsoid mask
    function haddfree_Callback(~,~)
        hadd_Callback(1)
    end
    function haddelip_Callback(~,~)
        hadd_Callback(0)
    end
    % Add a mask
    function hadd_Callback(isdraw)
        if zoommode==2
            r = getrect(hbw);
        else
            r = getrect(ha);
        end
        if isdraw
            tempfig = figure('Name','please draw the nucleus');
            locI = imcrop(Img(:,:,frame),r);
            locBW= double(imcrop(BW(:,:,frame),r));
            locsize = size(locI);
            ImBW = zeros(locsize(1),locsize(2),3);
            ImBW(:,:,1)=double(locI);
            ImBW(:,:,2)=double(locI);
            ImBW(:,:,3)=double(locI);
            ImBW(:,:,2)=ImBW(:,:,2).*imcomplement(locBW);
            ImBW(:,:,3)=ImBW(:,:,3).*imcomplement(locBW);
            imshow(uint8(ImBW));
            % Resize the frame 2x for better draw
            framesize=get(gca,'Position');zoomlv=2;
            set(gca,'Position',[0 0 framesize(3:4)*zoomlv]);
            % Draw a nucleus
            addingreg = imfreehand;
            BWtemp = createMask(addingreg);
            close(tempfig)
        else
            % Create an ellipse
            BWtemp = zeros(round(r(4)),round(r(3)));
            for i=1:r(4)
                for j=1:(r(3))
                    if ((i-r(4)/2)/r(4))^2+((j-r(3)/2)/r(3))^2<0.25
                        BWtemp(i,j)=1;
                    end
                end
            end
        end
        BWtemp=BWtemp>0;
        r = round(r);
        rs = size(BWtemp>0);
        r(3) = rs(2)-1;
        r(4) = rs(1)-1;
        newidx=max(max(BW(:,:,frame)))+1;
        bwi=BW(r(2):r(2)+r(4),r(1):r(1)+r(3),frame);
        bwi(BWtemp)=newidx;
        BW(r(2):r(2)+r(4),r(1):r(1)+r(3),frame) = bwi;
        showI;
    end

    % Realign the current frame
    function hrealign_Callback(~,~)
        ncnt=30;
        if frame>1
            cont=1;
            while cont
                Icrr=imadjust(imresize(Img(:,:,frame),alignratio));
                Ipre=imadjust(imresize(Img(:,:,frame-1),alignratio));
                
                % Perform some alignment
               
                [tform,~]=imregdemons(Ipre,Icrr,[ncnt ncnt ncnt ncnt ncnt],'AccumulatedFieldSmoothing',1,'DisplayWaitBar',true,'PyramidLevels',5);
                % Convert to polynomial surface projection
                [x,y]=meshgrid(1:size(Icrr,2),1:size(Icrr,1));
                x_=tform(:,:,1);y_=tform(:,:,2);
                tform_poly=struct;
                tform_poly.x = fit([x(:),y(:)],x_(:),'poly55');
                tform_poly.y = fit([x(:),y(:)],y_(:),'poly55');

                % Update the alignment file

                h=figure;
                subplot(211);
                imshowpair(uint8(Ipre),uint8(Icrr));
                title('Before alignment');
                subplot(212);
                [Itmp,rd]=imwarp_poly(Ipre,tform_poly);
                Inew=resetframe(Itmp,rd,size(Ipre));
                imshowpair(uint8(Inew),uint8(Icrr));
                set(h,'Position',[100 100 800 800]);
                title('After alignment');
                % Confirm the new alignment:
                pt=questdlg('Is this a good alignment?','Option',...
                    'Yes','No - Redo', 'Manual','Cancel');
                cont=0;
                switch pt
                    case 'Yes'
                        tform_rec{frame}=tform_poly;
                        cont=0;
                    case 'No - Redo'
                        cont=1;
                    case 'Manual'
                        tform_rec{frame}=mannual_align(Ipre,Icrr,h);
                        cont=0;
                    case ''
                        tform_rec{frame}=[];
                        cont=0;    
                end
                ncnt=ncnt+50;
                close(h);
            end
            showI;
        else
            msgbox('No previous frame');
        end
    end
    
    % Manual alignment:
    function [tform]=mannual_align(Ipre,Icrr,h)
        cont=1;
        prept=[];
        crrpt=[];
            while cont
                if numel(prept)>10
                    [prept,crrpt] =cpselect(Ipre,Icrr,prept,crrpt,'wait',true);
                else
                    [prept,crrpt] =cpselect(Ipre,Icrr,'wait',true);
                end
                [prept] = cpcorr(prept,crrpt,Ipre,Icrr);
                tform = fitgeotrans(prept,crrpt,'projective');
                figure(h);
                subplot(212);
                [Itmp,rd]=imwarp_poly(Ipre,tform);
                Inew=resetframe(Itmp,rd,size(Ipre));
                imshowpair(uint8(Inew),uint8(Icrr));
                pt=questdlg('Is this a good alignment?','Option',...
                    'Yes','No - Continue','Cancel','Cancel');
                cont=0;
                switch pt
                    case 'Yes'
                        cont=0;
                    case 'No - Continue'
                        cont=1;
                    case 'Cancel'
                        tform=[];
                        cont=0;
                    case ''
                        tform=[];
                        cont=0;
                end
            end
        
    end
    % Use the current mask to segment
    function hsegment_Callback(~,~)
        switch get(hsegchoice,'value')
            case 1  % Thresh1 - Threshold with small nuclei removal
                % Threshold method
                segment_Thresh1();
            case 2  % Thresh2 - Like Thresh1 but without small nuclei removal
                % Threshold method
                segment_Thresh1();
            case 3  % Circle - Detect circle
                % Find circle
                segment_findcircle();
                segment_Thresh1();
        end
    end

    function segment_Thresh1()
        BW(:,:,frame) = segment(BW(:,:,frame));
        bwi=imclearborder(BW(frametop:framebottom,:,frame));
        bwi=[zeros(frametop-1,size(bwi,2));bwi;zeros(sI(1)-framebottom,size(bwi,2))];
        ref_frame=0;
        switch seg
            case 0% Do nothing
            case 1% Use previous frame for reference - seek alignment if possible
                if frame>1
                    ref_frame=BW(:,:,frame-1); % Load previous mask
                    if numel(tform_rec{frame})
                        ref_frame=applyimwarp(ref_frame,tform_rec{frame},alignratio);
                    end
                    ref_frame=ref_frame>0;
                end
            case 2% Use next frame for reference
                if frame<Nf
                    ref_frame=BW(:,:,frame+1)>0;
                end
        end
        % Check overlap between ref_frame and current frame
        if numel(ref_frame)>1
            tmplb=unique(double(ref_frame).*double(bwi));
            tmplb_=unique(bwi);
            tmplb=setdiff(tmplb_,tmplb);
            for i=tmplb(:)'
                bwi(bwi==i)=0;
            end
        end
        BW(:,:,frame) = (bwi);
        % Increase the radius
%         tmp=radiusIncrease;
%         radiusIncrease=0;
%         for i=1:abs(tmp)
%             if tmp>0
%                 hbwplus_Callback();
%             else
%                 hbwminus_Callback();
%             end
%         end
        showI;
    end

function segment_findcircle(~,~)
        sharpenImg = imadjust(Img(:,:,frame));
        %sizeTemp=size(Img);

        %[x,y] = ginput(2);

        %    radius = int8(sqrt((x(1)-x(2)).^2+(y(1)-y(2)).^2))
        radius = uint8(round(get(hfilter,'Value')));

        [circleSpot,radiiSpot] = imfindcircles(sharpenImg,...
            [radius-radius/6,radius+radius/3],'Sensitivity',.92,...
            'Method','twostage','EdgeThreshold',.2);
        
        %showI_(Img(:,:,frame),sharpenImg);
        %viscircles(circleSpot,radiiSpot);
        
        % Convert circleSpot and radiiSpot to colored mask
        BWtmp=zeros(size(BW(:,:,frame)));
        [X,Y] = meshgrid(1:size(BWtmp,2),1:size(BWtmp,1));
        for i=1:numel(radiiSpot)
            Itmp = (X-circleSpot(i,1)).^2+(Y-circleSpot(i,2)).^2 < radiiSpot(i)^2;
            BWtmp=BWtmp+Itmp;
        end
        BW(:,:,frame)=double(BWtmp>0);
        showop(BW(:,:,frame));
        segmode=1;
    end

    % Expanding the size of the mask
    function hbwplus_Callback(~,~)
        %radiusIncrease=radiusIncrease+1;
        nhood=[1 1 1;1 1 1;1 1 1];
        se = strel('arbitrary',nhood);
        BWtmp=imdilate(BW(:,:,frame),se);
        % Check for overlap: do not invade if there exists cell before
        tmp1=BW(:,:,frame);
        tmp2=find(tmp1);
        BWtmp(tmp2)=tmp1(tmp2);
        BW(:,:,frame)=BWtmp;
        showI;
    end

    % Decrease the size of the mask
    function hbwminus_Callback(~,~)
        %radiusIncrease=radiusIncrease-1;
        nhood=[1 1 1;1 1 1;1 1 1];
        bwi=conv2(1-(BW(:,:,frame)>0),nhood,'same')==0;
        BW(:,:,frame)=BW(:,:,frame).*uint16(bwi);
        showI;
    end

    function htopframef_Callback(~,~)
        frametop = str2num(get(htopframef,'String'));
        showI;
    end

    function hbottomframef_Callback(~,~)
        framebottom = str2num(get(hbottomframef,'String'));
        showI;
    end

    function hsavebw_Callback(~,~)
        k = strfind(FileName,'.');
        maskfile=[PathName,FileName(1:k(end)-1),'_Mask_' num2str(Nf) '.mat'];
        alignfile=[PathName,FileName(1:k(end)-1),'_align.mat'];
        if ~exist(maskfile,'file')
            save(maskfile,'BW','-v7.3');
        else
            if strcmp(questdlg('Overwrite mask file?','File exists','Yes','No','Yes'),'Yes')
                save(maskfile,'BW','-v7.3');
                if strcmp(questdlg('Update the alignment file?','File exists','Yes','No','Yes'),'Yes')
                    save(alignfile,'tform_rec','alignratio','-v7.3');
                end
            end
        end
    end

    function hloadbw_Callback(~,~)
        k = strfind(FileName,'.');
        try
            temp = load([PathName,FileName(1:k(end)-1),'_Mask_' num2str(Nf) '.mat']);
            BW = temp.BW;
            showI
        catch
            msgbox('No mask file found');
        end
    end
    % Show cell count over time
    function hcellcount_Callback(~,~)
        cnt_all=zeros(1,Nf);
        cnt_wmother=zeros(1,Nf);
        h = waitbar(0,'Calculating cell count');
        for i=1:Nf
            waitbar(i/Nf,h);
            tmp1=unique(BW(:,:,i));
            cnt_all(i)=sum(tmp1>0);
            if (cnt_all(i)==0)||(i==1)
                cnt_wmother(i)=0;
            else
                % Perform the alignment of prev frame
                ref_frame=BW(:,:,i-1);
                if numel(tform_rec{i})
                    ref_frame=applyimwarp(ref_frame,tform_rec{i},alignratio);
                end
                [~,wmother]=cellcount(BW(:,:,i),ref_frame);
                cnt_wmother(i)=numel(wmother);
            end
        end
        close(h);
        figure;
        plot(cnt_all,'r');hold on;
        plot(cnt_wmother,'b');
        xlabel('Frame');
        ylabel('Nuclei count');
        h=legend('Nuclei count','Nuclei with mother count');
        set(h,'Location','NorthWest');
    end
    
    % Track nuclei identity
    function htrack_Callback(~,~)
        if ~strcmp('Yes',questdlg('You should view cell count first. Continue?','Tracking...','Yes','No','No'))
            return;
        end
        BW=uint16(BW);
        initcountnuc=Nf*6;  % Estimate nuclei number
        nuc=struct;
        % Initiating the result recorder
        nuc.frames = zeros(initcountnuc,Nf,'uint8');        % frame: existing or not
        nuc.x = zeros(initcountnuc,Nf);                     % Position x
        nuc.y = zeros(initcountnuc,Nf);                     % Position y
        nuc.ind = zeros(initcountnuc,Nf,'uint16');          % Id of nuclei
        nuc.parent=zeros(initcountnuc,Nf,'uint16');         % Parent of nuclei
        nuc.daughter1=zeros(initcountnuc,Nf,'uint16');      % Daughter 1 of nuclei
        nuc.daughter2=zeros(initcountnuc,Nf,'uint16');      % Daughter 2 of nuclei
        nuc.area = zeros(initcountnuc,Nf);                  % Area of nuclei        
        nuc.radius = zeros(initcountnuc,Nf);                % Radius (approximated by nuclei size)
        nuc.ecc = zeros(initcountnuc,Nf);                   % Eccentricity        
        nuc.ang = zeros(initcountnuc,Nf);                   % Angle
        nuc.cycle = zeros(initcountnuc,Nf);                 % Cell cycle

        % Create a panel to vísualize the process
        %tempfig = figure('Name','tracking...');
        locframe = 1;       % Identifier for currently processed frame
        
        CC = colorconncomp(BW(:,:,locframe));
        stats = regionprops(CC,'Centroid','Eccentricity','Orientation');
        Nnuc(locframe) = CC.NumObjects;
        countnuc = 1;
       

        % Creating nuclie information for the first panel
        bwtemp = zeros(sI);
        for j=1:Nnuc(locframe)
            nuc.frames(countnuc,locframe) = 1;
            nuc.x(countnuc,locframe) = stats(j).Centroid(1);
            nuc.y(countnuc,locframe) = stats(j).Centroid(2);
            nuc.ind(countnuc,locframe) = countnuc;
            nuc.parent(countnuc,locframe) = 0;
            nuc.daughter1(countnuc,locframe) = 0;
            nuc.daughter2(countnuc,locframe) = 0;
            nuc.area(countnuc,locframe) = numel(CC.PixelIdxList{j});            
            nuc.radius(countnuc,locframe) = sqrt(nuc.area(countnuc,locframe)/pi);
            nuc.ecc(countnuc,locframe) = stats(j).Eccentricity;
            nuc.ecc(countnuc,locframe) = stats(j).Orientation;
            bwtemp(CC.PixelIdxList{j}) =  nuc.ind(countnuc,locframe);
            countnuc = countnuc+1;
        end        
        BWL(:,:,locframe) =  bwtemp;    
        maxID=max(unique(bwtemp));
        % Update the list of active nuclei
        nucActive{locframe} = find(nuc.frames(:,locframe)==1);

        h=waitbar(1/Nf,['Frame: ' num2str(1) '/' num2str(Nf)]);
        % Processing following frames
        for i=2:Nf
            waitbar(i/Nf,h,['Frame: ' num2str(i) '/' num2str(Nf)]);
            locframe = i;
            %tempfig
            crr_frame=BW(:,:,locframe);
            CC = colorconncomp(crr_frame);
            stats = regionprops(CC,'Centroid','Eccentricity','Orientation');
            Nnuc(locframe) = CC.NumObjects;
            
            % Create list of prev cell - curr cell
            hits = zeros(Nnuc(locframe),1);
            parenthits = zeros(Nnuc(locframe),1);
            ref_frame=BWL(:,:,locframe-1);
            if numel(tform_rec{locframe})
                ref_frame=applyimwarp(ref_frame,tform_rec{locframe},alignratio);
            end
            Pre_edge = unique(ref_frame);   % Cell label in prev frame
            Crr_edge = double(unique(crr_frame));   % Cell label in current frame
            [Overlap_count]=histcounts2(crr_frame(:),ref_frame(:),[-1;Crr_edge]+0.5,[-1;Pre_edge]+0.5);
            
            for j=1:Nnuc(locframe)
                % Find the current index in Pre_edge
                cellid=find(Crr_edge==crr_frame(CC.PixelIdxList{j}(1)),1,'first');
                ident=find(Overlap_count(cellid,:));
                cntident=Overlap_count(cellid,ident);
                ident=Pre_edge(ident);
                cntident=cntident(ident>0);
                ident=ident(ident>0);
                % Find the current cell mask and check if overlap with previous frame
                    %bwtmp = (crr_frame==crr_frame(CC.PixelIdxList{j}(1))).*ref_frame;
                    % Check if overlap with two cells - choose the bigger one
                    %ident = sort(unique(bwtmp));
                    %cntident=hist(bwtmp(:),ident);
                    %cntident=cntident(ident>0);
                    %ident=ident(ident>0);
                if numel(ident)
                    if numel(ident)>1
                        [~,tmp]=max(cntident);
                    else
                        tmp=1;
                    end
                    % Assign the cell ancestor
                    hits(j)=ident(tmp(1));
                end
            end
            % Recreate new label mask
            bwtemp = zeros(sI);
            for j=1:Nnuc(locframe)
                % If there are no new parents
                if hits(j) == 0
                    maxID=maxID+1;
                    newID = maxID;
                    parenthits(j)= 0;
                end
                % If cell has a parent
                if (hits(j) >0)
                    % Check cell division
                    if sum(hits==hits(j))>1
                        maxID=maxID+1;
                        newID = maxID;
                        parenthits(j)= hits(j);
                    else
                        % Detect drastic change in cell size
                        if numel(CC.PixelIdxList{j})<numel(nuc.area(hits(j),locframe-1))*0.5
                            % If yes, then add new cell
                            maxID=maxID+1;
                            newID = maxID;
                            parenthits(j)= hits(j);
                        else
                            newID = hits(j);
                            parenthits(j)= nuc.parent(hits(j),locframe-1);
                        end
                    end
                end
                % Set attributes to nuclei: nuc(newframe)
                nuc.frames(newID,locframe) = 1;
                nuc.x(newID,locframe) = stats(j).Centroid(1);
                nuc.y(newID,locframe) = stats(j).Centroid(2);
                nuc.ind(newID,locframe) = newID;
                nuc.parent(newID,locframe) = parenthits(j);
                nuc.daughter1(newID,locframe) = 0;
                nuc.daughter2(newID,locframe) = 0;
                nuc.area(newID,locframe) = numel(CC.PixelIdxList{j});
                nuc.ecc(newID,locframe) = stats(j).Eccentricity;
                nuc.ang(newID,locframe) = stats(j).Orientation;
                nuc.radius(newID,locframe) = sqrt(nuc.area(newID,locframe)/pi);
                bwtemp(CC.PixelIdxList{j}) =  nuc.ind(newID,locframe);
            end
            % Create new label
            BWL(:,:,locframe) =  bwtemp;
            nucActive{locframe} = find(nuc.frames(:,locframe)==1);
        end
        close(h);
        % Find the daughter of nuclei
        for j=1:size(nuc.frames,1)
            % Find the cell ID
            indj=max(nuc.ind(j,:));
            if indj>0
                % Find its daughters
                daughter_true=unique(nuc.ind((nuc.parent==indj)));
                while numel(daughter_true)<2
                    daughter_true=[daughter_true 0];
                end
                nuc.daughter1((nuc.ind==indj))=daughter_true(1);
                nuc.daughter2((nuc.ind==indj))=daughter_true(2);
            end
        end
        % Find the cell cycle
        try
            numcell=sum(nuc.frames,1);  % Extract cell number per frame
            max_mitosis=4;              % Maximum of frame that mitosis can exists
            change_cc=[(numcell(1+max_mitosis:end))./numcell(1:end-max_mitosis)];
            change_cc=(change_cc<10)&(change_cc>1.5); % Find big rational changes in nuclei number, avoid blinking frames
            for j=numel(change_cc):-1:max_mitosis+1 % Clean big chunk of increases
                if change_cc(j)
                    change_cc(j-1:-1:j-max_mitosis)=0;
                end
            end
            change_bw=cumsum([zeros(1,max_mitosis) change_cc])+Cycle0;
            change_cc=[(numcell(1+max_mitosis:end))./numcell(1:end-max_mitosis)];
            change_cc=(change_cc<10)&(change_cc>1.5); % Find big rational changes in nuclei number, avoid blinking frames
            for j=1:1:numel(change_cc)-max_mitosis % Clean big chunk of increases
                if change_cc(j)
                    change_cc(j+1:1:j+max_mitosis)=0;            
                end
            end
            change_fw=cumsum([change_cc zeros(1,max_mitosis)])+Cycle0;
            for i=1:size(nuc.frames,1)
                nuc.cycle(i,:)=(change_bw+change_fw)/2;
            end
        catch
            msgbox('Cannot find cell cycle: check segmentation consistency');
        end
        
        % THE TRACKING IS FINNISHED. EXPORT THE DATA
        k = strfind(FileName,'.');
        trackfile=[PathName,FileName(1:k(end)-1),'_Track_' num2str(Nf) '.mat'];
        save(trackfile,'nuc','BWL','-v7.3');
        msgbox('Tracking done. Please set nuclei cycle next');
    end

    function hchecktrack_Callback(~,~)
        if (size(nuc.frames,2)==Nf)&&(frame<Nf)
            % Create a panel to vísualize the process
            tempfig = figure('Name','Visualizing...');
            clf

            subplot(211);
            locframe = frame;       % Identifier for currently processed frame
            nucActive{locframe} = find(nuc.frames(:,locframe)==1);
            % Prepare to draw the figures
            RGB = label2rgb(BWL(:,:,locframe),'lines','k');
            imshow(RGB,[]);
            % Writeout the label
            texts = cat(1,num2str(nuc.ind(nucActive{locframe},locframe)));
                posal = cat(1,[nuc.x(nucActive{locframe},locframe)-10 nuc.y(nucActive{locframe},locframe)]);
            text(posal(:,1),posal(:,2),texts,'FontSize',10,'FontWeight','bold','Color','w','HorizontalAlignment','center')
            text(-10,-10,strcat('frame ',num2str(locframe)),'FontSize',10,'FontWeight','bold','Color','k')
            drawnow
            hold off;

            subplot(212);
            locframe = frame+1;       % Identifier for currently processed frame
            nucActive{locframe} = find(nuc.frames(:,locframe)==1);
            if locframe <= Nf
                % Prepare to draw the figures
                RGB = label2rgb(BWL(:,:,locframe),'lines','k');
                imshow(RGB,[]);
                % Writeout the label
                texts = cat(1,num2str(nuc.ind(nucActive{locframe},locframe)));
                posal = cat(1,[nuc.x(nucActive{locframe},locframe)-10 nuc.y(nucActive{locframe},locframe)]);
                text(posal(:,1),posal(:,2),texts,'FontSize',10,'FontWeight','bold','Color','w','HorizontalAlignment','center')
                text(-10,-10,strcat('frame ',num2str(locframe)),'FontSize',10,'FontWeight','bold','Color','k')
                drawnow
                hold off;
            end
        else
            msgbox('Complete tracking first');
        end
    end

    function hvisualizetrack_Callback(~,~)
       % Load track file
       try
            k = strfind(FileName,'.');
            trackfile=[PathName,FileName(1:k(end)-1),'_Track_' num2str(Nf) '.mat'];
            
            load(trackfile,'nuc','BWL');
            if size(nuc.frames,2)==Nf
                % Contruct the tree:
                tree_parent=[];
                tree_time=[];
                for i=1:size(nuc.frames,1)
                    ttmp=find(nuc.frames(i,:),1,'first');
                    idcell=nuc.ind(i,ttmp);
                    tree_time(idcell)=ttmp;
                    tree_parent(idcell)=nuc.parent(i,ttmp);
                end
                figure;
                draw_tree(tree_parent,tree_time,nuc.cycle(1,:));
            end
       catch
            msgbox('Complete tracking first');
       end
    end
    
    

    function hsetcycle_Callback(~,~)
       % Load track file
       try
            k = strfind(FileName,'.');
            trackfile=[PathName,FileName(1:k(end)-1),'_Track_' num2str(Nf) '.mat'];
            load(trackfile,'nuc','BWL');
            % Roughly assign frame to nuclei cycle

            cellcycle_range=[8:14];
            for i=1:numel(cellcycle_range)
                Outtext{i}=['Range of nuclear cycle ' num2str(cellcycle_range(i)) ' (from to)'];
            end
            inpNf=inputdlg(Outtext,'Roughly assign frames to nuclear cycle');
            
            nuc.cycle(:,:)=0;
            for i=1:numel(cellcycle_range)
                tmp = str2num(inpNf{i});
                if numel(tmp)==2
                    nuc.cycle(:,tmp(1):tmp(2))=cellcycle_range(i);
                    nc_info(i,:)=[cellcycle_range(i) tmp(1) tmp(2)];
                else
                    nc_info(i,:)=[cellcycle_range(i) 1 -1];
                end
            end
            save(trackfile,'nuc','BWL','nc_info','-v7.3');
       catch
            msgbox('Input Error/Incomplete tracking first');
       end
    end
%% Auxiliary functions
    function Iout=addMarker(Iin)
        % Add the marker at position frametop and framebottom.
        Marker = [0 0 0 0 0 0 0 0 0 0 1 1;
                  0 0 0 0 0 0 0 1 1 1 1 1;...
                  0 0 0 0 1 1 1 1 1 1 1 1;...
                  0 0 1 1 1 1 1 1 1 1 1 1;...
                  1 1 1 1 1 1 1 1 1 1 1 1;...
                  0 0 0 1 1 1 1 1 1 1 1 1;...
                  0 0 0 0 0 0 1 1 1 1 1 1;...
                  0 0 0 0 0 0 0 0 1 1 1 1;...
                  0 0 0 0 0 0 0 0 0 0 1 1;...
                  ];
        % Inverse the marker
        Marker=Marker(:,end:-1:1);
        % Create marker image
        Itmp=zeros(size(Iin));
        Itmp(frametop,8)=1;
        Itmp(framebottom,8)=1;
        Iout=(conv2(Itmp,Marker,'same')==0);
    end
    
    % Count the number of cell and orphaned cell
    function [all,wmother,orphaned]=cellcount(bwi,bwpre)
        % Check cell count
        all=unique(bwi);
        all=all(all>0);
        % Check overlapping cell
        bwi=double(bwi).*(bwpre>0);
        wmother=unique(bwi);
        wmother=wmother(wmother>0);
        % No mother cell:
        orphaned=setdiff(all,wmother);
    end

    % Brightness adjustment
    function hbrightup_Callback(~,~)
        brightness=brightness*1.5;
        showI;
    end
    function hbrightdown_Callback(~,~)
        brightness=brightness/1.5;
        showI
    end
    function hbrightreset_Callback(~,~)
        brightness=1;
        showI
    end
    % Update mask and image frames
    function showI(img,bwi)
        if ~exist('img','var')
            img=imadjust(uint8(Img(:,:,frame)));
            bwi=uint16(BW(:,:,frame));
        end
        
        img=uint8(img*brightness);
        imgMarker=addMarker(img);
        % For IMAGE panel
        axes(ha);
        if ~exist('bwi','var')
            bwi=zeros(size(img));
        end
        bwi_=bwi>0;
        % Select between visual mode
        switch showBW
            case 0          % Show image only
                imagesc(img,[0 256]);
            case 1          % Show image with mask
                locsize = size(img);
                ImBW = zeros(locsize(1),locsize(2),3);
                ImBW(:,:,1)=double(img);
                ImBW(:,:,2)=double(img);
                ImBW(:,:,3)=double(img);
                ImBW(:,:,1)=ImBW(:,:,1).*imgMarker;
                ImBW(:,:,2)=ImBW(:,:,2).*imcomplement(bwi_).*imgMarker;
                ImBW(:,:,3)=ImBW(:,:,3).*imcomplement(bwi_).*imgMarker;
                imagesc(uint8(ImBW),[0 256]);
            case 2          % Show image with mask contour
                locsize = size(img);
                ImBW = zeros(locsize(1),locsize(2),3);
                ImBW(:,:,1)=double(img);
                ImBW(:,:,2)=double(img);
                ImBW(:,:,3)=double(img);
                bwi_=conv2(double(bwi),[0 -1 0;-1 4 -1;0 -1 0],'same')~=0;
                ImBW(:,:,1)=ImBW(:,:,1).*imgMarker+ bwi_*255;
                ImBW(:,:,2)=ImBW(:,:,2).*imgMarker+ bwi_*255;
                ImBW(:,:,3)=ImBW(:,:,3).*imgMarker+ bwi_*255;
                imagesc(uint8(ImBW),[0 256]);
            case 3          % Previous (aligned) mask (red) with current mask contour
                if frame==1
                    imagesc(img,[0 256]);
                else
                    ref_frame=BW(:,:,frame-1); % Load previous mask
                    if numel(tform_rec{frame})
                        ref_frame=applyimwarp(ref_frame,tform_rec{frame},alignratio);
                    end
                    ref_frame=ref_frame>0;
                    % Find cells with no mother
                    [~,wmother,orphaned]=cellcount(BW(:,:,frame),ref_frame);
                    % Create image of cell with no mother:
                    bw_wmother=bwi.*uint16(ismember(BW(:,:,frame),wmother));
                    bw_wmother=conv2(double(bw_wmother),[0 -1 0;-1 4 -1;0 -1 0],'same')~=0;
                    % Create image of cell with mother:
                    bw_orphaned=bwi.*uint16(ismember(BW(:,:,frame),orphaned));
                    bw_orphaned=conv2(double(bw_orphaned),[0 -1 0;-1 4 -1;0 -1 0],'same')~=0;
                    
                    locsize = size(img);
                    ImBW = zeros(locsize(1),locsize(2),3);
                    ImBW(:,:,1)=double(img);
                    ImBW(:,:,2)=double(img);
                    ImBW(:,:,3)=double(img);
                    ImBW(:,:,1)=ImBW(:,:,1).*imgMarker;
                    ImBW(:,:,2)=ImBW(:,:,2).*imcomplement(ref_frame).*imgMarker;
                    ImBW(:,:,3)=ImBW(:,:,3).*imcomplement(ref_frame).*imgMarker;
                    ImBW(:,:,1)=ImBW(:,:,1).*imgMarker+ bw_wmother*255;
                    ImBW(:,:,2)=ImBW(:,:,2).*imgMarker+ bw_wmother*255 + bw_orphaned*255;
                    ImBW(:,:,3)=ImBW(:,:,3).*imgMarker+ bw_wmother*255 + bw_orphaned*255;
                    imagesc(uint8(ImBW),[0 255]);
                end
        end
        str=['Frame ',num2str(frame)];
        set(hmessages,'String',[str '. ' num2str(numel(unique(bwi(bwi(:)>0)))) ' cells']);
        % For MASK panel
        axes(hbw);
        bwodd=zeros(size(bwi));
        RGB = label2rgb(bwi,'lines','k');
        % If we need to check cell eccentricity
        
        if (get(hcheckecc,'Value') ==1)&&(numel(unique(bwi))>1)
            % Find odd looking cells
            CC = colorconncomp(bwi);
            stats = regionprops(CC,'Eccentricity','Area');
            sizerec = cellfun(@numel,CC.PixelIdxList);
            % Alert on cell too small, odd looking
            oddnuclei=find(([stats.Eccentricity]>0.8)|(sizerec<100)|([stats.Area]<median([stats.Area])*0.3));
            for i=oddnuclei
                bwodd(CC.PixelIdxList{i})=1e5;
            end
            if numel(oddnuclei)
                RGB=RGB*0.3+uint8(cat(3,bwodd,bwodd,bwodd));
            end
        end
        imshow(RGB);
    end

    function hzoom_Callback(~,~)
        %ha
        hapos=get(ha,'Position');
        hbwpos=get(hbw,'Position');
        ds=hbwpos(2)+hbwpos(4)-hapos(2); % Height of the two axes
        zoommode=mod(zoommode+1,3);
        zoomscope=1.5;
        switch zoommode
            case 0
                % Set equal axis
                hapos(3:4)=hapos(3:4)/(2-zoomscope);
                hbwpos(3:4)=hbwpos(3:4)/zoomscope; 
            case 1
                % Zoom in ha
                hapos(3:4)=hapos(3:4)*zoomscope;
                hbwpos(3:4)=hbwpos(3:4)*(2-zoomscope);
            case 2
                % Zoom in hbw
                hapos(3:4)=hapos(3:4)/zoomscope*(2-zoomscope);
                hbwpos(3:4)=hbwpos(3:4)*zoomscope/(2-zoomscope);
        end
        hbwpos(2)=hapos(2)+hapos(4);
        set(ha,'Position',hapos);
        set(hbw,'Position',hbwpos);
        %hbw
    end
    function showop(opimage)
        axes(hbw);
        imagesc(opimage,[0 1]);
    end

    function bwo = segment(bwi)
        bw=double(bwi);
        % Segmentation algorithm
        %l = graythresh(uint8(bw));
        %bw = im2bw(uint8(bw),l);
        D = bwdist(~bw);
        l = graythresh(uint8(D));
        bw = im2bw(uint8(D),l);
        % Smoothen the masks
        bw = imdilate(bw,strel('diamond',2));
        bw = imclose(bw,strel('disk',2));
        %bw = imfill(bw,'holes');
        % Remove small regions
        bwl = bwconncomp(bw);
        S = regionprops(bwl,'Area');
        Sm = mean([S.Area]);
        if ~isnan(Sm)
            switch get(hsegchoice,'value')
                case 2
                    bw = bwareaopen(bw,round(Sm/5));
                otherwise
                    bw = bwareaopen(bw,round(Sm/3));
            end
        end
        D = -bwdist(~bw);
        mask = imextendedmin(D,2);
        D2 = imimposemin(D,mask);
        LD = watershed(D2);
        bw(LD == 0) = 0;
        bwo = imopen(bw,strel('diamond',3));
        
        % Convert from binary image to something:
        bwo=bwlabel(bwo);
    end


    % Propagate to the prev frame
    function hpropminus_Callback(~,~)
        frame = frame - 1;
        if frame<1
            frame=1;
        end
        if frame>Nf
            frame = Nf;
        end
        set(hframe,'Value',frame);
        hfilter_Callback();
        hsegment_Callback();
        autosize();
    end

    % Propagate to the next frame
    function hpropplus_Callback(~,~)
        frame = frame + 1;
        if frame<1
            frame=1;
        end
        if frame>Nf
            frame = Nf;
        end
        set(hframe,'Value',frame);
        hfilter_Callback();
        hsegment_Callback();
        autosize();
    end

    function autosize()
        sizeadjust = str2double(get(hautosizef,'String'));
        if numel(sizeadjust)
            if sizeadjust>0
                for i=1:sizeadjust
                    hbwplus_Callback();
                end
            else
                for i=1:abs(sizeadjust)
                    hbwminus_Callback();
                end
            end
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HOTKEY HOTKEY HOTKEY HOTKEY HOTKEY HOTKEY HOTKEY HOTKEY HOTKEY %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function keypress(~, eventdata, ~)
        switch eventdata.Key
            case 'a'
                if numel(eventdata.Modifier)
                    haddelip_Callback;
                else
                    haddfree_Callback;
                end
            case 'r'
                hremove_Callback();
            case 's'
                hsegment_Callback();
            case 'x'
                hshow_Callback();
            case 'space'
                if numel(eventdata.Modifier)
                    hpropminus_Callback;
                else
                    hpropplus_Callback;
                end
            case 'o'
                hbwplus_Callback();
            case 'p'
                hbwminus_Callback();
            case 'f1'
                show_help();
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELP HELP HELP HELP HELP HELP HELP HELP HELP HELP HELP HELP HELP %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function show_help()
        Msg={'Welcome to LiveFly - Nuclei segmentation user interface',...
            '',...
            'Brief procedure:',...
            '      1.Load the maximum projection movies',...
            '      2.Select a proper filter size for the image',...
            '      3.Segment the frame: ',...
            '            + Use Automatic segmentation',...
            '            + Add or Remove the nuclei mask if needed',...
            '            + Adjust the size of the mask (with BW+, BW-) if needed',...
            '            + Move on to the next frame',...
            '      4.Track nuclei and lineages',...
            '      5.Plot lineages tree to check nuclei lineage and accuracy of division timing',...
            '      6.Readjust cell cycle duration if needed ',...
            '',...
            'Useful hotkeys:',...
            '      S: automatic segmentation of the current frame',...
            '      A: add a mask (using Zommed panel)',...
            '           A: Free hand drawing ',...
            '           Shfit+A: Draw an eclipse in selected retangular box',...
            '      R: remove a mask (using Zommed panel)',...
            '      O: increase mask size by 1',...
            '      P: decrease mask size by 1',...
            '      X: Switch visual mode: ',...
            '      Space: Mask the next frame using the same setting',...
            '      Shift+Space: Mask the previous frame using the same setting',...
            '',...
            'Hints',...
            '      -Always save the mask before Track',...
            '      -Backup the mask often (?)',...
            '      -When correcting a nuclei mask using Add, make sure to delete the mask first (using Remove)',...
            };
        msgbox(Msg,'Help','help')
    end
end
