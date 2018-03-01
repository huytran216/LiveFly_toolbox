clear all; %clear the memory from all the variables
close all; %close all the open figures
warning off;
clc;       %reset the writing space
addpath('..\Tool\bfmatlab\');
%%  LOAD THE MULTI-TIFF FILE FOR THE NUP CHANNEL
mov_folder='D:\Matlab\HUY TESTBOX FOR LIVEFLY\Test';                       %indicate the full path movies file
main_mov='RAWMovie.tif';                            %title of ensemble movie (with extension)
nuclei_mov='RED';                                   %title of segmented nuclei movie (with/without extension)

%frame for segmentation file
seg_init = 1;
seg_end  = 220;
%select the frame range you want to analyze
it_start=1;
it_end=220;
%% Parameters:
ConfigName=0;
if strcmp(questdlg('Choose a configuration file?'),'Yes')
    [ConfigName,PathName,FilterIndex] = uigetfile(fullfile(mov_folder,'*.m'),'Select the configuration file');
    run(fullfile(PathName,ConfigName));
else
    % --------MODIFY THE PARAMETERS HERE----------
    
    % Movie frame and resolution parameters
    shift_left=786;         % Length of the cat portion on the left
    shift_right=805;        % Length of the cat portion on the right
    dt=14.73;                % Time step in seconds that you read from the metafile
    A_pole=1;               % Anterior pole position 1=left
    channel=0;              % Channel for Nuclei image (0 or 1)
    x_resolution=0.171;     % Pixel size (micron/pixels) on XY axis (optional)
    z_resolution=0.5;       % Pixel size (micron/pixels) on Z axis

    % IMPORTANT PARAMETERES FOR SPOT DETECTION %% 
    th1=57;                 % Minimal threshold for spot detection, set to 0 for debug mode
    th2=29;                 % Minimal threshold for median filter, set to 0 for debug mode
    th=[th1 th2];
    % Spot detection parameters:
    averaging_radius=3;     % there is an average filter for the images. averaging_radius=1 means no average. About the spot size
    voxels_min=9;          % minimal number of voxels for a spot
    voxels_max=60;          % Maximum number of voxels for a spot
    fact_r=1.2;             % tolerance radius (if fact_r=1) the distance is equal to the radius. fact_r for histone is generally bigger than 1.
end
%% Ask some question about debug mode
if strcmp(questdlg('Enter debug mode to find th1, th2?'),'Yes')
    th=[0 0];
    answer=inputdlg('Enter a frame number with spots');
    it_start=str2double(answer);
    it_end=it_start;
end
%% Brightness adjustment after spot detection
brightness_ratio=[1.2 1.2]; % Amplify the (1) green and (2) red intensity for better visualization
%% Load basic information on movies
[~,nuclei_mov,~] = fileparts(nuclei_mov);
ms2_mov = fullfile(mov_folder,main_mov);
seg_mov = fullfile(mov_folder,nuclei_mov);

I=imread([seg_mov '.tif'],1);         %read the frame 1

Lx=size(I,1);
Ly=size(I,2);
%% Load the segmentation image
load([seg_mov '_Track_' num2str(seg_end) '.mat'],'nuc','BWL');

%% Create the reader for tif file
% Construct an empty Bio-Formats reader
reader = bfGetReader();
% Decorate the reader with the Memoizer wrapper
reader = loci.formats.Memoizer(reader);
reader.setId(ms2_mov);
%% Get metadata
omeMeta = reader.getMetadataStore();
z_max = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
%% SPOTS INFORMATIONS (compulsory)
%load the movie and store some variables regarding the images
Lx2=Lx;
Ly2=Ly;

window_around_spot=30;  %pixels, this is the window around the spot to do the gaussian fit, 30 pixels should be fine for the usual resolution 2048x2048

if it_end>seg_end   % If spot movie is longer than segmentation movie 
    it_end=seg_end;
end

Ispot_new=zeros(Lx2,Ly2,it_end-it_start+1);
Ispot_raw=zeros(Lx2,Ly2,it_end-it_start+1);
it14=0;

tic
for it=it_start:it_end
    
    s_start=(it-1)*z_max+1;
    s_end=(it-1)*z_max+z_max;    
    
    fprintf(1,'time %d\n',it);
 
    [spot,Ispot,raw_spot]=find_ms2_spots_1spot(reader,channel,it,nuc,th,z_max,voxels_min,voxels_max,fact_r,window_around_spot,averaging_radius);

    %figure,imshow(Ispot);
    
    fprintf(1,'\n');
    Ispot_new(:,:,it)=Ispot;
    raw_spot_t(:,:,it)=raw_spot;
    raw_nuclei_t(:,:,it)=BWL(:,:,it);
    infos_spot{it}=spot;
end
toc
%% Export the movies
mkdir(fullfile(mov_folder,'output_images'));
for it=it_start:it_end
    if it==it_start
        writemod='overwrite';
    else
        writemod='append';
    end
    I_mask=imread([seg_mov '.tif'],it);
    % Combined the two channels RED and GREEN
    I_total=cat(3,uint8(brightness_ratio(1)*I_mask),uint8(brightness_ratio(2)*raw_spot_t(:,:,it)),uint8(zeros(size(I_mask))));
    imwrite(I_total,fullfile(mov_folder,'output_images','img_combined.tif'),'tif','compression','none','writemode',writemod);
    % Create label on image:
    ax=imshow(zeros([Lx Ly]),'InitialMagnification',100);
    proportionXY=get(gcf,'Position').*get(gca,'Position');proportionXY=proportionXY([4 3])./[Lx Ly];
    for i=1:numel(nuc.frames(:,it))
        if nuc.frames(i,it)
            postmp=[nuc.x(i,it) nuc.y(i,it)];
            labeltmp=num2str(nuc.ind(i,it));
            text(gca,postmp(1)-6,postmp(2),labeltmp,'Color','white','FontSize',8);
        end
    end
    F = getframe(gca);
    I_lb = F.cdata;I_lb=I_lb(:,:,1);
    I_lb = imresize(I_lb,[Lx Ly],'method','bilinear');
    cla;
    Iadd=(raw_nuclei_t(:,:,it)>=1)*0.3 +(I_lb>100)*0.3;
    cla;
    %save the segmented images in the folder segmented_spots
    I_res=double(Ispot_new(:,:,it)+Iadd);
    imwrite([double(I_total)/256;cat(3,I_res,I_res,I_res)],fullfile(mov_folder,'output_images','img_check.tif'),'tif','compression','none','writemode',writemod); %save the segmented images in the folder segmented_spots
end

%% connect and allocate the infos regarding nuclei and spots (compulsory)

id_max=max(nuc.ind(:));

spots_number=1;   % maximal number of spots allowed
n_spots_infos=11; % number of variables for each spot

infos_spot_nuc=zeros(id_max,it_end,n_spots_infos,spots_number);

for it=it_start:it_end
    if ~isempty(infos_spot{it})
        for in=1:size(nuc.frames,1)
            id=nuc.ind(in,it);
            count_spot=0;
            
            for iis=1:size(infos_spot{it}.id_n,2)  %%%loop on the spots
                if infos_spot{it}.id_n(iis)==id    %%%to identify the nuclei
                    count_spot=count_spot+1;
                    ind_spot(count_spot)=iis;
                end
            end
            
            if count_spot>0
                ord_spot=zeros(count_spot,1);
                ord_spot_I=sort([infos_spot{it}.I(ind_spot(:))],'descend');
                ord_spot_A=sort([infos_spot{it}.A(ind_spot(:))],'descend');
                
                for ic=1:count_spot
                    ord_spot(ic)=find(infos_spot{it}.A(:)==ord_spot_A(ic));
                end
                
                for i=1:count_spot
                    for j=1:count_spot
                        jj=ord_spot(i);
                        infos_spot_nuc(id,it,1,i)=infos_spot{it}.x(jj);
                        infos_spot_nuc(id,it,2,i)=infos_spot{it}.y(jj);
                        infos_spot_nuc(id,it,3,i)=infos_spot{it}.z(jj);
                        infos_spot_nuc(id,it,4,i)=infos_spot{it}.size(jj);
                        infos_spot_nuc(id,it,5,i)=infos_spot{it}.I(jj);
                        infos_spot_nuc(id,it,6,i)=infos_spot{it}.A(jj);
                        infos_spot_nuc(id,it,7,i)=infos_spot{it}.ssx(jj);
                        infos_spot_nuc(id,it,8,i)=infos_spot{it}.ssy(jj);
                        infos_spot_nuc(id,it,9,i)=infos_spot{it}.I2d(jj);
                        infos_spot_nuc(id,it,10,i)=infos_spot{it}.bkg(jj);
                        infos_spot_nuc(id,it,11,i)=infos_spot{it}.resid(jj);
                    end
                end
            end
        end
    end
end
%% print a table with all the infos for each nucleus in each time step (compulsory)
mkdir(fullfile(mov_folder,'table_summary'));

filename = fullfile(mov_folder,'table_summary',['th' num2str(th1) '_' num2str(th2) '.txt']);
configname = fullfile(mov_folder,'table_summary',['th' num2str(th1) '_' num2str(th2) '_config.m']);

fp = fopen(filename,'w+');

id_max=max(nuc.ind(:));
% ir: index of nuclei
% id: id of nuclei
% it: time

for ir=1:size(nuc.frames,1)
    for it=it_start:it_end 
        id=nuc.ind(ir,it);
        if nuc.frames(ir,it)
            %if id<=size(v_frame(it).vs.ind,2) %in the case a nucleus was lost
                % Print the data
                daughter=[nuc.daughter1(ir,it) nuc.daughter2(ir,it)];
                daughter_=zeros(1,2);
                for i=1:numel(daughter)
                    daughter_(i)=daughter(i);
                end
                %fprintf(1,'%d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g ',id,it,nuc.x(ir,it),nuc.y(ir,it),nuc.radius(ir,it),nuc.ecc(ir,it),nuc.ecc(ir,it),nuc.ind(ir,it),nuc.parent(ir,it),daughter_(1),daughter_(2),nuc.cycle(ir,it),dt,Lx2,Ly2,shift_left,shift_right,th(end),A_pole);
                fprintf(fp,'%d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g ',id,it,nuc.x(ir,it),nuc.y(ir,it),nuc.radius(ir,it),nuc.ecc(ir,it),nuc.ecc(ir,it),nuc.ind(ir,it),nuc.parent(ir,it),daughter_(1),daughter_(2),nuc.cycle(ir,it),dt,Lx2,Ly2,shift_left,shift_right,th(end),A_pole);
                for j=1:spots_number
                    for i=1:n_spots_infos
                        %fprintf(1,'%g ',infos_spot_nuc(id,it,i,j));
                        fprintf(fp,'%g ',infos_spot_nuc(id,it,i,j));
                    end
                end                
                %fprintf(1,'\n');
                fprintf(fp,'\n');
            %end
        end
    end
end
fclose(fp);
%% Choose the config file:
if ~ConfigName
    if (th(1)>0)&&(th(2)>0)
        A={};   % Create introduction cell
        A_{1}=['Config file created on ' date];
        A_{2}='THIS IS THE CONFIGURATION FILE FOR SPOT DETECTION PROGRAM: ';
        A_{3}='If you need to reanalyze the data, change the params here';
        A_{4}='Then rerun the Spot detection with this file loaded';
        matlab.io.saveVariablesToScript(configname,{'A_','main_mov','nuclei_mov','shift_left','shift_right','dt','A_pole','channel','x_resolution','z_resolution','th1','th2','th','averaging_radius','voxels_min','voxels_max','fact_r'});
    end
end