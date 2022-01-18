addpath('..\..\Tool\bfmatlab\');
addpath('AddTextToImage\');
headfolder = 'Y:\equipe_dostatni\y_diaw\RAW';
add_timestamp = 1;
add_embryo_axis = 1;

margin_frame = 10;  % Add 10 frame before and after the interested phase
%% Specific params
mov_folder = '200129_YD1';
nc_range = [11 12 13];
export_folder = 'E:\Synthetic_paper_data\Z2B6';

rename_movie = 'Z2B6_movie1';

movie_label = 'Z2B6_movie1';

mkdir(export_folder);

%% Get file name and config data
listing = dir(fullfile(headfolder,mov_folder,'table_summary'));
for i=1:numel(listing)
    if strfind(listing(i).name,'config.m')
        ConfigName=listing(i).name;
    end
end
ConfigName
run(fullfile(headfolder,mov_folder,'table_summary',ConfigName));
filename=[ConfigName(1:end-9) '_fixed.txt'];
if ~exist(fullfile(headfolder,mov_folder,'table_summary',filename),'file')
    filename=[ConfigName(1:end-9) '.txt'];
end
%% Create interested frame:
matdata = dlmread(fullfile(headfolder,mov_folder,'table_summary',filename));
fr_range_=[];
nc_range_ = [];
for nc = nc_range
    fr_range_ = [fr_range_;matdata(matdata(:,12)==nc,2)];
    nc_range_ = [nc_range_;matdata(matdata(:,12)==nc,12)];
end
fr_range = [min(fr_range_)-15 max(fr_range_)+15];
fr_range(fr_range<=0)=0;
fr_range(fr_range>max(matdata(:,2)))=max(matdata(:,2));

ncfrom=zeros(1,1000);
ncto=zeros(1,1000);
for fr=fr_range_'
    nc_possible = unique(nc_range_(fr_range_==fr));
    if numel(nc_possible)==1
        ncfrom(fr) = nc_possible;
        ncto(fr) = 0;
    end
    if numel(nc_possible)==2
        ncfrom(fr) = min(nc_possible);
        ncto(fr) = max(nc_possible);
    end
end
%% Parameters for 3D projection (generally unchanged)
fsize=1;                 % Size of the filter before creating maximum projection, should be small (~1 or 2) for nc14.
max_projection=5;        % number of maximum layer for the projection - approximating half size of a nuclei in z stack
alignratio=0.3;          % Resize the images before alignment for speed
brightness=[1.5 1];      % Increase in brightness if needed for green (1) and red (2) channel
%% Create the reader for tif file
% Construct an empty Bio-Formats reader
reader = bfGetReader();
% Decorate the reader with the Memoizer wrapper
reader = loci.formats.Memoizer(reader);
reader.setId(fullfile(headfolder,mov_folder,main_mov));
%% Get metadata
omeMeta = reader.getMetadataStore();
z_max = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices

%% Create the maximum projection file
n_frame = fr_range(1):fr_range(2);
tic
for frame=n_frame
    display(frame)
    if frame==n_frame(1)
        writemode='overwrite';
    else
        writemode='append';
    end
    % Load the images
    Ifinal = [];
    I = [];
    for channel=[1 0]
        for z=1:z_max
            iPlane = reader.getIndex(z - 1, channel, frame - 1) + 1;
            Itmp = bfGetPlane(reader, iPlane)*brightness(channel+1);
            if fsize>1
                Itmp_=medfilt2(Itmp,[fsize fsize]);                
            else
                Itmp_=Itmp;
            end
            I(:,:,z)=round(Itmp_);
        end
    % Create projection %1
        %[I,~]=sort(I,3,'descend');
        %Iout=mean(I(:,:,1:max_projection),3);
    % Create projection %2
        Iout=max(I,[],3);
        % Save projection
        Iout=uint8(Iout);
        Ifinal = cat(3,Ifinal,Iout);
    end
    Ifinal = cat(3,Ifinal,Iout*0);
    % Make label frame:
    Isize = size(Iout);
    t = (frame - n_frame(1))*dt;
    Ifinal = AddTextToImage(Ifinal,[movie_label '. T=' num2str(round(t),'%d') 's'],[20 Isize(2)-270],[1 1 1],'FontSize',20);
    if ncfrom(frame)
        %if ncto(frame)
        %    Ifinal = AddTextToImage(Ifinal,['nc' num2str(ncfrom(frame)) '-' 'nc' num2str(ncto(frame))],[40 Isize(2)-250],[1 1 1],'FontSize',20);
        %else
        Ifinal = AddTextToImage(Ifinal,['nc' num2str(ncfrom(frame))],[45 Isize(2)-270],[1 1 1],'FontSize',20);
        %end
    end
    % Draw AP axis:
    for i=-5:1:5
        lwd = 1;
        if i==0
            llen=20;
        else
            llen=10;
        end
        APlen = shift_left+shift_right+Isize(2);
        y0 = APlen/2 - shift_left;
        ypos = round(y0 + i*APlen/10);
        if (ypos>1)&&(ypos<Isize(2)-1)
            Ifinal(Isize(1)-(0:llen-1),ypos-lwd:ypos+lwd,1:3)=255;
        end
    end
    % Save all
    imwrite(Ifinal,fullfile(export_folder,[rename_movie '.tif']),'WriteMode',writemode);
end
display('Done');
tic