% Follow spots with specific ID:
DatasetName = 'B6-near';
id = 7124;
addpath('..\..\Tool\bfmatlab\');
addpath('AddTextToImage\');
addpath '../../Data analysis';
dttmp = load(fullfile('../../Data analysis/final_dataset/',DatasetName));

export_movie = ['spot_' DatasetName '_' num2str(id) '.tif'];
%% Get params
datamat = dttmp.datamat(id);
tscnt = datamat.tscnt;
mov_folder = fullfile(alternative_path(dttmp.DatasetList(tscnt).Path),'..');

%% Get file name and config data
listing = dir(fullfile(mov_folder,'table_summary'));
for i=1:numel(listing)
    if strfind(listing(i).name,'config.m')
        ConfigName=listing(i).name;
    end
end
ConfigName
run(fullfile(mov_folder,'table_summary',ConfigName));
filename=[ConfigName(1:end-9) '_fixed.txt'];
if ~exist(fullfile(mov_folder,'table_summary',filename),'file')
    filename=[ConfigName(1:end-9) '.txt'];
end
%% Parameters for 3D projection (generally unchanged)
fsize=3;                 % Size of the filter before creating maximum projection, should be small (~1 or 2) for nc14.
max_projection=5;        % number of maximum layer for the projection - approximating half size of a nuclei in z stack
alignratio=0.3;          % Resize the images before alignment for speed
brightness=[1.5 1.5];      % Increase in brightness if needed for green (1) and red (2) channel
img_size = 40;           % Make image of img_size x img_size pixel
%% Create the reader for tif file
% Construct an empty Bio-Formats reader
reader = bfGetReader();
% Decorate the reader with the Memoizer wrapper
reader = loci.formats.Memoizer(reader);
reader.setId(fullfile(mov_folder,main_mov));
%% Get metadata
omeMeta = reader.getMetadataStore();
z_max = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices

%% Create the maximum projection file
n_frame = sort(round(datamat.time/dt));
tic
cntfr = 0;
for cntfr=1:numel(n_frame)
    frame = n_frame(cntfr);
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
            xtmp = datamat.xrec(cntfr);
            ytmp = datamat.yrec(cntfr);
            Itmp = Itmp(round(ytmp)+[-img_size/2:img_size/2],round(xtmp)+[-img_size/2:img_size/2]);
            if (fsize>1)&(channel>0)
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
        % Save projection . color
        Iout=uint8(Iout);
        Ifinal = cat(3,Ifinal,Iout);
    end
    Ifinal = cat(3,Ifinal,Iout*0);
    Ifinal_ = cat(3,[Ifinal(:,:,1);Ifinal(:,:,3)],...
        [Ifinal(:,:,3);Ifinal(:,:,2)],...
        [Ifinal(:,:,3);Ifinal(:,:,3)]);
    % Make label frame:
    Isize = size(Iout);
    t = (frame - n_frame(1))*dt;
    %Ifinal = AddTextToImage(Ifinal,['t=' num2str(round(t),'%d') 's'],[20 Isize(2)-270],[1 1 1],'FontSize',20);
    
    % Save all
    imwrite(Ifinal_,export_movie,'WriteMode',writemode);
end
display('Done');
tic
%% Plot data
plot(datamat.Adjustedtime - datamat.Adjustedtime(1),datamat.AdjustedIntensity);
xlim([0 1100]);
set(gca,'XTick',[0 200 400 600 800 1000]);