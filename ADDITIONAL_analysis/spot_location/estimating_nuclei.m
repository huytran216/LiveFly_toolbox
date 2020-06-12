addpath('..\..\Tool\bfmatlab\');

%% Input folder
folder='\\zserver\umr3664\equipe_dostatni\c_perez-romero\RAW\';           % Original movie folder
TSname = '170525_CP1';
TStype ='H6B6';
mov_in='H6B6HisRFPMCPnoNLS1.czi';               % Original movie name with extension

channel=1;               % Channel for Nuclei image (0 or 1)
n_frame=[1:382];         % Range of frame to be analyzed
%n_frame=330;
%% Parameters for 3D projection (generally unchanged)
fsize=3;                 % Size of the filter before creating maximum projection, should be small (~1 or 2) for nc14.
quantile_threshold=[0.5];% Quantile 
max_projection=5;        % number of maximum layer for the projection - approximating half size of a nuclei in z stack
alignratio=0.3;          % Resize the images before alignment for speed
brightness=1.2;          % Increase in brightness if needed
%% Create the reader for tif file
% Construct an empty Bio-Formats reader
reader = bfGetReader();
% Decorate the reader with the Memoizer wrapper
reader = loci.formats.Memoizer(reader);
reader.setId(fullfile(folder,TSname,mov_in));
%% Get the data file
run(fullfile(folder,TSname,'/table_summary/correction.m'));
matfile = dlmread(fullfile(folder,TSname,'table_summary',filename1),'\t');
%run(fullfile(folder,'/table_summary/',[filename1(1:end-10) '_config.m' ]));
[x_resolution,z_resolution] = get_resolution_config(['_' TSname '.txt']);
zratio = z_resolution/x_resolution;
%% Get metadata
omeMeta = reader.getMetadataStore();
z_max = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
%% Create the maximum projection file
tform_rec=cell(1,max(n_frame));    % Storage for transformation
I=[];
Ispot=[];
for frame=n_frame
    display(frame)
    if frame==n_frame(1)
        writemode='overwrite';
    else
        writemode='append';
    end
    % Load the images
        for z=1:z_max
            iPlane = reader.getIndex(z - 1, channel, frame - 1) + 1;
            Itmp = bfGetPlane(reader, iPlane)*brightness;
            if fsize>1
                Itmp_=medfilt2(Itmp,[fsize fsize]);
                %for j=1:numel(quantile_threshold)
                %    Itmp_=Itmp_+double(ordfilt2(Itmp,round(quantile_threshold(j)*(fsize^2)),true(fsize)));
                %end
            else
                Itmp_=Itmp;
            end
            I(:,:,z)=round(Itmp_/numel(quantile_threshold));
            
            iPlane = reader.getIndex(z - 1, 1-channel, frame - 1) + 1;
            Itmp = bfGetPlane(reader, iPlane)*brightness;
            Ispot(:,:,z) = Itmp;
        end
    % Create projection %1
        I = convn(I,ones(1,1,max_projection)/max_projection,'same');
        [Lx,Ly,Lz] = size(I);
        
        matid_range = find(matfile(:,2)==frame);
        
        for id=matid_range'
            matfile(id,32)=-1;
            % get nuclei position:
            x = round(matfile(id,4));
            y = round(matfile(id,3));
            s = matfile(id,5);
            rim_size=10;
            if (x-s-rim_size>1)&(x+s+rim_size<=Lx-1)&(y-s-rim_size>1)&(y+s+rim_size<=Ly-1)
                % Find the nuclei position in z plan
                Ilocal = I(round(x-s-rim_size:x+s+rim_size),round(y-s-rim_size:y+s+rim_size),:);
                Ispotlocal = Ispot(round(x-s-rim_size:x+s+rim_size),round(y-s-rim_size:y+s+rim_size),:);
                % estimating z position:
                Imax= max(Ilocal,[],3); thresh = median(Imax(:));
                Ilocal_ = Ilocal > thresh;
                Imean = mean(mean(Ilocal_,1),2);
                w = gausswin(s*2/zratio);
                Imean = conv(Imean(:),w,'same');
                [tmp,z] = max(Imean(:));
                zspot = matfile(id,22);
%                 if (zspot)&&(id==10504)
%                     subplot(321);
%                     plot(Imean(:));hold on;
%                     plot([z-s/zratio z+s/zratio],[tmp tmp],'LineWidth',2);
%                     hold off;
%                     subplot(322)
%                     imagesc(Ispotlocal(:,:,zspot));
%                     title(num2str(zspot));
%                     subplot(323);
%                     imagesc(Ilocal_(:,:,z));
%                     title(num2str(z));
%                     subplot(324);
%                     imagesc(Ilocal_(:,:,zspot));
%                     title(num2str(zspot));
%                     subplot(325);
%                     imagesc(Ilocal(:,:,z));
%                     title(num2str(z));
%                     subplot(326)
%                     imagesc(Ilocal(:,:,zspot));
%                     title(num2str(zspot));
%                     if abs(zspot-z)>10
%                         pause;
%                         id
%                     end
%                 end
                % Find right_size
                %if ((z-s*zratio-1)>0)&&((z+s*zratio+1)<=Lz)
                %    zbt= floor(z-s*zratio-1);
                %    Ilocal = Ilocal(:,:,zbt+1:ceil(z+z*zratio+1));
                %end
                matfile(id,32)=z;
                
             end
        end
        plot3(matfile(matid_range,3),matfile(matid_range,4),matfile(matid_range,32),'o');
end
%% save
dlmwrite([TStype '_' TSname '.txt'],matfile);