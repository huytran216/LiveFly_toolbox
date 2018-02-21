addpath('..\Tool\bfmatlab\');

%% Input folder
folder='C:\Matlab\Test\';     % Original movie filder
mov_in='RAWMovie.tif';              % Original movie name with extension
mov_out='RED';                      % Maximum projection movie [NO EXTENSION]

% Additional parameters
channel=1;               % Channel for Nuclei image (0 or 1)
n_frame=[1:220];         % Range of frame to be analyzed
fsize=2;                 % Size of the filter before creating maximum projection, should be small (~1 or 2) for nc14.
quantile_threshold=[0.5];% Quantile 
max_projection=5;        % number of maximum layer for the projection - approximating half size of a nuclei in z stack
alignratio=0.3;          % Resize the images before alignment for speed
brightness=1.2;          % Increase in brightness if needed
%% Create the reader for tif file
% Construct an empty Bio-Formats reader
reader = bfGetReader();
% Decorate the reader with the Memoizer wrapper
reader = loci.formats.Memoizer(reader);
reader.setId(fullfile(folder,mov_in));
%% Get metadata
omeMeta = reader.getMetadataStore();
z_max = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
%% Create the maximum projection file
tform_rec=cell(1,max(n_frame));    % Storage for transformation
I=[];
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
                Itmp_=zeros(size(Itmp));
                for j=1:numel(quantile_threshold)
                    Itmp_=Itmp_+double(ordfilt2(Itmp,round(quantile_threshold(j)*(fsize^2)),true(fsize)));
                end
            else
                Itmp_=Itmp;
            end
            I(:,:,z)=round(Itmp_/numel(quantile_threshold));
        end
    % Create projection %1
        %[I,~]=sort(I,3,'descend');
        %Iout=mean(I(:,:,1:max_projection),3);
    % Create projection %2
        I = convn(I,ones(1,1,max_projection)/max_projection,'same');
        Iout=max(I,[],3);
        % Save projection
        Iout=uint8(Iout);
        imwrite(Iout,fullfile(folder,[mov_out '.tif']),'WriteMode',writemode);
    % Calculate alignment:
        tform_rec{frame}=[];
        Iout=imadjust(imresize(Iout,alignratio));
        
        if frame>n_frame(1)
        % Calculate the geometric transformation
            % Normal fit
            [tform,~]=imregdemons(Ipre,Iout,[50 50 50 50 50],'AccumulatedFieldSmoothing',2,'DisplayWaitBar',false,'PyramidLevels',5);
            % Convert to polynomial surface projection
            [x,y]=meshgrid(1:size(Iout,2),1:size(Iout,1));
            x_=tform(:,:,1);y_=tform(:,:,2);
            tform_poly=struct;
            tform_poly.x = fit([x(:),y(:)],x_(:),'poly44');
            tform_poly.y = fit([x(:),y(:)],y_(:),'poly44');
            tform_rec{frame}=tform_poly;
        % Plot verification - Compare
            if numel(n_frame)<3
                % Create checkerboard matrix
                sqrsize=20;
                Isqr=checkerboard(sqrsize,ceil(size(Iout,1)/sqrsize),ceil(size(Iout,2)/sqrsize));
                Isqr=Isqr(1:size(Iout,1),1:size(Iout,2));

                Isqr=Ipre;
                % Original set
                    subplot(311);imshowpair(uint8(Isqr),uint8(Iout));
                    title('Original');
                % Test poly method
                [Itmp,rd]=imwarp_poly(Isqr,tform_poly);
                    Inew=resetframe(Itmp,rd,size(Ipre));
                    subplot(313);imshowpair(uint8(Isqr),uint8(Inew));
                    title('Poly');
                % Test daemon method
                [Itmp,rd]=imwarp(Isqr,tform);
                    Inew=resetframe(Itmp,rd,size(Ipre));
                    subplot(312);imshowpair(uint8(Isqr),uint8(Inew));
                    title('Daemon');
            end
        end
        Ipre=Iout;
end
save(fullfile(folder,[mov_out '_align.mat']),'tform_rec','alignratio','-v7.3');