function findzposition(folder,main_mov,datamat,outfile,zratio)
    zPosColumn = 6;     % Replacing the eccentricity field
    findBGIntensity = true;
        
    %% Begin analyzing
    n_frame = unique(datamat(:,2))';
    channel=1;               % Channel for Nuclei image (0 or 1)
    idcs   = strfind(folder(1:end-1),filesep);
    upperfolder = folder(1:idcs(end)-1);
    %% Parameters for 3D projection (generally unchanged)
    fsize=1;                 % Size of the filter before creating maximum projection, should be small (~1 or 2) for nc14.
    quantile_threshold=[0.5];% Quantile 
    max_projection=5;        % number of maximum layer for the projection - approximating half size of a nuclei in z stack
    brightness=1.2;          % Increase in brightness if needed
    rim_size=10;             % Margin (x,y scale) 
    %% Create the reader for tif file
    % Construct an empty Bio-Formats reader
    reader = bfGetReader();
    % Decorate the reader with the Memoizer wrapper
    reader = loci.formats.Memoizer(reader);
    reader.setId(fullfile(upperfolder,main_mov));
    %% Get metadata
    omeMeta = reader.getMetadataStore();
    z_max = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
    y_max = omeMeta.getPixelsSizeX(0).getValue(); % number of X slices
    x_max = omeMeta.getPixelsSizeY(0).getValue(); % number of X slices
    %% Create the maximum projection file
    % Record for spot intensity and background intensity
    Ispot_rec = zeros(1,1);
    Ibg_rec = zeros(1,2);
    I=zeros(x_max,y_max,z_max);
    Ispot=zeros(x_max,y_max,z_max);
    % Begin scanning
    for frame=n_frame
        display([num2str(frame) '/' num2str(max(n_frame))])
        % Load the images
            for z=1:z_max
                % Load nuclei channel
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
                if findBGIntensity
                    iPlane = reader.getIndex(z - 1, 1-channel, frame - 1) + 1;
                    Itmp = bfGetPlane(reader, iPlane)*brightness;
                    Ispot(:,:,z) = Itmp;
                end
            end
        % Create projection %1
            I = convn(I,ones(1,1,max_projection)/max_projection,'same');
            [Lx,Ly,~] = size(I);

            matid_range = find(datamat(:,2)==frame);

            for id=matid_range'
                datamat(id,zPosColumn)=-1;
                % get nuclei position:
                x = round(datamat(id,4));
                y = round(datamat(id,3));
                s = datamat(id,5);

                if (x-s-rim_size>1)&(x+s+rim_size<=Lx-1)&(y-s-rim_size>1)&(y+s+rim_size<=Ly-1)
                    % Find the nuclei position in z plan
                    Ilocal = I(round(x-s-rim_size:x+s+rim_size),round(y-s-rim_size:y+s+rim_size),:);
                    % estimating z position:
                    Imax= max(Ilocal,[],3); thresh = median(Imax(:));
                    Ilocal_ = Ilocal > thresh;
                    Imean = mean(mean(Ilocal_,1),2);
                    w = gausswin(s*2/zratio);
                    Imean = conv(Imean(:),w,'same');
                    [~,z] = max(Imean(:));
                    % Get the MCP-GFP intensity inside nuclei
                    if findBGIntensity
                        if (z>2)&&(z<z_max-1)
                            if datamat(id,28)>0
                                'spot here';
                                Ispotexact = Ispot(round(datamat(id,21)),round(datamat(id,20)),round(datamat(id,22)));
                                Ispot_rec = [Ispot_rec;Ispotexact];
                                if Ispotexact==255
                                    'alert'
                                end
                            else
                                'no spot here';
                                Ispotlocal = Ispot(round(x-s-rim_size/2:x+s+rim_size/2),round(y-s-rim_size/2:y+s+rim_size/2),z-1:z+1);
                                Ibg_rec = [Ibg_rec;mean(Ispotlocal(:)) sqrt(var(Ispotlocal(:)))];
                            end
                        end
                    end
    %                 Check for specific nuclei (if you know the id)                
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
                    
                    % Update the z position
                    datamat(id,zPosColumn)=z_max*1000+z;
                 end
            end
            %plot3(datamat(matid_range,3),datamat(matid_range,4),datamat(matid_range,zPosColumn),'o');
    end
    dlmwrite(outfile,datamat,'precision',6,'Delimiter','\t');
    if findBGIntensity
        save([outfile '_tmp.mat'],'Ispot_rec','Ibg_rec');
    end
end