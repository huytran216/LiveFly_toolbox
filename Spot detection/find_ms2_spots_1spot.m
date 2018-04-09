function [detected_spot,Ispot,raw2d]=find_ms2_spots_1spot(reader,channel,it,nuc,th,z_max,voxels_min,voxels_max,fact_r,window,averaging_radius)
% Replace v_frame with nuc, at frame it.
% The use of r (rnuc) should be obsolete
% Input:
%       -reader: reader handle of the main movie file
%       -it: frame to analyze
%       -nuc: nuclei information (from Mathieu tool)
%       -th: threshold for spot detection
%       -z_max: number of z stack
%       -voxels_min: minimum spot volume
%       -fact_r: tolerence of spot position (outside nuclei mask)
%       -window: size of the cube around the spot to do gaussian fitting
%       -avaraging radius:

% Read the first frame and extract the image size
iPlane = reader.getIndex(0, 0, 0) + 1;
II = bfGetPlane(reader, iPlane);

Lx=size(II,1);
Ly=size(II,2);
I3=zeros(size(II,1),size(II,2),z_max);
raw3=zeros(size(II,1),size(II,2),z_max);
filtered3=zeros(size(II,1),size(II,2),z_max);

Ispot=zeros(size(II,1),size(II,2));


clear detected_spot;
detected_spot=struct('id_n',[],'id_s',[],'x',[],'y',[],'z',[],'size',[],'I',[],'A',[],'ssx',[],'ssy',[],'I2d',[],'bkg',[],'resid',[]);
%% Load all 3D image in the frame
    Irec=zeros(Lx,Ly,z_max);
    for zs=1:z_max  %loop on z
        iPlane = reader.getIndex(zs-1, channel, it-1) + 1;
        Irec(:,:,zs) = bfGetPlane(reader, iPlane);
    end
    Imax=max(Irec,[],3);
    h = fspecial('average', averaging_radius);
    Imax=imfilter(Imax,h);
    Imax=Imax-medfilt2(Imax, [averaging_radius*5 averaging_radius*5]);% Background intensity
%% Try spot detection on the 3D
    for zs=1:z_max  %loop on z
        I = Irec(:,:,zs);


        F=imfilter(I,h);
        
        if exist('ax1','var')
            xl=get(gca,'xlim');
            yl=get(gca,'ylim');
        end

        % Use a median filter
        Fbg=medfilt2(F, [averaging_radius*5 averaging_radius*5]);% Background intensity
        F_ = (F-Fbg);
        raw3(:,:,zs)=F_;
        filtered3(:,:,zs)=F_;
        
        if any(th==0)
            ax1=subplot(3,2,1:2);
            imagesc(F);
            colorbar;
            title(['Original image. z=' num2str(zs)]);
            ax2=subplot(3,2,3:4);
            imagesc(F_);
            colorbar;
            title(['Filtered image. z=' num2str(zs)]);
            ax3=subplot(3,2,5:6);
            imagesc(Imax);
            colorbar;
            title('Maximum projection');
            linkaxes([ax1,ax2,ax3],'xy');
            if exist('yl','var')
                xlim(xl);
                ylim(yl);
            end
            pause;
            I3(:,:,zs)=F*0;
        else
            I3(:,:,zs)=(F>th(1))&(F_>th(2));
        end
    end
    clear Irec;
%% Labeling the spot intensity
    L3D=bwlabeln(I3,8);                         % Label the spots 3D image
    loc3 = regionprops(L3D,'Area','Centroid');  % Extract the region property
    no=max(max(max(L3D)));                      % Number of spot detected
    raw2d=max(raw3,[],3);
    if numel(loc3)>500
        display([num2str(numel(loc3)) ' spot detected. Considering increasing threshold value']);
    end
%% Eliminate the clusters too far from the nuclei centres
    inner_clusters=[];
    d=1000*ones(no,size(nuc.frames,1));
    for i=1:no
        goodspot=false;
        for j=1:size(nuc.frames,1)
            if nuc.frames(j,it)
                dx=loc3(i).Centroid(1)-nuc.x(j,it);
                dy=loc3(i).Centroid(2)-nuc.y(j,it);
                d(i,j)=sqrt( dx^2 +dy^2 );
                if d(i,j)<nuc.radius(j,it)*fact_r  % If found a near nuclei
                    goodspot=true;
                end
            end
        end
        if goodspot
            inner_clusters=[inner_clusters i];
        end
    end
    %% eliminate the ones that are smaller than a certain number of voxels
    count_bigger_clusters=0;
    if numel(inner_clusters)
        clear bigger_clusters_tmp;
        bigger_clusters_tmp=find(([loc3(inner_clusters(:)).Area]>=voxels_min)&([loc3(inner_clusters(:)).Area]<=voxels_max));
        clear bigger_clusters;
        bigger_clusters=inner_clusters(bigger_clusters_tmp(:));
        if isempty(bigger_clusters)==0
            count_bigger_clusters=size(bigger_clusters,2);
        end
    end

    %%  For the survived clusters compute the total intensity
    if count_bigger_clusters>0  
        LL3D=zeros(size(II,1),size(II,2),z_max);
        Itot=zeros(size(bigger_clusters,2),1);
        %for i=1:size(inner_clusters,2)
        for i=1:size(bigger_clusters,2)
            ind=bigger_clusters(i);
            %ind=inner_clusters(i);
            ii=(L3D(:)==ind);
            LL3D(ii)=ind;
            loc3(ind).Area;
            %Itot(i)=sum(raw3(ii));
            Itot(i)=sum(filtered3(ii));
        end
        Itot(Itot<0)=0;
        Ispot=max(LL3D,[],3);
        %%  assign spots to nuclei
        %%%%%%%%%%%%%%%%%%###############################################

        %% Attribute nuclei to spots
        assigned_nucleus=zeros(size(nuc.frames,1),1);
        assigned_spot=zeros(size(bigger_clusters,2),1);
        count_spot=zeros(size(nuc.frames,1),1);
        for i=1:size(bigger_clusters,2)
            ind=bigger_clusters(i);
            [mindist,j]=min(d(ind,:));
            count_spot(j)=count_spot(j)+1;
            i_spot(j,count_spot(j))=i;
            dmin(bigger_clusters(i))=mindist;
        end
       %% Attribute spots to nuclei
        for j=1:size(nuc.frames,1)
            nu=nuc.ind(j,it);
            if count_spot(j)==1
                if nu>0
                    ind=bigger_clusters(i_spot(j,1));
                    assigned_nucleus(j)=ind;
                    assigned_spot(i_spot(j,1))=j;
                end
            end
            if count_spot(j)>1
                clear vtmp;
                vtmp=[Itot( i_spot(j,1:count_spot(j)) )];
                max_s=max(vtmp);
                max_is=find(vtmp==max_s(1));
                ispot=i_spot(j,max_is(1));
                ind=bigger_clusters(ispot);
                if nu>0 
                    assigned_nucleus(j)=ind;
                    assigned_spot(ispot)=j;
                end
            end
        end
        %%  fill the vector detected_spot 
        ind_spot=0;
        for i=1:size(bigger_clusters,2)
            ind=bigger_clusters(i);
            inuc=assigned_spot(i); %temporary index of nucleus
            if inuc>0
                nu=nuc.ind(inuc,it);
                if nu>0
                    %fprintf(1,'---- spot %d  belong to nucleus %d  at distance  %g\n',bigger_clusters(i),nu,d(bigger_clusters(i),inuc));
                    ind_spot=ind_spot+1;    
                    detected_spot.id_n(ind_spot)=nu;
                    detected_spot.id_s(ind_spot)=ind;
                    detected_spot.x(ind_spot)=loc3(ind).Centroid(1);
                    detected_spot.y(ind_spot)=loc3(ind).Centroid(2);
                    detected_spot.z(ind_spot)=loc3(ind).Centroid(3);
                    detected_spot.size(ind_spot)=loc3(ind).Area;
                    detected_spot.I(ind_spot)=Itot(i);   %the index i is not a mistake, it must be i and not ind
                    detected_spot.dist(ind_spot)=dmin(ind);
                    
                    % do the gaussian fit for each spot - not verified yet
                    [A,sx,sy,II2d,bk,res] = do_fit_gauss_2d_spot_projection(reader,channel,it,loc3(ind).Centroid(1),loc3(ind).Centroid(2),loc3(ind).Centroid(3),z_max,Lx,Ly,window);
                    detected_spot.A(ind_spot)=A;
                    detected_spot.ssx(ind_spot)=sx;
                    detected_spot.ssy(ind_spot)=sy;
                    detected_spot.I2d(ind_spot)=II2d;
                    detected_spot.bkg(ind_spot)=bk;
                    detected_spot.resid(ind_spot)=res;
                end
            end
        end        
        fprintf(1,'%d spot(s) found',sum(assigned_spot>0));
    end
end

