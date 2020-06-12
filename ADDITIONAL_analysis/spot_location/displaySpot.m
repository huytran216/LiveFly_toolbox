%% Choose movie
movie = 'H6B6C6_';
movieidx = 1;
nc=13;
%% Load data from movie file
matfile=[];
listing = dir(pwd);
cnt=0;
for i = 1:numel(listing)
    if strfind(listing(i).name,'WT_')&strfind(listing(i).name,'_.txt')
        cnt=cnt+1;
        if (movieidx==0)||(cnt==movieidx)
            tmp=dlmread(listing(i).name,',');
            [x_resolution,z_resolution]=get_resolution_config(listing(i).name);
            tmp(:,[3 4 5 14 15 16 17 20 21 26 27 33 34])=tmp(:,[3 4 5 14 15 16 17 20 21 26 27 33 34])*x_resolution;
            tmp(:,[22 32 35])=tmp(:,[22 32 35])*z_resolution;
            matfile=[matfile;tmp];
        end
    end
end
%% Find the embryo surface
n_frame = unique(matfile(:,2));
%n_frame =100;

figure;
for frame=n_frame(:)'
    frame
    matid_range = find((matfile(:,2)==frame)&(matfile(:,32)>0)&(matfile(:,20)>0)&(matfile(:,12)==nc));
    if any(matfile(matid_range,20)>0)
        if numel(n_frame)==1
            % Get average nuclei size
            matnc_range = find((matfile(:,2)==frame)&(matfile(:,32)>0));
            fsize = median(matfile(matnc_range,5));
            % Draw spots
            plot3(matfile(matid_range,33),matfile(matid_range,34),matfile(matid_range,35),'o'); hold on;
            % Draw circles:
            plot3([0 0],[0 0],[-fsize fsize],'LineWidth',2);
            plot3([0 0],[-fsize fsize],[0 0],'LineWidth',1);
            plot3([-fsize fsize],[0 0],[0 0],'LineWidth',1);
            a=0:0.1:2*pi;
            for ztmp = linspace(-fsize,fsize,11)
                fn = sqrt(fsize^2 - ztmp^2);
                plot3(fn*cos(a),fn*sin(a),a*0+ztmp);
            end
            axis square;
        else
            % Normalize
            % Get average nuclei size
            matnc_range = find((matfile(:,2)==frame)&(matfile(:,32)>0));
            fsize = median(matfile(matnc_range,5));
            % Draw spots
            plot3(matfile(matid_range,33)/fsize,matfile(matid_range,34)/fsize,matfile(matid_range,35)/fsize,'o'); hold on;

            % Draw circles:
            plot3([0 0],[0 0],[-1 1],'LineWidth',2);
            plot3([0 0],[-1 1],[0 0],'LineWidth',1);
            plot3([-1 1],[0 0],[0 0],'LineWidth',1);
            a=0:0.1:2*pi;
            for ztmp = linspace(-0.95,0.95,11)
                fn = sqrt(1^2 - ztmp^2);
                plot3(fn*cos(a),fn*sin(a),a*0+ztmp,'color','k');
            end
            xlim([-2 2]);
            ylim([-2 2]);
            zlim([-2 2]);
        end
    end
end
%%
idselect= find((matfile(:,12)==nc)&(matfile(:,20)>0)&(matfile(:,32)>0));
figure;
% z frame vs y position
    subplot(231);
    plot(matfile(idselect,4),matfile(idselect,35),'o');
    xlabel('y position (micron)');
    ylabel('Distance in z plane (micron)');
    mtmp=[];
    stmp=[];
    [Y,E] = discretize(matfile(idselect,4),20);
    for i=1:numel(E)
        idtmp = find(Y==i);
        mtmp(i)= mean(matfile(idselect(idtmp),35));
        stmp(i)=sqrt(var(matfile(idselect(idtmp),35)));
    end
    hold on;
    errorbar(E,mtmp,stmp);
% z frame vs x position
    subplot(234);
    plot(matfile(idselect,3),matfile(idselect,35),'o');
    xlabel('x position (micron)');
    ylabel('Distance in z plane (micron)');
    mtmp=[];
    stmp=[];
    [Y,E] = discretize(matfile(idselect,3),20);
    for i=1:numel(E)
        idtmp = find(Y==i);
        mtmp(i)= mean(matfile(idselect(idtmp),35));
        stmp(i)=sqrt(var(matfile(idselect(idtmp),35)));
    end
    hold on;
    errorbar(E,mtmp,stmp);
% z frame vs frame
    subplot(232);
    plot(matfile(idselect,2),matfile(idselect,35),'o');
    xlabel('Frame');
    ylabel('Distance in z plane (micron)');
    mtmp=[];
    stmp=[];
    unique_pos=unique(matfile(idselect,2));
    for i=1:numel(unique_pos)
        idtmp = find(matfile(idselect,2)==unique_pos(i));
        mtmp(i)= mean(matfile(idselect(idtmp),35));
        stmp(i)=sqrt(var(matfile(idselect(idtmp),35)));
    end
    hold on;
    errorbar(unique_pos,mtmp,stmp);
% z frame vs intensity
    subplot(233);
    plot(matfile(idselect,28),matfile(idselect,35),'o');
    xlabel('Spot intensity');
    ylabel('Distance in z plane (micron)');
    mtmp=[];
    stmp=[];
    [Y,E] = discretize(matfile(idselect,28),20);
    for i=1:numel(E)
        idtmp = find(Y==i);
        mtmp(i)= mean(matfile(idselect(idtmp),35));
        stmp(i)=sqrt(var(matfile(idselect(idtmp),35)));
    end
    hold on;
    errorbar(E,mtmp,stmp);
% distance to center vs intensity
    subplot(235);
    plot(matfile(idselect,28),sqrt(matfile(idselect,33).^2+matfile(idselect,34).^2),'o');
    xlabel('Spot intensity');
    ylabel('Distance from z-axis (micron)');
    mtmp=[];
    stmp=[];
    [Y,E] = discretize(matfile(idselect,28),20);
    for i=1:numel(E)
        idtmp = find(Y==i);
        mtmp(i)= mean(sqrt(matfile(idselect(idtmp),33).^2+matfile(idselect(idtmp),34).^2));
        stmp(i)=sqrt(var(sqrt(matfile(idselect(idtmp),33).^2+matfile(idselect(idtmp),34).^2)));
    end
    hold on;
    errorbar(E,mtmp,stmp);
% Distance from sphere center
    subplot(236);
    [tmp,ax]=hist(sqrt(matfile(idselect,33).^2+matfile(idselect,34).^2 + matfile(idselect,35).^2),20);
    plot(ax,tmp/sum(tmp));
    % Get uniform distribution
    