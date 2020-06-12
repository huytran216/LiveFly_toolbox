filename = 'H6B6_170525_CP1';
matfile=dlmread([filename '.txt'],',');
[x_resolution,z_resolution] = get_resolution_config(filename);
plot_report=1;
%% Plot 
nc=13;
idselect= find((matfile(:,12)==nc)&(matfile(:,20)>0)&(matfile(:,32)>0));
figure(1);
subplot(121);
plot(matfile(idselect,21),matfile(idselect,32)-matfile(idselect,22),'o');
xlabel('y position');
ylabel('Distance in z plane');
title('Unoriented spot position');
subplot(122);
plot(matfile(idselect,2),matfile(idselect,32)-matfile(idselect,22),'o');
xlabel('Frame');
ylabel('Distance in z plane');
title('Unoriented spot position');
mtmp=[];
stmp=[];
unique_pos=unique(matfile(idselect,2));
for i=1:numel(unique_pos)
    idtmp = find(matfile(idselect,2)==unique_pos(i));
    mtmp(i)= mean(matfile(idselect(idtmp),32)-matfile(idselect(idtmp),22));
    stmp(i)=sqrt(var(matfile(idselect(idtmp),32)-matfile(idselect(idtmp),22)));
end
hold on;
errorbar(unique_pos,mtmp,stmp);
%% Find the embryo surface
n_frame = unique(matfile(:,2));
%n_frame =100;
figure;
for frame=n_frame(:)'
    frame
    matid_range = find((matfile(:,2)==frame)&(matfile(:,32)>0));
    if any(matfile(matid_range,20)>0)
        x=matfile(matid_range,3)*x_resolution;
        y=matfile(matid_range,4)*x_resolution;
        z=matfile(matid_range,32)*z_resolution;
        % Fit embryo membrane
        f = fit( [x, y], z, 'poly44');
        figure(2);
        ax=linspace(min(x),max(x),20);
        ay=linspace(min(y),max(y),20);
        [tmpax,tmpay]=meshgrid(ax,ay);
        az=arrayfun(@(x,y) f(x,y),tmpax,tmpay);
        plot3(x,y,z,'o');hold on;
        mesh(ax,ay,az);
        xlabel('x');
        ylabel('y');
        title('Nuclei position');
        % For each nuclei, find the orientation and find the spot:
        for i=1:numel(matid_range)
            id=matid_range(i);
            % Find closest distance to surface:
            x0=matfile(id,3)*x_resolution;
            y0=matfile(id,4)*x_resolution;
            z0=matfile(id,32)*z_resolution;
            fun =@(cin) (cin(1)-x0).^2 + (cin(2)-y0).^2 + (f(cin(1),cin(2))-z0).^2;
            beta = fminsearch(fun,[x0+1e-1 y0+1e-1]);
            x1 = beta(1);
            y1 = beta(2);
            z1 = f(x1,y1);
            % Find the orientation:
            if (x1-x0~=0)
                if (matfile(id,20))
                    ez = [x1-x0 y1-y0 z1-z0];
                    if f(x0,y0)<z0
                        ez = -ez;
                    end
                    ey = [0 -(z1-z0) y1-y0];
                    ex = [(-(y1-y0)^2-(z1-z0)^2)/(x1-x0)  y1-y0 z1-z0];
                end
            else
                'ops'
            end
            % Find the spot location:
            if matfile(id,20)
                ez = ez./sqrt(sum(ez.^2));
                ex = ex./sqrt(sum(ex.^2));
                ey = ey./sqrt(sum(ey.^2));
                E = [ex;ey;ez];
                cout = (-matfile(id,[20 21 22]).*[x_resolution x_resolution z_resolution]+...
                    matfile(id,[3 4 32]).*[x_resolution x_resolution z_resolution])*E;
                [cout;(-matfile(id,[20 21 22]).*[x_resolution x_resolution z_resolution]+...
                    matfile(id,[3 4 32]).*[x_resolution x_resolution z_resolution])]
                matfile(id,[33 34 35])=cout;
                % Plot the results
                figure(3);
                plot3(cout(1),cout(2),cout(3),'o'); hold on;
                title('Calibrated spot position');
            end
        end
    end
end
dlmwrite([filename '_.txt'],matfile);