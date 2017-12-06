function [nuc,xfinrec,tfinrec,pfinrec]=reassign_cycle(nuc,dt,reso,drift_thresh,patch_crit,nc_range,nc_info,delimiter)
% Find t0 of the nuclei cycle:
    % nuc: nuclear information or the output file
    % dt: time resolution (in second)
    % res: xy pixel resolution (in micrometer)
    % drift_thresh: 1x2 vector containing threshold values for separation detection
        % drift_thresh(1): threshold for first large speed detection
        % drift_thresh(2): threshold for end of large speed detection
    % patch_crit: patch criteria
        % patch_crit(1): min time to patch before trace
        % patch_crit(2): min time to patch after trace
        % patch_crit(3): min trace length/interphase for correction
    % nc_range: range of nuclei cycle we want to reassign
    % nc_info (optional): Nx3 nuclear cycle info:
        % nc_info(:,1): nuclear cycle range
        % nc_info(:,2): begining of nuclear cycle
        % nc_info(:,3): end of nuclear cycle
    % delimiter (optional): delimiter in the output file
% Output:
    % cycle_nuc: [nuclei x time] specifying the nuclei information
    
    % Reconstruct nuc if the input is a file name:
        if ischar(nuc)
            filename=nuc;
            % Load the data file specified in filname
            nuc=struct;
            position_to_row=[];
            datamat=dlmread(filename,'\t',0,0);
            if size(datamat,2)<10
                datamat=dlmread(filename,' ',0,0);
            end
            dt=datamat(1,13);
            for i=1:size(datamat,1)
                % Set frames
                nuc.frames(datamat(i,1),datamat(i,2))=1;
                % Set ID
                nuc.ind(datamat(i,1),datamat(i,2))=datamat(i,1);
                % Set x
                nuc.x(datamat(i,1),datamat(i,2))=datamat(i,3);
                % Set y
                nuc.y(datamat(i,1),datamat(i,2))=datamat(i,4);
                % Set cell cycle
                nuc.cycle(datamat(i,1),datamat(i,2))=datamat(i,12);
                % Set mother
                nuc.mother(datamat(i,1),datamat(i,2))=datamat(i,9);
                % Set daughter1
                nuc.daughter1(datamat(i,1),datamat(i,2))=datamat(i,10);
                % Set daughter2
                nuc.daughter2(datamat(i,1),datamat(i,2))=datamat(i,11);
                % Set index from nuc.frames to row index in datamat
                position_to_row(datamat(i,1),datamat(i,2))=i;
            end
        else
            filename=[];
        end
        
        Nf=size(nuc.frames,2);
        % If nc_info is provided then reapply the nuclear cycle:
        if exist('nc_info','var')
            for i=1:size(nc_info,1)
                nuc.cycle(:,nc_info(i,2):nc_info(i,3))=nc_info(i,1);
            end
        else
            nc_info=[];
        end
        % There delimiter is specified:
        if ~exist('delimiter','var')
            delimiter='\t';
        end
        % Find nuclei cycle mean feature:
            function obj = fun1(x,fun2,isround)
                if isround
                    x(x~=round(x))=0;
                end
                obj=fun2(x(x>0));
                if isnan(obj)
                    obj=0;
                end
            end
    % Extract the features of nuclei
        x_nuc=arrayfun(@(i) fun1(nuc.x(i,:),@mean,0),1:size(nuc.x,1)); % Mean x position
        nuc_tmp=nuc.cycle.*(nuc.frames>0);
        nuc_tmp(nuc_tmp~=round(nuc_tmp))=0;
        cycle_nuc= arrayfun(@(i) fun1(nuc_tmp(i,:),@mode,1),1:size(nuc.x,1)); % General cell cycle
        daughter1_nuc = arrayfun(@(i) fun1(nuc.daughter1(i,:),@mode,0),1:size(nuc.x,1)); % Daughter 1 ID
        daughter2_nuc = arrayfun(@(i) fun1(nuc.daughter2(i,:),@mode,0),1:size(nuc.x,1)); % Daughter 2 ID
        daughter1_nuc(daughter1_nuc<=1)=0;
        daughter2_nuc(daughter2_nuc<=1)=0;
        life_nuc=arrayfun(@(i) fun1(nuc.frames(i,:)>0,@sum,0),1:size(nuc.x,1));  % Life time of nuclei        
        
    % Convert from nuclei ID to index:
        id_nuc=arrayfun(@(i) fun1(nuc.ind(i,:),@mode,0),1:size(nuc.x,1)); % ID
        for i=1:size(nuc.x)
            if id_nuc(i)
                idx_nuc(id_nuc(i))=i;
            end
        end
    % Find the separation point:    
        % Find nuclei with two daughter2
            valid_nuc=find(x_nuc&cycle_nuc&daughter1_nuc&daughter2_nuc);
        % Find every pair of mother - daughters:
            max_frame=round(120/dt);             % Track cells 30s after mitosis
            sep_rec=zeros(1,numel(valid_nuc));  % Record separation time
            div_rec=zeros(1,numel(valid_nuc));  % Record division time
            nc_rec=zeros(1,numel(valid_nuc));   % Record of nuclei cycle
            x_rec=zeros(1,numel(valid_nuc));   % Record of nuclei position
            dis_rec=zeros(numel(valid_nuc),max_frame);
            
        % Find the distance and drift-away speedbetween the daughter cells 
        % following separation (detected by segmentation tool)
            cnt=0;
            for i=valid_nuc
                % If a cell has two daughters in the next cell cycle
                idx_d1=idx_nuc(daughter1_nuc(i));
                idx_d2=idx_nuc(daughter2_nuc(i));
                if cycle_nuc(idx_d1)==cycle_nuc(idx_d2)
                    if cycle_nuc(idx_d1)==cycle_nuc(i)+1
                        cnt=cnt+1;
                        nc_rec(cnt)=cycle_nuc(i);
                        x_rec(cnt)=x_nuc(i);
                        % Find the time/position when they divide:
                        div_time=find(nuc.frames(i,:),1,'last');
                        % Find the time/position when they are farthest
                        for ttime=1:max_frame
                            if div_time+ttime<=Nf
                                if nuc.frames(idx_d1,div_time+ttime)&nuc.frames(idx_d2,div_time+ttime)
                                    dis_rec(cnt,ttime)=sqrt( ...
                                        (nuc.x(idx_d1,div_time+ttime) ...
                                        -nuc.x(idx_d2,div_time+ttime))^2 +...
                                        (nuc.y(idx_d1,div_time+ttime) ...
                                        -nuc.y(idx_d2,div_time+ttime))^2);
                                else
                                    break;
                                end
                            end
                        end
                        div_rec(cnt)=div_time;
                    end
                end
            end
            
            % Plot separation distance and speed
            figure;
            subplot(131);
            plot(dis_rec'*reso);    % Plot separation speed
            xlabel('Frame from separation');
            ylabel('Nuclei distance (\mum)');
            subplot(132);
            sep_speed=diff(dis_rec'*reso/dt);
            plot(sep_speed);    % Plot separation speed
            xlabel('Frame from separation');
            ylabel('Drift speed (\mum/s)');
            % Find sepration time
            first_sep=arrayfun(@(x) find(sep_speed(:,x)>drift_thresh(1),1,'first'),1:cnt,'UniformOutput',false);
            first_sep(cellfun(@isempty,first_sep))={0};
            last_sep=arrayfun(@(x) find(sep_speed(first_sep{x}+1:end,x)<drift_thresh(2),1,'first'),1:cnt,'UniformOutput',false);
            last_sep(cellfun(@isempty,last_sep))={0};
            for i=1:cnt
                if ~(first_sep{i})
                    sep_rec(i)=div_rec(i);
                else
                    if last_sep{i}==0
                        sep_rec(i)=NaN;
                    else
                        sep_rec(i)=div_rec(i)+first_sep{i}-1+last_sep{i};
                    end
                end
            end
            subplot(133);       % Plot separation time
            plot(x_rec,sep_rec,'xb','DisplayName','t0');hold on;
            plot(x_rec,div_rec,'xr','DisplayName','Sep time');
            
     % Fit and find frame for t0 at each given position:
        % Set the position axis
        ax=sort(x_rec);
        ax=ax(ax>0);
        % Set recorders
        pfinrec=cell(1,15);               % For polynomial fit
        end_cycle=zeros(1,15);      % For cycle starting time
        xfinrec=cell(1,15);
        tfinrec=cell(1,15);
        for cycle=1:15
            idselect1=(nc_rec==cycle)&(~isnan(sep_rec));
            idselect=(abs(sep_rec-mode(sep_rec(idselect1)))*dt<100)&idselect1;
            if cycle&(sum(idselect)>=8)
                % Fit the sep_time of same cycle
                xfinrec{cycle}=x_rec(idselect);
                tfinrec{cycle}=sep_rec(idselect);
                % Convert frame, pixel to second, %EL
                tfinrec{cycle}=tfinrec{cycle}*dt;
                xfinrec{cycle}=(datamat(1,16)+xfinrec{cycle})/(datamat(1,16)+datamat(1,17)+datamat(1,15));
                xfeed=xfinrec{cycle};
                yfeed=tfinrec{cycle};
                pfinrec{cycle} = fminsearchbnd(@ffit,[mean(yfeed) mean(xfeed) 1],[0 0.2 0],[1e5 0.8 1e5]);
                tmp=mitotic_val(pfinrec{cycle},[0:0.01:1]);
                plot([0:0.01:1]*(datamat(1,16)+datamat(1,17)+datamat(1,15))-datamat(1,16),tmp/dt,'DisplayName',['nc' num2str(cycle) '-nc' num2str(cycle+1)],'LineWidth',2,'LineStyle','--');
                end_cycle(cycle)=mean(tmp);                
            end
        end
        legend show;
    
    % Recalibrating the trace length and export the file    
    if numel(filename)
        k = strfind(filename,'.');
        outfile=[filename(1:k(end)-1) '_fixed.txt'];
        
        datamat_=datamat;

        % Begin realigning the traces
        for cycle=nc_range
             % Find the interphase duration:
            if ~end_cycle(cycle-1)
                end_cycle(cycle-1)=1;
            end
            if end_cycle(cycle)
                interphase_duration=end_cycle(cycle)-end_cycle(cycle-1);
            else
                interphase_duration=size(nuc.frames,2)-end_cycle(cycle-1);
            end
            % Select traces in the cycle that is neither too long or short
            idselect=find((cycle_nuc==cycle)&...
                (life_nuc>=interphase_duration*patch_crit(3)/100)...
                &((life_nuc<interphase_duration*1.2)));
            % Begin adjusting the nuclei traces:
            for i=idselect
                % Find predicted begin and end of traces:
                    % Begin
                        if pfinrec{cycle-1}
                            start_it=mitotic_val(pfinrec{cycle-1},x_nuc(i));
                        else
                            start_it=1e5;
                        end
                    % End
                        if pfinrec{cycle}
                            stop_it=mitotic_val(pfinrec{cycle},x_nuc(i));
                        else
                            stop_it=0;
                        end
                    % If traces end too early:
                    t_bg=find(nuc.frames(i,:),1,'first');
                    t_end=find(nuc.frames(i,:),1,'last');
                    if ((t_bg-start_it)*dt>patch_crit(1))||((stop_it-t_end)*dt>patch_crit(2))
                        continue;
                    end
                % Creates a blank row
                    empty_row=datamat(1,:)*0;
                    % Create blank row to account for time before real trace
                        empty_row(1:19)=datamat(position_to_row(i,t_bg),1:19);
                        % Put new column at the end of datamat matrix:
                        for j=ceil(start_it):t_bg-1
                            empty_row(2)=j;
                            empty_row(12)=cycle;
                            datamat_(end+1,:)=empty_row;
                            nuc.frames(i,j)=1;
                        end                        
                    % Create blank row to account for time after real trace
                        empty_row(1:19)=datamat(position_to_row(i,t_end),1:19);
                        % Put new column at the end of datamat matrix:
                        for j=t_end+1:floor(stop_it)
                            empty_row(2)=j;
                            empty_row(12)=cycle;
                            datamat_(end+1,:)=empty_row;
                            nuc.frames(i,j)=1;
                        end
                    % Set the remaining cycle to the right value:
                        for j=t_bg:t_end
                            datamat_(position_to_row(i,j),12)=cycle;                            
                        end
            end
        end
        
        % Export the output file if the input file is a string:
        dlmwrite(outfile,datamat_,'precision',6,'Delimiter',delimiter);
        
        % Compare nuclei life time before and after:
        life_nuc_=arrayfun(@(i) fun1(nuc.frames(i,:)>0,@sum,0),1:size(nuc.x,1));  % Life time of nuclei after calibrating
        figure;
        plot(life_nuc_*dt);hold on;
        plot(life_nuc*dt);
        ylabel('Nuclei life time (s)');
        xlabel('Nuclei');
    end

    function obj=ffit(p)
        % p: function parameters
            % p(1): max position for y
            % p(2): x that maximize y
            % p(3): speed that y is decreasing
        obj=sum((yfeed-mitotic_val(p,xfeed)).^2);
    end
end
