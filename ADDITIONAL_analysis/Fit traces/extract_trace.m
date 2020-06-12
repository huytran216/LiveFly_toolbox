%% Re-extract all features in all strains before running the script
% The data should be in Data_analysis/
%% Params
addpath '..\..\Tool\hmm_dont_edit';
addpath '..\..\Tool\hmm_dont_edit\utilities\';
fld='../../Data analysis/datasets'; % Location of the dataset

%% Select data 
Dataset='Hb-M1-New';                   % Name of the dataset
cycleno = [13];                     % Interphase duration

% Steady state windows
if cycleno==12
    time_SS = [450 550];            % Begin _ end time (s)
end
if cycleno==13
    time_SS = [450 800];            % Begin _ end time (s)
end

%time_SS = [0 800];

% Extract position
position_SS = [-35 -20];            % Begin _ end position (%EL) for B6
%position_SS = [-35 -20 ];            % Begin _ end position (%EL) for WT-M1-New
%position_SS = [-35 -30];            % Begin _ end position (%EL) for B6
%% MS2 configuration and memory length
w=5;

Nbs=[24];         % ms2 Configuration 1 (number of MS2 stemloops)
Padsize=[600];    % ms2 Configuration 2 (length of random sequence after loops in bp)

Lgen = Nbs*52+Padsize;

dt= Lgen/w/45;
wL=get_wL(Nbs,Padsize,dt);
scale_time = 1;                     % Scaling interphase by mean interphase duration

filename = ['trace_' Dataset '_nc' num2str(cycleno) ];
%% Load the dataset
load(fullfile(fld,Dataset),'Nmov','datamat','DatasetFeature','DatasetList');
%% Define
subindex = @(x,y) x(y);
%% Color order:
 corder=    [0    0.4470    0.7410
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560;...
    0.4660    0.6740    0.1880;...
    0.3010    0.7450    0.9330;...
    0.6350    0.0780    0.1840];
%% Refine interphase
ts_spec=find(arrayfun(@(x) getfield(DatasetList,{x},['nc' num2str(cycleno)]),1:Nmov));
trace_rec={};
time_all={};
trace_all={};
figure;
if numel(ts_spec)
    % Get interphase duration
    tphase=DatasetFeature(cycleno-8).vborder_rec(10,:)*2;
    if scale_time
        tphase_set = ones(size(tphase))*mean(tphase(tphase>0));
    else
        tphase_set = tphase;
    end
    total=0;
    tmax=0;
    cntempty=0; 
    minpositive=100000;
    cnt=0;
    for tsidx = ts_spec
        cnt=cnt+1;
        % Get nuclei id that match cycleno and time series and position
        idselect=find(([datamat(:).cycle]==cycleno)&([datamat(:).tscnt]==tsidx)&(100*([datamat(:).x]-0.5)>=position_SS(1))&(100*([datamat(:).x]-0.5)<=position_SS(2)));
        ison=arrayfun(@(x) subindex(datamat(x).Feature,1),idselect); % Is cell valid
        idselect = idselect(ison>=0);
        time_ax = time_SS(1):dt:time_SS(2);
        cnt=0;
        for id=idselect
            total=total+1;
            cnt=cnt+1;
            % Get traces
            tr = interp1([datamat(id).time-datamat(id).time(1) 1e5],[datamat(id).Intensity -1e10],time_ax*tphase(tsidx)/tphase_set(tsidx));
            trace_rec{tsidx,cnt}=tr(tr>=0);
            trace_all{total}=tr(tr>=0);
            time_all{total}=sum(tr>=0);
            tmax=max(time_all{total},tmax);
            if sum(tr>0)==0
                cntempty=cntempty+1;
            end
            if minpositive>min(tr(tr>0))
                if min(tr(tr>0))>0.01
                    minpositive = min(tr(tr>0));
                end
            end
            %if (total>10)&(total<15)
                plot(time_ax(tr>=0),tr(tr>=0),'color',corder(mod(cnt,7)+1,:));hold on;
            %end
        end
    end
    % Plot mean curve
    Irec=cell(1,tmax);
    for i=1:total
        for j=1:time_all{i}
            Irec{j}=[Irec{j} trace_all{i}(j)];
        end
    end
    mIrec=[];
    sIrec=[];
    for j=1:numel(Irec)
        mIrec(j)=mean(Irec{j});
        sIrec(j)=sqrt(var(Irec{j}));
    end
end
errorbar(time_ax,mIrec,sIrec,'Display','Dataset','LineWidth',2,'color','k');
display([num2str(total) ' traces extracted']);
display([num2str(cntempty),' empty traces']);
%% Save the data
mkdir('trace_data');
save(['trace_data/' filename],'trace_rec','trace_all','dt','total','time_all','mIrec','sIrec','time_SS','position_SS','w','wL','minpositive');
