function [Irec,AdjustedIntensity,trec,Adjustedtime,fearec,x,y,xrec,yrec,sizerec,xspotrec,yspotrec,id,mother,daughter1,daughter2,cyclerec,dt,Imax,BG,ShiftL,ShiftR,APole]=extract_movie(fullpath,tinterphase,BG_man)
    % Extract the features (1:8) from the time series
    % Input:
    %    % datamat: input matrix or path to data file
    %    % isfix: Is a fix file required or not (to remove unwanted spot detection)
    % Output:
    %    % Irec: record of cell intensity over time
    %    % fearec: list of feature extracted
    %    % xrec: position along the x-axis
    %    % yrec: position along the y-axis
    %    % id: cell idAdjustedIntensity
    %    % mother: mother id
    %    % daughter1: daughter 1 id
    %    % daughter2: daughter 2 id
    %    % cycle: cell cycle
    %    % dt: time resolution (s)
    %    % Imax: Embryo maximum intensity 
    %% Junk correction file:
    % Actually not needed, but needed anyway for the stupid static
    % workspace problem
    delimiter='a';
    xlen=0;
    ylen=0;
    xlim_left=0;
    xlim_right=0;
    ylim_up=0;
    ylim_down=0;
    Icolumn=28;
    rm={};
    nclist={};
    pixel_ignore_posterior = 1e10;
    
    range1=1:10001;filename1='Result_file1';
    range2=10002:30000;filename2='Result_file2';
    border13=10000;
    keepid13=[];
    leftcensored=0;
    rightcensored=0;
    startcutcycle=[0 0 0 0 0 0];
    endcutcycle=[10000 10000 10000 10000 10000 10000];
%% Load the fix file
    run(fullfile(fullpath,'correction.m'));
    if ~exist('rm','var')
        rm={};
    end
    if ~exist('nclist','var')
        nclist={};
    end
    if ~exist('border13','var')
        border13=10e5;
    end
    if ~exist('keepid13','var')
        keepid13=[];
    end
    if (~exist('tinterphase'))||(tinterphase==0)
        tinterphase=ones(1,6);
    end
    if ~exist('BG_man')
        BG_man=1;
    end
    % Remove cell within 1%EL near frame anyway
    xlim_left=25;
    xlim_right=25;
    ylim_up=25;
    ylim_down=25;
%% Load the data file
    datamat1=[];
    delimiter_list={delimiter,' ','\t',','};
    cnt=0;
    while size(datamat1,2)<=10
        cnt=cnt+1;
        if cnt<=numel(delimiter_list)
            delimiter=delimiter_list{cnt};
        else
            msgbox(['Corupted or non-extant ''' fullpath '''']);
            break;
        end
        try
            datamat1=dlmread(fullfile(fullpath,[filename1]),delimiter,0,0);
        catch
            datamat1=[];
        end
    end
    datamat1=datamat1(ismember(datamat1(:,2),range1),:);
    if ~strcmp(filename1,filename2)
        datamat2=dlmread(fullfile(fullpath,[filename2]),delimiter,0,0);
        datamat2=datamat2(ismember(datamat2(:,2),range2),:);
        datamat=[datamat1;datamat2];
    else
        datamat=datamat1;
    end
    clear datamat1 datamat2;
    % Find the time resolution
    dt=datamat(1,13);
    % Find ShiftL, ShiftR, APole
    ShiftL=datamat(1,16);
    ShiftR=datamat(1,17);
    APole=datamat(1,19);
    % Extract the background intensity from gaussian fit
    BG=datamat(:,29);
    BG=BG(BG>0);
    BGsigma = sqrt(var(BG));
    BG = mean(BG);
    % Remove spots with weak intensity
    datamat(datamat(:,25)<=2*BGsigma,[20:30])=0;
%% Flip the x value if A pole is 1
    if (xlen~=datamat(1,14))&&(ylen~=datamat(1,15))
        display(['Warning about Lx, Ly in movie ' fullpath]);
    end
    if (xlen==0)||(ylen==0)
        msgbox(['Check xlen ylen in ' fullpath]);
    end
    xpos_range=[3 20];
    if ~datamat(1,19)
        for xpos=xpos_range
            datamat(:,xpos)=xlen-datamat(:,xpos)+1;
        end
        % Flip Shift_L and Shift_R
        tmp=datamat(:,16);
        datamat(:,16)=datamat(:,17);
        datamat(:,17)=tmp;
    end
%% Ignore spots beyond pixel_ignore_posterior
    if ~datamat(1,19)
        flttmp = datamat(:,xpos_range(1))<pixel_ignore_posterior;
    else
        flttmp = datamat(:,xpos_range(1))>=pixel_ignore_posterior;
    end
    datamat(flttmp,Icolumn)=0;
%% Reconstruct cell lineage - See Header.txt for reference
    id=datamat(:,1);
    mother=datamat(:,9);
    daughter1=datamat(:,10);
    daughter2=datamat(:,11);
    [id,uniqueid]=unique(id);
    missingid=find(arrayfun(@(x) ~any(id==x),1:max(id)));
    
    mother=mother(uniqueid);mother(id)=mother;mother(missingid)=0;
    daughter1=daughter1(uniqueid);daughter1(id)=daughter1;daughter1(missingid)=0;
    daughter2=daughter2(uniqueid);daughter2(id)=daughter2;daughter2(missingid)=0;
%% Gather the cell intensity traces, position and lineage information
    Irec=cell(max(id),1);       % Record spot intensity
    trec=cell(max(id),1);       % Record frame
    cyclerec=cell(max(id),1);   % Record cell cycle
    fearec=cell(max(id),1);     % Record features
    xrec=cell(max(id),1);       % Record nuclei position
    yrec=cell(max(id),1);       % Record nuclei position
    xspotrec=cell(max(id),1);   % Record spot position
    yspotrec=cell(max(id),1);   % Record spot position
    sizerec=cell(max(id),1);    % Record nuclei size

    Imax=[];                    % Record maximum intensity
    tmax=max([datamat(:,2)])-1;
    tmin=dt;
    
    cycletmp=cell(1,14);                % Record interphase durations
    for i=1:size(datamat,1)
        trec{datamat(i,1)}=[trec{datamat(i,1)} datamat(i,2)];
        cyclerec{datamat(i,1)}=[cyclerec{datamat(i,1)} datamat(i,12)];
        
        Ival=datamat(i,Icolumn)/BG_man;
        % Peform correcting the data based on correction.m file
        if numel(rm)>=i
            if ismember(datamat(i,2),rm{i})
                Ival=0;
            end
        end
        if numel(nclist)
            if ismember(datamat(i,1),rm{i})
                Ival=0;
            end
        end
        if (datamat(i,3)>=border13)&&(~ismember(i,keepid13))&&(datamat(i,12)==13)
            Ival=0;
        end
        Irec{datamat(i,1)}=[Irec{datamat(i,1)} Ival];
        % Calculate the cell position (%) along the AP axis (x)
        xrec{datamat(i,1)}=[xrec{datamat(i,1)} datamat(i,3)];
        yrec{datamat(i,1)}=[yrec{datamat(i,1)} datamat(i,4)];
        x(datamat(i,1))=mean((datamat(i,3)+datamat(i,16))./(xlen+datamat(i,16)+datamat(i,17)));
        y(datamat(i,1))=mean(datamat(i,4)./((ylen+datamat(i,16)+datamat(i,17))));
        % Calculate the relative spot position from the nuclei center
        xspotrec{datamat(i,1)}=[xspotrec{datamat(i,1)} datamat(i,3)-datamat(i,20)];
        yspotrec{datamat(i,1)}=[yspotrec{datamat(i,1)} datamat(i,4)-datamat(i,21)];
        % Calculate the spot size
        sizerec{datamat(i,1)}=[sizerec{datamat(i,1)} datamat(i,5)];
        % Check if the nuclei is located in the edge of the images
        if (mean(xrec{datamat(i,1)})<xlim_left) ...
                &&(mean(xrec{datamat(i,1)})>xlen-xlim_right) ...
                &&(mean(yrec{datamat(i,1)})>ylen-ylim_down) ...
                &&(mean(yrec{datamat(i,1)})<ylim_up)
            censored(datamat(i,1))=-1;
        else
            censored(datamat(i,1))=0;
        end
        if numel(max(Irec{datamat(i,1)}))
            Imax(datamat(i,1))=max(Irec{datamat(i,1)});
        else
            Imax(datamat(i,1))=0;
        end
    end
%% Post process and extract the features
    iodd=[]; % Storing odd traces (to be removed?)
    for i=1:max(id)
        [~,tmp]=sort(trec{i});
        Irec{i}=Irec{i}(tmp);
        trec{i}=trec{i}(tmp);
        cyclerec{i}=cyclerec{i}(tmp);
        % Select only interphase period to analysize
        tmp= (cyclerec{i}-floor(cyclerec{i}))==0;
        if (sum(tmp)>0)
            cyclerec{i}=mode(cyclerec{i}(tmp));
            Irec{i}=Irec{i}(tmp);
            trec{i}=trec{i}(tmp);
            xrec{i}=xrec{i}(tmp);
            yrec{i}=yrec{i}(tmp);
            xspotrec{i}=xspotrec{i}(tmp);
            yspotrec{i}=yspotrec{i}(tmp);
            sizerec{i}=sizerec{i}(tmp);
        else
            cyclerec{i}=9;  % unclassified nuclei
            Irec{i}=0;
            trec{i}=0;
            xrec{i}=0;
            yrec{i}=0;
            xspotrec{i}=0;
            yspotrec{i}=0;
            sizerec{i}=0;
        end
        if ~(cyclerec{i}>9)
            cyclerec{i}=9;
        end
        % Remove spots that survive mitosis
        Irec{i}(trec{i}<startcutcycle(cyclerec{i}-8))=0;
        Irec{i}(trec{i}>endcutcycle(cyclerec{i}-8))=0;
        % Calculate the interphase duration
        cycletmp{cyclerec{i}}=[cycletmp{cyclerec{i}} numel(trec{i})];
        
        Imax(i)=-1;
        if (x(i)<0.4)&&numel(Irec{i})
            Imax(i)=max(Irec{i});
        end        
    end
    
    [ncnt,ax]=hist(Imax(Imax>100));
    [~,pos]=max(ncnt);
    Imax=ax(pos);
    
    for i=1:max(id)
        % Remove irregular cells or cells in the end of time series
        if numel(trec{i})
            if (trec{i}(end)>=tmax)
                if rightcensored
                    cencored(i)=-1;     % Invalid cell,
                else
                    censored(i)=1;      % Right censored cell
                end
            end
            if (trec{i}(1)<=tmin)
                if leftcensored
                    censored(i)=-1;         % invalid cell
                else
                    censored(i)=2;          % Left censored cell
                end
            end
            if numel(trec{i})<median(cycletmp{cyclerec{i}})/2
                censored(i)=-1;    % Invalid cell - too short interphase
            end
        else
            censored(i)=-1;        % Cell does not exist
        end
        % Extract the features - no trimming this time
        threshold=0;
        [fearec{i},AdjustedIntensity{i},Adjustedtime{i}]=extract_feature(Irec{i},trec{i},threshold,dt,Imax,censored(i));
%         if cyclerec{i}==12
%             'pause'
%             fearec{i}
%             length(Irec{i})*dt
%         end
        
        trec{i}=trec{i}*dt;
%         % ALERT IF WEIRD FEATURES ARE DECTED
%             % Early activating cell (most likely due to aggregrate)
%             if (fearec{i}(2)<0.1)&&(fearec{i}(2)>0)&&(cyclerec{i}>9)
%                 %display(['Weird activation time (' num2str(fearec{i}(2)) '), id: ' num2str(i) ', x: ' num2str(x(i)) ', y: ' num2str(y(i))]);
%                 %display(['cycle ' num2str(cyclerec{i}),', Path: ' fullpath]);
%                 plot(trec{i},Irec{i});hold on;
%                 iodd=[iodd i];
%             end
%             title('Weird traces detected');
%             xlabel('time (s)');
%             % Cell with no spot at the anterior
%              if (fearec{i}(1)==0)&&(x(i)<0.45)
%                  fearec{i}=extract_feature;
%              end
            if numel(fearec{i})==0
                woap
            end
    end
    cyclerec(missingid)={9};
    id=1:max(id);
    iodd
    cyclerec(iodd)
end 