function [datamat,DatasetList]=local_extract(DatasetList)
    f=waitbar(0,'Loading');
    
    % Extracting data from the movies and combine them into one
    datamat=struct('id',[],'Intensity',{},'Feature',[],'x',[],'y',[],'border',[],'mother',[],'daughter1',[],'daughter2',[],'tscnt',[],'cycle',[],'tlen',[],'start',[]);
    
    % Calculate shifting level per embryo for zero mean shift
    shift_EL=zeros(1,numel(DatasetList));
    for idx=1:numel(DatasetList)
        shift_EL(idx) = DatasetList(idx).shift_EL;
    end
    shift_EL=shift_EL-mean(shift_EL);
    % Load the movies
    cnt=0;
    for idx=1:numel(DatasetList)
        waitbar((idx-1)/numel(DatasetList),f);
        fullpath=DatasetList(idx).Path;
        [Itmp,Itmp2,trec,trec2,fea,x,y,z,xrec,yrec,zrec,sizerec,xspotrec,yspotrec,id,mother,daughter1,daughter2,cycle,dt,Imax,BG,BGsigma,ShiftL,ShiftR,APole]=extract_movie(fullpath,0,DatasetList(idx).BG_man);
        DatasetList(idx).shift_EL=shift_EL(idx);
        DatasetList(idx).BG=BG;
        DatasetList(idx).BGsigma=BGsigma;
        DatasetList(idx).dt=dt;
        DatasetList(idx).APole=APole;
        DatasetList(idx).ShiftL=ShiftL;
        DatasetList(idx).ShiftR=ShiftR;
        id=id+cnt;
        mother(mother>0)=mother(mother>0)+cnt;
        daughter1(daughter1>0)=daughter1(daughter1>0)+cnt;
        daughter2(daughter2>0)=daughter2(daughter2>0)+cnt;
        tscnt=idx*ones(size(id));
        cnt0=0;
        for i=1:numel(id)
            if cycle{i}>9                                          % Only care about nc10
                cnt0=cnt0+1;
                datamat(cnt+cnt0).cycle=cycle{i};                  % Cell cycle 1
                datamat(cnt+cnt0).id=id(i);                        % nuclei ID (note, different index)
                datamat(cnt+cnt0).oldid=id(i)-cnt;                 % old nuclei ID, as in the unmerged data file
                datamat(cnt+cnt0).Intensity=Itmp{i};               % nuclei's spot intensity
                datamat(cnt+cnt0).time=trec{i};                    % time corresponding to intensity
                datamat(cnt+cnt0).Imax=Imax;                       % embryo maximum spot intensity
                datamat(cnt+cnt0).dt=dt;                           % time resolution
                datamat(cnt+cnt0).AdjustedIntensity=Itmp2{i};      % nuclei's adjusted spot intensity
                datamat(cnt+cnt0).Adjustedtime=trec2{i};           % time corresponding to intensity
                datamat(cnt+cnt0).Feature=fea{i};                  % extracted features
                datamat(cnt+cnt0).Feature_store=fea{i};            % extracted features
                datamat(cnt+cnt0).x=x(i)+shift_EL(idx)/100;        % nuclei position (mean position along AP axis)
                datamat(cnt+cnt0).y=y(i);                          % nuclei position (mean position along Y axis)
                datamat(cnt+cnt0).z=z(i);                          % nuclei position (position over time)
                datamat(cnt+cnt0).xrec=xrec{i}+shift_EL(idx)/100;  % nuclei position (position over time)
                datamat(cnt+cnt0).yrec=yrec{i};                    % nuclei position (position over time)
                datamat(cnt+cnt0).zrec=zrec{i};                    % nuclei position (position over time)
                datamat(cnt+cnt0).sizerec=sizerec{i};              % nuclei size (absolute)
                datamat(cnt+cnt0).xspotrec=xspotrec{i};            % spot position (relative to nuclei position over time)
                datamat(cnt+cnt0).yspotrec=yspotrec{i};            % spot position (relative to nuclei position over time)
                datamat(cnt+cnt0).mother=mother(i);                % Mother ID
                datamat(cnt+cnt0).daughter1=daughter1(i);          % Daughter1 ID
                datamat(cnt+cnt0).daughter2=daughter2(i);          % Daughter2 ID
                datamat(cnt+cnt0).tscnt=idx;                       % Index of time series
            end
        end
        cnt=cnt+cnt0;
    end
    close(f);
end