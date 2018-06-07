function [datamat,DatasetList]=local_extract(DatasetList)

    % Extracting data from the movies and combine them into one
    datamat=struct('id',[],'Intensity',{},'Feature',[],'x',[],'y',[],'border',[],'mother',[],'daughter1',[],'daughter2',[],'tscnt',[],'cycle',[],'tlen',[],'start',[]);
    outtab=[];
    cnt=0;
    % Load the movies
    for idx=1:numel(DatasetList)
        fullpath=DatasetList(idx).Path;
        [Itmp,Itmp2,trec,trec2,fea,x,y,xrec,yrec,sizerec,xspotrec,yspotrec,id,mother,daughter1,daughter2,cycle,dt,Imax,BG,ShiftL,ShiftR,APole]=extract_movie(fullpath,0,DatasetList(idx).BG_man);
        DatasetList(idx).BG=BG;
        DatasetList(idx).dt=dt;
        DatasetList(idx).APole=APole;
        DatasetList(idx).ShiftL=ShiftL;
        DatasetList(idx).ShiftR=ShiftR;
        id=id+cnt;
        mother(mother>0)=mother(mother>0)+cnt;
        daughter1(daughter1>0)=daughter1(daughter1>0)+cnt;
        daughter2(daughter2>0)=daughter2(daughter2>0)+cnt;
        tscnt=idx*ones(size(id));
        for i=1:numel(id)
            datamat(id(i)).id=id(i);                        % nuclei ID
            datamat(id(i)).oldid=id(i)-cnt;                 % old nuclei ID, as in the unmerged data file
            datamat(id(i)).Intensity=Itmp{i};               % nuclei's spot intensity
            datamat(id(i)).AdjustedIntensity=Itmp2{i};      % nuclei's adjusted spot intensity
            datamat(id(i)).Imax=Imax;                       % embryo maximum spot intensity
            datamat(id(i)).dt=dt;                           % time resolution
            datamat(id(i)).time=trec{i};                    % time corresponding to intensity
            datamat(id(i)).Adjustedtime=trec2{i};           % time corresponding to intensity
            datamat(id(i)).Feature=fea{i};                  % extracted features
            datamat(id(i)).Feature_store=fea{i};            % extracted features
            datamat(id(i)).x=x(i);                          % nuclei position (relative along AP axis)
            datamat(id(i)).y=y(i);                          % nuclei position (relative along Y axis)
            datamat(id(i)).xrec=xrec{i};                    % nuclei position (absolute)
            datamat(id(i)).yrec=yrec{i};                    % nuclei position (absolute)
            datamat(id(i)).sizerec=sizerec{i};              % nuclei size (absolute)
            datamat(id(i)).xspotrec=xspotrec{i};            % spot position (relative to nuclei)
            datamat(id(i)).yspotrec=yspotrec{i};            % spot position (relative to nuclei)
            datamat(id(i)).mother=mother(i);                % Mother ID
            datamat(id(i)).daughter1=daughter1(i);          % Daughter1 ID
            datamat(id(i)).daughter2=daughter2(i);          % Daughter2 ID
            datamat(id(i)).cycle=cycle{i};                  % Cell cycle 1
            datamat(id(i)).tscnt=idx;                       % Index of time series
        end
        cnt=cnt+numel(id);
    end
end