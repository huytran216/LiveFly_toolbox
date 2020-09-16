function [res,x,time]=extract_feature(x,time,threshold,dt,Imax,censored,trange,z,tburst)
    % Feature extraction
    % Input:
    %   x: spot intensity over time during phase S.
    %   time: time to reach steady state
    %   type: type of feature to be extracted
    %       % 1. ON/OFF
    %       % 2. Relative First Activation time
    %       % 3. Relative Last Activation time
    %       % 4. Relative Total active duration
    %       % 5. Relative Total spot duration of existence 
    %       % 6. Relative Time to max intensity
    %       % 7. Brightest spot intensity
    %       % 8. Mean Intensity produced (absolute time) (including non-expressing nuclei and non-expressing frame)
    %       % 9. Total Intensity produced (including non-expressing nuclei and non-expressing frame)
    %       % 10. Total interphase duration
    %       % 11. Absolute First Activation time
    %       % 12. Absolute Last Activation time
    %       % 13. Absolute Total active duration
    %       % 14. Absolute Total spot duration of existence
    %       % 15. Absolute Time to max intensity
    %       % 16. Mean Spot Intensity (including non-expressing nuclei and non-expressing frame)
    %       % 17. Mean Spot Intensity (only expressing nuclei (full traces) including non-expressing frames)
    %       % 18. Mean Spot Intensity (only expressing nuclei (trimmed traces) including non-expressing frames)
    %       % 19. Mean Spot Intensity (only expressing nuclei and expressing frames)
    %       % 20. Mean PSpot (including non-expressing nuclei)
    %       % 21. Mean PSpot (only non-expressing nuclei (full traces))
    %       % 22. Mean PSpot (only non-expressing nuclei (trimmed traces))
    %       % SEE CREATE_FEATURE_LABEL.M FOR UPDATES
    %   threshold: threshold to define if x is significant or not
    %   dt: time resolution
    %   Imax: maximum intensity
    %   censored: is the time trace sensored or not?
    %       % -1: invalid time traces
    %       % 0: not censored
    %       % 1: right censored
    %       % 2: left censored
    %   trange: range of valid time (used only for trimming)
    % Output:
        % res: array of features to be extracted
    %% Set initial parameters if missing
    NFeature = 24;
    res=-ones(1,NFeature);
    if nargin==0
        res=-ones(1,NFeature);
        return;
    end
    if ~exist('threshold','var')
        threshold=0;
    end
    if ~exist('dt','var')
        dt=1;
    end
    if ~exist('Imax','var')
        Imax=1;
    end
    if ~exist('censored','var')
        censored=0;
    end
    if ~exist('tburst','var')
        tburst = 300;
    end
    % Fit invalid cell then don't do anything
    if (censored<0)
        res=-ones(1,NFeature);
        return;
    end
    % Remove too dim spots
    x(x<threshold)=0;
    % Remove single-frame spots
        xtmp=[0 x>0 0];
        for i=1:numel(x)
            if norm(xtmp([i i+1 i+2]) - [0 1 0])==0
                x(i)=0;
            end
        end
    %% Trim trace if necessary
    tinterphase=length(x);
    if ~exist('z','var')
        z = ones(1,numel(x));
    end
    if ~exist('trange','var')
        trange = [0 10000];
    end
    lower_valid = trange(1);
    upper_valid = trange(2);
    % if time window is non-trivial
    x_ori = x;
    time_ori = time;
    if trange(2)<=tinterphase*dt
        x=x(round(trange(1)/dt)+1:round(trange(2)/dt));
        time=time(round(trange(1)/dt)+1:round(trange(2)/dt));
    else
        upper_valid = tinterphase*dt-100;
    end
    if trange(1)==0
        lower_valid = 200;
    end
    %If nuclei out of bound >> ignore
    if any(z((time>lower_valid)&(time<upper_valid))<0)
        res=-ones(1,NFeature);
        return;
    end
    %% Extract features
        
        % P_ON/P_OFF
                pthresh=0.01;    % We consider the gene is ON if the spot exists for more than 1% of the observation window
                nthresh=3;       % We consider the gene is ON if the spot exists for more 3 frames
                if (double(sum(x_ori>0)/numel(x_ori)>=pthresh))&(sum(x_ori>0)>=nthresh)
                    on_ori=1;
                else
                    on_ori=0;
                    x_ori=x_ori*0;
                end
                if (double(sum(x>0)/numel(x)>=pthresh))&(sum(x>0)>=nthresh)
                    res(1)=1;
                else
                    res(1)=0;
                    x=x*0;
                end
        % T_0
                if res(1)
                    res(2)=find(x,1,'first');
                else
                % If not found, we have a censored value of t0
                    res(2)=length(x);
                end
        % T_END
                if res(1)
                    res(3)=(length(x)-find(x,1,'last'));
                else
                    res(3)=length(x);
                end
        % T_ACTIVE
                res(4)=length(x)-res(3)-res(2)+1;
                if res(4)<0
                    res(4)=0;
                end
        % T_EXIST Total spot duration of existence 
                res(5)=sum(x>0);
        % Time to max intensity
                Ith=Imax*2/3;
                tmp2=find(x>=Ith,1,'first');
                % If not found, we have a censored value of t0
                if tmp2
                    res(6)=tmp2;
                else
                    res(6)=length(x);
                end
        % Brightest spot intensity
                res(7)=max(x);
        % Mean RNA synthesis rate (relative time)
                if res(4)>0
                    res(8)=sum(x)/res(4)*length(x);
                else
                    res(8)=0;
                end
        % Total RNA synthesized
                if res(4)>0
                    res(9)=sum(x)*dt;
                else
                    res(9)=0;
                end
        % Total duration of interphase
                res(10)=tinterphase*dt;
        % Normalize the time feature
                if res(4)>=0
                    res(11:15)=res(2:6)*dt;
                    res(2:6)=res(2:6)/length(x);
                end
        % Mean Spot intensity (including non-expressing nuclei and non-expressing frame)
                if res(4)>0
                    res(16)=mean(x);
                else
                    res(16)=0;
                end
        % Mean Spot Intensity (only expressing nuclei (full trace)and non-expressing frames)
                if on_ori>0
                    t0 = find(x_ori>0,1,'first');
                    tend = min(round(trange(2)/dt),numel(x_ori));
                    [~,ia] = intersect(time,time_ori(t0:tend));
                    if res(4)>0
                        res(17)=mean(x(ia));
                    else
                        res(17)=0;
                    end
                else
                    res(17)=NaN;
                end
        % Mean Spot Intensity (only expressing nuclei (trimmed trace)and non-expressing frames)
                if res(4)>0
                    res(18)=mean(x);
                else
                    res(18)=NaN;
                end
        % Mean Spot Intensity (only expressing nuclei and expressing frames)
                if res(4)>0
                    res(19)=mean(x(x>0));
                else
                    res(19)=NaN;
                end
        % Mean PSpot:
                if res(4)>0
                    res(20)=mean(x>0);
                else
                    res(20)=0;
                end
        % Mean PSpot (only expressing nuclei (full trace)and non-expressing frames)
                if on_ori>0
                    t0 = find(x_ori>0,1,'first');
                    tend = min(round(trange(2)/dt),numel(x_ori));
                    [~,ia] = intersect(time,time_ori(t0:tend));
                    if res(4)>0
                        res(21)=mean(x(ia)>0);
                    else
                        res(21)=0;
                    end
                else
                    res(21)=NaN;
                end
        % Mean PSpot (only expressing nuclei (trimmed trace)and non-expressing frames)
                if res(4)>0
                    res(22)=mean(x>0);
                else
                    res(22)=NaN;
                end
        % 23. Mean Spot Intensity (after 1st spot appearance, during the 1st burst only)
        % 24. Burst's fraction of time (after 1st spot appearance, duration of 1st burst divided by observation time window)
                if on_ori>0
                    % Extract first spot appearance
                    t0 = find(x_ori>0,1,'first');
                    tend = min(round(trange(2)/dt),numel(x_ori));
                    if (tend-t0)>=tburst/dt
                        % Last spot appearance
                        toff = find(x_ori(t0:tend)==0,1,'first')-1;
                        if toff*dt<tburst
                            x_trim = x_ori(t0:t0+toff-1);
                            res(24) = toff/round(tburst/dt);
                        else
                            x_trim = x_ori(t0:t0+round(tburst/dt));
                            res(24) = 1;
                        end
                        res(23) = mean(x_trim);
                    else
                        % Not enough remaining time to analyze
                        res(23)=NaN;
                        res(24)=NaN;
                    end
                else
                    res(23)=NaN;
                    res(24)=NaN;
                end   
    end