function [res,x,time]=extract_feature(x,time,threshold,dt,Imax,censored,trange)
    % Feature extraction
    % Input:
    %   x: spot intensity over time during phase S.
    %   type: type of feature to be extracted
    %       % 1. ON/OFF
    %       % 2. Relative First Activation time
    %       % 3. Relative Last Activation time
    %       % 4. Relative Total active duration
    %       % 5. Relative Total spot duration of existence 
    %       % 6. Relative Time to max intensity
    %       % 7. Brightest spot intensity
    %       % 8. Mean RNA produced (absolute time)
    %       % 9. Total RNA produced
    %       % 10. Total interphase duration
    %       % 11. Absolute First Activation time
    %       % 12. Absolute Last Activation time
    %       % 13. Absolute Total active duration
    %       % 14. Absolute Total spot duration of existence
    %       % 15. Absolute Time to max intensity    
    %       % 16. Mean Spot Intensity 
    %       % 17. Mean PSpot (not ON)
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
    res=-ones(1,17);
    if nargin==0
        res=-ones(1,17);
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
    
    % Fit invalid cell then don't do anything
    if (censored<0)
        res=-ones(1,17);
        return;
    end
    % Remove too dim spots
    x(x<threshold)=0;
    %% Trim trace if necessary
    tinterphase=length(x);
    if exist('trange','var')
        if trange(2)<=tinterphase*dt
            x=x(round(trange(1)/dt)+1:round(trange(2)/dt));
            time=time(round(trange(1)/dt)+1:round(trange(2)/dt));
        end
    end
    %% Extract features
        % Remove single-frame spots
        xtmp=[0 x>0 0];
        for i=1:numel(x)
            if norm(xtmp([i i+1 i+2]) - [0 1 0])==0
                x(i)=0;
            end
        end
        % P_ON/P_OFF
                pthresh=0.01;    % We consider the gene is ON if the spot exists for more than 1% of the observation window
                nthresh=3;       % We consider the gene is ON if the spot exists for more 3 frames
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
        % Mean Spot intensity
                if res(4)>0
                    res(16)=mean(x);
                else
                    res(16)=0;
                end
        % Mean PSpot:
                if res(4)>0
                    res(17)=mean(x>0);
                else
                    res(17)=0;
                end
    end