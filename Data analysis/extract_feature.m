function [res,x,time]=extract_feature(x,time,threshold,dt,Imax,censored,tinterphase)
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
    %       % 8. Mean RNA produced
    %       % 9. Total RNA produced
    %       % 10. Total interphase duration
    %       % 11. Absolute First Activation time
    %       % 12. Absolute Last Activation time
    %       % 13. Absolute Total active duration
    %       % 14. Absolute Total spot duration of existence 
    %       % 15. Absolute Time to max intensity
    %   threshold: threshold to define if x is significant or not
    %   dt: time resolution
    %   Imax: maximum intensity
    %   censored: is the time trace sensored or not?
    %       % -1: invalid time traces
    %       % 0: not censored
    %       % 1: right censored
    %       % 2: left censored
    % Output:
        % res: array of features to be extracted
    %% Set initial parameters if missing
    res=-ones(1,15);
    if nargin==0
        res=-ones(1,15);
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
    if ~exist('tinterphase','var')
        tinterphase=length(x);
    end
    % Fit invalid cell then don't do anything
    if (censored<0)||(tinterphase<=0)
        res=-ones(1,15);
        return;
    end
    % If the trace is right censored, pad x with zero on the right
    if censored==1
        x=[x zeros(1,round(tinterphase)-numel(x))];
        time=[time time(end)+dt*[1:(round(tinterphase)-numel(x))]];
    end
    if censored==2
        x=[zeros(1,round(tinterphase)-numel(x)) x];
        time=[time(1)+dt*[-(round(tinterphase)-numel(x)):-1]  time];
    end
    x(x<threshold)=0;
    %% Extract features
        % P_ON/P_OFF
                pthresh=0.01;    % We consider the gene is ON if the spot exists for more than 1% of the observation window
                if double(sum(x>0)/numel(x)>pthresh)
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
        % Total spot duration of existence 
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
        % Mean RNA synthesis rate
                if res(4)>0
                    res(8)=sum(x)/res(4);
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
                res(10)=(length(x))*dt;
        % Normalize the time feature
                if res(4)>=0
                    res(11:15)=res(2:6);
                    res(2:6)=res(2:6)/length(x);
                end
    end