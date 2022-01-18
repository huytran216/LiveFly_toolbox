%% Instruction
%   Load the data file directly to matlab before this analysis
%   Modify the cell cycle to be analyzed: cycle_m
%   Modify the feature to be tested
load('../../Data analysis/feature_label.mat');

datafilename_list={'Z2B6-near','H6B6-near','hb-vk33','B9-near','B6-near','B12-near','mounia_data1'};
%   Run the script
%% Load 
for datafilename_idx = 7
    datafilename = datafilename_list{datafilename_idx};
    load(['../../Data analysis/final_dataset/' datafilename '.mat'],'datamat','DatasetFeature');
    %% Define the plot parameters:
    cycle_m=13;             % Cycle of the mother cell.
                            % Cycle of daughter will be cycle_m+1
    align_border=0;         % Align the border or not

    fea_m_range=[1:24];      % Range of daughter feature to be compared - see feature label for reference
    fea_d_range=[1:24];      % Range of mother feature to be compared - see feature label for reference
    trimmed_Feature = datafilename_idx==7;
    %% Determine the samples for the question and features to look at
    wd=5;               % Set width of the scanning windows
    pos_start=[-15:15];
    question=1;         % Type of question to ask
        % Question 1. How to I look like my mother?
        % Question 2. How does my (better) daughter resemble me?
        % Question 3. How does my (worse) daughter resemble me?
        % Question 4. How do I determine the difference between my daughter ?
    plot_example=0;
    plot_sliding_window = 0;
    Nmix_embryo = 50;
    Nmix_all = 50;
    %% Extract the cycles and features of interest
    % Find the embryo of interest:
    ts_range = intersect(DatasetFeature(cycle_m-8).tsrec,DatasetFeature(cycle_m-8+1).tsrec);
    % Find all nuclei in two nuclei cycles of interest
    pos_m = find(ismember([datamat(:).tscnt],ts_range)&([datamat(:).cycle]==cycle_m));
    id_m = [datamat(pos_m).id];
    pos_d = find(ismember([datamat(:).tscnt],ts_range)&([datamat(:).cycle]==cycle_m+1));
    id_d = [datamat(pos_d).id];
    % Create the output folder:
    foldername=[datafilename '_' num2str(cycle_m) 'vs' num2str(cycle_m+1)];
    mkdir(foldername);
    %% Scan for feature pairs
    for fea_m=fea_m_range
        for fea_d=fea_d_range
            filename_img=[foldername '/img_' num2str(fea_m) '_' num2str(fea_d) '_' num2str(wd)];
            filename_data=[foldername '/data_' num2str(fea_m) '_' num2str(fea_d) '_' num2str(wd)];
            display([datafilename ': ' num2str(fea_m) ' ' num2str(fea_d)]);
            %% Select samples to address the question
            id1 = [];
            switch question
                case 1 % Question 1. How to I look like my mother?
                    id2=id_d;
                    id1= [datamat(pos_d).mother];
                case 2 % Question 2. How does my (better) daughter resemble me?
                    id1=[id_m];
                    id21=[datamat(id1).daughter1];
                    id22=[datamat(id1).daughter2];
                    tmp=(id21>0)&(id22>0);          % Filtering cells with less than two daughters
                    id1=id1(tmp);
                    id21=id21(tmp);
                    id22=id22(tmp);
                    fea21=arrayfun(@(x) subindex(datamat(x).Feature_store,fea_d),id21);
                    fea22=arrayfun(@(x) subindex(datamat(x).Feature_store,fea_d),id22);
                    id2=[id21;id22];
                    order=(fea21>=fea22)+1;
                    id2=arrayfun(@(x) id2(order(x),x),1:numel(order));
                case 3 % Question 3. How does my (worse) daughter resemble me?
                    id1=[id_m];
                    id21=[datamat(id1).daughter1];
                    id22=[datamat(id1).daughter2];
                    tmp=(id21>0)&(id22>0);          % Filtering cells with less than two daughters
                    id1=id1(tmp);
                    id21=id21(tmp);
                    id22=id22(tmp);
                    fea21=arrayfun(@(x) subindex(datamat(x).Feature_store,fea_d),id21);
                    fea22=arrayfun(@(x) subindex(datamat(x).Feature_store,fea_d),id22);
                    id2=[id21;id22];
                    order=(fea21<fea22)+1;
                    id2=arrayfun(@(x) id2(order(x),x),1:numel(order));
                case 4
                    id1=[id_m];
                    id21=[datamat(id1).daughter1];
                    id22=[datamat(id1).daughter2];
                    tmp=(id21>0)&(id22>0);          % Filtering cells with less than two daughters
                    id1=id1(tmp);
                    id21=id21(tmp);
                    id22=id22(tmp);
                    fea21=arrayfun(@(x) subindex(datamat(x).Feature_store,fea_d),id21);
                    fea22=arrayfun(@(x) subindex(datamat(x).Feature_store,fea_d),id22);
                    fea2=abs(fea21-fea22);
                    id2=id21;
            end

            % Making sure that mother/daughter cell exist            
            flt=(id1).*(id2)>0;
            id1=id1(flt);
            id2=id2(flt);
            ia=ismember(id1,id_m); % Check whether they are null cells
            id2=id2(ia);
            id1=id1(ia);
            ia=ismember(id2,id_d); % Check whether they are null cells
            id1=id1(ia);
            id2=id2(ia);
            % Convert to position
            [~,pos_id1] = ismember(id1,[datamat(:).id]);
            [~,pos_id2] = ismember(id2,[datamat(:).id]);
            % Extracting features for mother cells
            x1=[datamat(pos_id1).x]*100-50;
            y1=[datamat(pos_id1).y]*100;
            if trimmed_Feature
                fea1=arrayfun(@(x) subindex(datamat(x).Feature,fea_m),pos_id1);
            else
                fea1=arrayfun(@(x) subindex(datamat(x).Feature_store,fea_m),pos_id1);
            end
            ts1=[datamat(pos_id1).tscnt];
            % Extracting features for daughter cells - NO ALIGNMENT FOR NOW
            x2=[datamat(pos_id2).x]*100-50;
            y2=[datamat(pos_id2).y]*100;
            if question~=4
                if trimmed_Feature
                    fea2=arrayfun(@(x) subindex(datamat(x).Feature,fea_d),pos_id2);
                else
                    fea2=arrayfun(@(x) subindex(datamat(x).Feature_store,fea_d),pos_id2);
                end
            end
            ts2=[datamat(pos_id2).tscnt];
    %         %% Plot cell displacement analysis
    %         figure(1);
    %         idxtmp=(x1>-10)&(x1<-0);
    %         x1_=x1(idxtmp);x2_=x2(idxtmp);
    %         y1_=y1(idxtmp);y2_=y2(idxtmp);
    %         %subplot(121);
    %         plot(x2_-x1_,y2_-y1_,'x'); title('Nuclei displacement');
    %         xlabel('AP axis (%EL)');
    %         ylabel('side axis (%EL)');
            %% Probing for the relationship
            flt=(fea1>=0)&(fea2>=0);
            if plot_example
                example_relationship(x1(flt),y1(flt),fea1(flt),fea2(flt),wd,ts2(flt),[feature_label{fea_m} ' mother'],[feature_label{fea_d} ' daughter']);
            else
                [pos_range,tau,pval,Nsample,Nsample_ts,F,tau_mix_within,pval_mix_within,tau_mix_all,pval_mix_all]=probe_relationship(x1(flt),y1(flt),fea1(flt),fea2(flt),wd,ts2(flt),[feature_label{fea_m} ' mother'],[feature_label{fea_d} ' daughter'],pos_start,plot_sliding_window,Nmix_embryo, Nmix_all);
                % Plotting the results of memoryo analysis

                if numel(pos_range)>1
                    %figure(100);
                        subplot(131);   % Plot the Spearman tau between features
                        set(gca, 'ColorOrder', [1 0 0; 0 1 0; 0 0 1], 'NextPlot', 'replacechildren');
                        plot(pos_range',tau');
                        if numel(pos_range)>1
                            axis([pos_range(1) pos_range(end) -1 2]);
                        end
                        ylabel('Spearman''s rho');
                        title([feature_label{fea_m} ' vs ' feature_label{fea_d}]);
                        hold on;
                        for i=pos_range((pval(1,:)<=0.05)&(pval(2,:)>0.05)&(pval(3,:)>0.05))
                            plot(i,1.5,'r*');
                        end

                        h=legend('F_m vs F_d','F_d vs x','F_d vs y','Sign. Corr.');
                        set(h,'Location','NorthEast');
                        hold off;
        %                     subplot(132);   % Plot the Spearman tau`s p-value
        %                     set(gca, 'ColorOrder', [1 0 0; 0 1 0; 0 0 1], 'NextPlot', 'replacechildren');
        %                     plot(pos_range',pval');hold on;
        %                     plot(pos_range,pos_range*0+0.05,'--k');hold off;
        %                     axis([pos_range(1) pos_range(end) 0 1]);
                        title('rho''s pvalue');
                        subplot(132);   % Plot the sample numbers
                        plot(pos_range,Nsample);
                        title('Sample no.');
                        hold off;
                        subplot(133);   % Plot the information
                        plot(pos_range,mean(tau_mix_within,1),'r'); hold on;
                        plot(pos_range,mean(tau_mix_all,1),'b');
                        legend('Mix within','Mix all');
                        ylim([-1  2]);
                        hold off;
                        %plot(pos_range,Hm); hold on; plot(pos_range,Hd);plot(pos_range,Imd);
                        %plot(pos_range,Ibgmd,'--');
                        %for i=pos_range(pImd<=0.05)
                        %    plot([i i],[0 1],'r--');
                        %end
                        %legend('Hm','Hd','Imd','Ibg');
                        %hold off;
                    % Save the results as a movie
                    set(gcf, 'Position', [0 300 1200 300]);
                    savefig(gcf,[filename_img '.fig']);
                    if plot_sliding_window
                        F2=getframe(gcf);
                        imwrite(F2.cdata,[filename_img '.tif']);
                    end
                end
                for i=1:numel(F)
                    imwrite(F(i).cdata,[filename_img '.tif'],'writemode', 'append');
                end
                save([filename_data '.mat'],'flt','x1','fea1','fea2','ts2',...
                    'pos_range','tau','pval','Nsample','Nsample_ts','tau_mix_all','tau_mix_within','pval_mix_all','pval_mix_within');
            end
        end
    end
end