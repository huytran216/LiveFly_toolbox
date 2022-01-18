addpath 'C:\Users\bvhutr\Dropbox (UMR3664)\backup\Memory_analysis\sample_calculation'
addpath 'C:\Users\bvhutr\Dropbox (UMR3664)\backup\Memory_analysis\sample_calculation\julia'
addpath 'D:\Users\HuyT\Dropox_curie\Dropbox (UMR3664)\backup\Memory_analysis\sample_calculation'
addpath 'D:\Users\HuyT\Dropox_curie\Dropbox (UMR3664)\backup\Memory_analysis\sample_calculation\julia'
datafilename = 'mounia_data1';
H=0.01;      % Hill coefficient
x0 = 0;   % Boundary position
wd=5
datafolder = [datafilename '_13vs14'];
fea_range = [1:24];
min_alpha = [1:24];
alpha_prob_rec = [];
alpha_prob_rec_mix_within = [];
alpha_prob_rec_mix_embryo = [];
to_plot = numel(fea_range)==1;
for fea1 = fea_range
    for fea2 = fea_range
        load([datafolder '/data_' num2str(fea1) '_' num2str(fea2) '_' num2str(wd) '.mat'],...
            'pval_mix_within','pval_mix_all','pval','Nsample','pos_range');
        Nmix_within = size(pval_mix_within,1);
        Nmix_all = size(pval_mix_all,1);
        [probability_alpha,alpha_range,p_rec,pos_range] = probability_profile(...
            [pval(1,:)<0.05;...
            pval_mix_within(:,:)<0.05;...
            pval_mix_all(:,:)<0.05],Nsample,'gamma',pos_range,H,20,x0);
        alpha_prob_rec(fea1,fea2,:) = probability_alpha(1,:);
        alpha_prob_rec_mix_within(fea1,fea2,:) = mean(probability_alpha(2:Nmix_within+1,:),1);
        alpha_prob_rec_mix_all(fea1,fea2,:) = mean(probability_alpha(Nmix_within+2:end,:),1);
        if find(alpha_prob_rec(fea1,fea2,:)>0.05,1,'last')
            min_alpha(fea1,fea2)=alpha_range(find(alpha_prob_rec(fea1,fea2,:)>0.05,1,'last'));
        else
            min_alpha(fea1,fea2)=0;
        end
        
        if find(alpha_prob_rec_mix_within(fea1,fea2,:)>0.05,1,'last')
            min_alpha_mix_within(fea1,fea2)=alpha_range(find(alpha_prob_rec_mix_within(fea1,fea2,:)>0.05,1,'last'));
        else
            min_alpha_mix_within(fea1,fea2)=0;
        end
        if find(alpha_prob_rec_mix_all(fea1,fea2,:)>0.05,1,'last')
            min_alpha_mix_all(fea1,fea2)=alpha_range(find(alpha_prob_rec_mix_all(fea1,fea2,:)>0.05,1,'last'));
        else
            min_alpha_mix_all(fea1,fea2)=0;
        end
        [fea1 fea2 min_alpha(fea1,fea2)]
        % Plot difference in alpha inference:
        if to_plot
            plot(alpha_range, probability_alpha(2:Nmix_within+1,:),'color','r','LineWidth',1);  hold on;
            plot(alpha_range, probability_alpha(Nmix_within+2:end,:),'color','g','LineWidth',1); 
            plot(alpha_range, probability_alpha(1,:),'color','b','LineWidth',2);
            xlabel('alpha');
            ylabel('Probability of data');
            hold off;
        end
    end
end
%% Save the data
if numel(fea_range)>10
    save([datafolder '/min_alpha'],'min_alpha','alpha_prob_rec',...
        'min_alpha_mix_within','alpha_prob_rec_mix_within',...
        'min_alpha_mix_all','alpha_prob_rec_mix_all',...
        'alpha_range','fea_range','x0','H');
end
%% Plot the new found map
figure;
subplot(121);
surf(min_alpha);
colorbar
caxis([0 1]);
ylabel('Feature of mother');
xlabel('Feature of daughter');
tmp = sort(unique(min_alpha(:))); tmp(end-1)
subplot(132);
surf(min_alpha_mix_within);
colorbar
caxis([0 1]);
ylabel('Feature of mother');
xlabel('Feature of daughter');
tmp = sort(unique(min_alpha_mix_within(:))); tmp(end-1)
subplot(133);
surf(min_alpha_mix_all);
colorbar
caxis([0 1]);
ylabel('Feature of mother');
xlabel('Feature of daughter');
tmp = sort(unique(min_alpha_mix_all(:))); tmp(end-1)