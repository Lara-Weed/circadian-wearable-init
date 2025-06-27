%% Process Days and Inits Data
% Lara Weed
% 5 Apr 2025

%% Load Data
basepath = '/Users/lara/Library/CloudStorage/OneDrive-Stanford/Research/Projects/Light Circadian/Data/InitsDays/stH';

rawData = load('/Users/lara/Library/CloudStorage/OneDrive-Stanford/Research/Projects/Light Circadian/Data/Preprocessed/subsetData_20241003.mat');

d = dir(basepath);

fn = {d.name}';
fn = fn(contains(fn,'.mat'));

tt = [];
for i = 1:length(fn)

    if strcmp(basepath(end-3:end),'impM')
        gg = 13;
        ss = 21; 
    elseif strcmp(basepath(end-3:end),'/Exp')
        gg = 13;
        ss = 19;
    elseif strcmp(basepath(end-3:end),'H10K')
        gg = 13;
        ss = 22;
    else
        gg = 13;
        ss = 19;
    end
    fprintf('   %s - %d/%d\n',fn{i}(ss:end-gg),i,length(fn))
    
    load(fullfile(basepath,fn{i}))

    
    TrueDLMO = rawData.T.(fn{i}(ss:end - gg)).DLMO;
    TrueDLMO_hr = hour(TrueDLMO) + minute(TrueDLMO)./60 + second(TrueDLMO)./3600;

    T = extractTinitDays(DT);

    T = rmmissing(T);

    uDays = unique(T.DayfromDLMO);

    %T.DLMO = T.DLMO + hours(1/3); % correction to 6.667 hour offset
    %betweeen CBTmin and DLMO instead of 7 hours

    T.DLMOrel = hour(T.DLMO) + minute(T.DLMO)./60 + second(T.DLMO)./3600 - TrueDLMO_hr;
    T.DLMOrel(T.DLMOrel>12) = T.DLMOrel(T.DLMOrel>12) - 24; 
    T.DLMOrel(T.DLMOrel<-12) = T.DLMOrel(T.DLMOrel<-12) + 24; 

    for j = 1:length(uDays)
        [f1,xf1] = kde(T.DLMOrel(T.DayfromDLMO == uDays(j)),"EvaluationPoints",linspace(-12,12,97));
        [f2,xf2] = kde(T.DLMOrel(T.DayfromDLMO == uDays(j)));
        NLL = -log(f1(xf1 == 0)./(sum(f1)) + 10e-6);
        MLE = xf2(f2 == max(f2));
        RMSE = sqrt(MLE.^2);
        [mu, K, Entropy] = estimateVonMisesParamsWithEntropy(2*pi*T.DLMOrel(T.DayfromDLMO == uDays(j))./24, 288);
        
        Mu = 24*mu./(2*pi);

        tt = [tt;i,uDays(j),length(T.DLMOrel(T.DayfromDLMO == uDays(j))),NLL,MLE,RMSE,Mu,K,Entropy];
    end

end


%% Look at tt


uSubs = unique(tt(:,1));


TT.NLL = nan(length(uSubs),round(max(tt(:,2))));
TT.MLE = nan(length(uSubs),round(max(tt(:,2))));
TT.RMSE = nan(length(uSubs),round(max(tt(:,2))));
TT.Mu = nan(length(uSubs),round(max(tt(:,2))));
TT.K = nan(length(uSubs),round(max(tt(:,2))));
TT.Entropy = nan(length(uSubs),round(max(tt(:,2))));

figure
for i = 1:length(uSubs)

    if i <57
        color = 'k';
    else
        color = 'r';
    end

    ind = tt(:,1) == uSubs(i);

    ax(1) = subplot(6,1,1);
    plot(tt(ind,2),tt(ind,4),'o-',Color=color)
    hold on
    title('NLL')

    ax(2) = subplot(6,1,2);
    plot(tt(ind,2),tt(ind,5),'o-',Color=color)
    hold on
    title('MLE')

    ax(3) = subplot(6,1,3);
    plot(tt(ind,2),tt(ind,6),'o-',Color=color)
    hold on
    title('RMSE')

    ax(4) = subplot(6,1,4);
    plot(tt(ind,2),tt(ind,7),'o-',Color=color)
    hold on
    title('Mu')

    ax(5) = subplot(6,1,5);
    plot(tt(ind,2),tt(ind,8),'o-',Color=color)
    hold on
    title('K')

    ax(6) = subplot(6,1,6);
    plot(tt(ind,2),tt(ind,9),'o-',Color=color)
    hold on
    title('Entropy')

    for j = 1:round(max(tt(:,2)))

        subtt = tt(ind & tt(:,2)==j,4:9);

        if ~isempty(subtt)
            TT.NLL(i,j) = subtt(1);
            TT.MLE(i,j) = subtt(2);
            TT.RMSE(i,j) = subtt(3);
            TT.Mu(i,j) = subtt(4);
            TT.K(i,j) = subtt(5);
            TT.Entropy(i,j) = subtt(6);
        end
    end

end




%% Plots


% Define x-axis (assuming 45 days)
x = 1:45;

metrics = fields(TT);

figure('Renderer','painters','Position',[100 100 700 750])
for jj = 1:length(metrics)
    ax(jj) = subplot(3,2,jj);
    hold on
    
    metric = TT.(metrics{jj});
    
    RS_ind = false(length(uSubs),1);
    RS_ind(uSubs<=56) = true;

    % Group 1: black line
    mean1 = nanmedian(metric(RS_ind,:));
    p25_1 = prctile(metric(RS_ind,:), 25);
    p75_1 = prctile(metric(RS_ind,:), 75);

    keep_ind = sum(~isnan(metric(RS_ind,:)))>10;
    
    fill([x(keep_ind) fliplr(x(keep_ind))], [p25_1(keep_ind) fliplr(p75_1(keep_ind))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    plot(mean1(keep_ind), 'k', 'LineWidth', 2)
    
    % Group 2: red line
    mean2 = nanmedian(metric(~RS_ind,:));
    p25_2 = prctile(metric(~RS_ind,:), 25);
    p75_2 = prctile(metric(~RS_ind,:), 75);

    keep_ind2 = sum(~isnan(metric(~RS_ind,:)))>10;
    
    fill([x(keep_ind2) fliplr(x(keep_ind2))], [p25_2(keep_ind2) fliplr(p75_2(keep_ind2))], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    plot(mean2(keep_ind2), 'r', 'LineWidth', 2)
    
    % Labels and grid
    ylabel(metrics{jj})
    if jj >= 5 
        xlabel('Recording Days')
    end
    grid on
    % if jj ==    1
    %     legend({'', 'RS', '', 'SW'},'Location','best')
    % end
    set(gca,'fontweight','bold','fontsize',15)
    xlim([1 29])

    if jj == 1
        ylim([2 10])
    elseif jj== 2
         ylim([-5 5])
    elseif jj== 3
        ylim([0 7])
    elseif jj== 4
        ylim([-5 5.2])
    elseif jj== 5
        ylim([0 86])
    elseif jj== 6
        ylim([3.75 9])
    end

    
end
linkaxes(ax,'x');

xlim([1 14])

