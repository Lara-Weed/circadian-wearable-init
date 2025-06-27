%% Compute light Diet Features & Performance
% Lara Weed
% 21 Apr 2025

%% Load Data
basepath = '/Users/lara/Library/CloudStorage/OneDrive-Stanford/Research/Projects/Light Circadian/Data/InitsDays/stH';

rawData = load('/Users/lara/Library/CloudStorage/OneDrive-Stanford/Research/Projects/Light Circadian/Data/Preprocessed/subsetData_20241003.mat');

d = dir(basepath);

fn = {d.name}';
fn = fn(contains(fn,'.mat'));

tt = [];
numDays = 14;
fourier_avg = nan(length(fn),1440);

fourier_day = nan(length(fn)*14,2880);

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
    uInit = flipud(unique(T.InitTime));

    T.DLMO = T.DLMO + hours(1/3); % correction to 6.667 hour offset
    %betweeen CBTmin and DLMO instead of 7 hours

    T.DLMOrel = hour(T.DLMO) + minute(T.DLMO)./60 + second(T.DLMO)./3600 - TrueDLMO_hr;
    T.DLMOrel(T.DLMOrel>12) = T.DLMOrel(T.DLMOrel>12) - 24; 
    T.DLMOrel(T.DLMOrel<-12) = T.DLMOrel(T.DLMOrel<-12) + 24; 

    LightDiet = nan(length(uDays),7); 
    PartData = rawData.T.(fn{i}(ss:end-gg));


    for j = 1:length(uDays)
        % Compute Circadian Model Results
        [f1,xf1] = kde(T.DLMOrel(T.DayfromDLMO == uDays(j)),"EvaluationPoints",linspace(-12,12,97));
        [f2,xf2] = kde(T.DLMOrel(T.DayfromDLMO == uDays(j)));
        NLL = -log(f1(xf1 == 0)./(sum(f1)) + 10e-6);
        MLE = xf2(f2 == max(f2));
        RMSE = sqrt(MLE.^2);
        [mu, K, Entropy] = estimateVonMisesParamsWithEntropy(2*pi*T.DLMOrel(T.DayfromDLMO == uDays(j))./24, 288);
        
        Mu = 24*mu./(2*pi);

        % Compute Light Diet Features

        day_ind =  PartData.t>=uInit(j) & PartData.t<= PartData.DLMO;

        LightCon = PartData.I(day_ind);
        tCon = PartData.t(day_ind);

        Luxmean = mean(LightCon);% Mean lux
        Luxmedian = median(LightCon);% Median Lux
        LuxAboveZero = sum(LightCon>0)/length(LightCon);% Percent Above zero lux
        LuxAbove1000 = sum(LightCon>1000)/length(LightCon);% Percent Above 1000 lux
        Lux95 = prctile(LightCon,95); % 95th Percentile lux
        
        [LuxIV,LuxIS] = calcIVIS(LightCon,PartData.t(day_ind));

        % table
        tt = [tt;i,uDays(j),length(T.DLMOrel(T.DayfromDLMO == uDays(j))),NLL,MLE,RMSE,Mu,K,Entropy,Luxmean,Luxmedian ,LuxAboveZero,LuxAbove1000 ,Lux95,LuxIV,LuxIS];


    end

end


%%

fT = readtable('/Users/lara/Library/CloudStorage/OneDrive-Stanford/Research/Projects/Light Circadian/Data/FPCA/fpca output (1).xlsx');

uSubs = unique(tt(:,1));


TT.NLL = nan(length(uSubs),round(max(tt(:,2))));
TT.MLE = nan(length(uSubs),round(max(tt(:,2))));
TT.RMSE = nan(length(uSubs),round(max(tt(:,2))));
TT.Mu = nan(length(uSubs),round(max(tt(:,2))));
TT.K = nan(length(uSubs),round(max(tt(:,2))));
TT.Entropy = nan(length(uSubs),round(max(tt(:,2))));

TT.Luxmean = nan(length(uSubs),round(max(tt(:,2))));
TT.Luxmedian = nan(length(uSubs),round(max(tt(:,2))));
TT.LuxAboveZero = nan(length(uSubs),round(max(tt(:,2))));
TT.LuxAbove1000 = nan(length(uSubs),round(max(tt(:,2))));
TT.Lux95 = nan(length(uSubs),round(max(tt(:,2))));
TT.LuxIV = nan(length(uSubs),round(max(tt(:,2))));
TT.LuxIS = nan(length(uSubs),round(max(tt(:,2))));

TT.PC1 = nan(length(uSubs),round(max(tt(:,2))));
TT.PC2 = nan(length(uSubs),round(max(tt(:,2))));
TT.PC3 = nan(length(uSubs),round(max(tt(:,2))));
TT.PC4 = nan(length(uSubs),round(max(tt(:,2))));


%figure
for i = 1:length(uSubs)

    % % if i <57
    % %     color = 'k';
    % % else
    % %     color = 'r';
    % % end
    % % 
    ind = tt(:,1) == uSubs(i);
    % % 
    % % ax(1) = subplot(6,1,1);
    % % plot(tt(ind,2),tt(ind,4),'o-',Color=color)
    % % hold on
    % % title('NLL')
    % % 
    % % ax(2) = subplot(6,1,2);
    % % plot(tt(ind,2),tt(ind,5),'o-',Color=color)
    % % hold on
    % % title('MLE')
    % % 
    % % ax(3) = subplot(6,1,3);
    % % plot(tt(ind,2),tt(ind,6),'o-',Color=color)
    % % hold on
    % % title('RMSE')
    % % 
    % % ax(4) = subplot(6,1,4);
    % % plot(tt(ind,2),tt(ind,7),'o-',Color=color)
    % % hold on
    % % title('Mu')
    % % 
    % % ax(5) = subplot(6,1,5);
    % % plot(tt(ind,2),tt(ind,8),'o-',Color=color)
    % % hold on
    % % title('K')
    % % 
    % % ax(6) = subplot(6,1,6);
    % % plot(tt(ind,2),tt(ind,9),'o-',Color=color)
    % % hold on
    % % title('Entropy')

    for j = 1:round(max(tt(:,2)))

        subtt = tt(ind & tt(:,2)==j,4:end);

        if ~isempty(subtt)
            TT.NLL(i,j) = subtt(1);
            TT.MLE(i,j) = subtt(2);
            TT.RMSE(i,j) = subtt(3);
            TT.Mu(i,j) = subtt(4);
            TT.K(i,j) = subtt(5);
            TT.Entropy(i,j) = subtt(6);

            TT.Luxmean(i,j) = subtt(7);
            TT.Luxmedian(i,j) = subtt(8);
            TT.LuxAboveZero(i,j) = subtt(9);
            TT.LuxAbove1000(i,j) = subtt(10);
            TT.Lux95(i,j) = subtt(11);
            TT.LuxIV(i,j) = subtt(12);
            TT.LuxIS(i,j) = subtt(13);

            if j == 14

                ind2 = fT.Var1>= numDays*(uSubs(i) - 1)+1 & fT.Var1<= numDays*(uSubs(i) - 1)+numDays;

                TT.PC1(i,j) = mean(fT.Var2(ind2)); %I_score(i,1);
                TT.PC2(i,j) = mean(fT.Var3(ind2));%I_score(i,2);
                TT.PC3(i,j) = mean(fT.Var4(ind2));%I_score(i,3);
                TT.PC4(i,j) = mean(fT.Var5(ind2));%I_score(i,4);
            end
        end
    end

end


%%

% Define x-axis (assuming 45 days)
x = 1:45;

Metrics = fields(TT);

figure
for jj = 1:13%length(Metrics)
    ax(jj) = subplot(7,2,jj);
    hold on
    
    metric = TT.(Metrics{jj});
    
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
    ylabel(Metrics{jj})
    if jj >= 5 
        xlabel('Recording Days')
    end
    grid on
    if jj ==    1
        legend({'', 'RS', '', 'SW'},'Location','best')
    end
    set(gca,'fontweight','bold','fontsize',12)
    xlim([1 29])

    if jj == 1
        ylim([0 12])
    elseif jj== 2
         ylim([-5 5])
    elseif jj== 3
        ylim([0 7])
    elseif jj== 4
        ylim([-10 10])
    %elseif jj== 5
    %    ylim([0 7])
    elseif jj== 6
        ylim([0 9])
    end

    
end
linkaxes(ax,'x');


%% Scatter with all days

% Define x-axis (assuming 45 days)
x = 1:45;

Metrics = fields(TT);

for jj = 7:length(Metrics)
    figure
    ax(1) = subplot(1,2,1);
    hold on
    
    metric = TT.(Metrics{jj});
    
    RS_ind = false(length(uSubs),1);
    RS_ind(uSubs<=56) = true;

    % Group 1: black line
    plot(metric(RS_ind,14),TT.RMSE(RS_ind,14),'ro')
    hold on
    plot(metric(~RS_ind,14),TT.RMSE(~RS_ind,14),'ko')
    
    % Labels and grid
    xlabel(Metrics{jj})
    
    ylabel('RMSE')

    grid on

    legend({ 'RS', 'SW'},'Location','best')

    set(gca,'fontweight','bold','fontsize',12)

  
    ax(2) = subplot(1,2,2);
    hold on
    
    metric = TT.(Metrics{jj});
    
    RS_ind = false(length(uSubs),1);
    RS_ind(uSubs<=56) = true;

    % Group 1: black line
    plot(metric(RS_ind,14),TT.Entropy(RS_ind,14),'ro')
    hold on
    plot(metric(~RS_ind,14),TT.Entropy(~RS_ind,14),'ko')
    
    % Labels and grid
    xlabel(Metrics{jj})
    
    ylabel('Entropy')

    grid on
    set(gca,'fontweight','bold','fontsize',12)
    
    linkaxes(ax,'x');
end





%% 

fieldsN = fieldnames(TT);
numFields = length(fieldsN);
nSubjects = size(TT.NLL, 1);
nModels = size(TT.NLL, 2);
nTotal = nSubjects * nModels;

% % Create a table with all variables vectorized
% dataTable = table;
% for i = 1:numFields
%     fieldName = fields{i};
%     dataTable.(fieldName) = reshape(TT.(fieldName), nTotal, 1);
% end
% 
% % Add subject and model identifiers
% [subjectGrid, modelGrid] = ndgrid(1:nSubjects, 1:nModels);
% dataTable.Subject = reshape(subjectGrid, nTotal, 1);
% dataTable.ModelID = reshape(modelGrid, nTotal, 1);
% 
% numericVars = dataTable(:, setdiff(fields, {'Subject', 'ModelID'}));
% [X, P] = corr(table2array(numericVars), 'rows', 'pairwise');
% %mask = triu(true(size(X)), 1);
% sigMask = P < 0.05;
% R_sig = X;
% R_sig(~sigMask) = NaN;
% %R_sig(mask) = NaN;
% 
% figure
% h = heatmap(fields, fields, R_sig, ...
%     'Colormap', jet, ...
%     'ColorLimits', [-1 1], ...
%     'Title', 'Significant Correlations Only');
% 
% 
% mdl = fitlm(dataTable, 'RMSE ~ LuxIV + Luxmean');
% disp(mdl);
% 
% lme = fitlme(dataTable, 'RMSE ~ LuxIV + Luxmean + (1|Subject)');
% disp(lme);
% 
% 
% 
% % Upper triangle mask (excluding diagonal)
% mask = triu(true(size(R)), 1);
% 

%%
% Define the performance metrics and light features of interest
perfMetrics = {'RMSE', 'Entropy'};
lightFeatures = {'Luxmean', 'Luxmedian', 'LuxAboveZero', 'LuxAbove1000', ...
                 'Lux95', 'LuxIV', 'LuxIS','PC1','PC2','PC3','PC4'};


% Combine selected fields
selectedFields = [perfMetrics, lightFeatures];

% Create reduced table
dataTable1 = table;
for i = 1:length(selectedFields)
    dataTable1.(selectedFields{i}) = reshape(TT.(selectedFields{i}), nTotal, 1);
end

% Compute correlation matrix and p-values
[R, P] = corr(table2array(dataTable1), 'rows', 'pairwise');

% Significance mask (only keep p < 0.05 in upper triangle)
sigMask =  triu(true(size(R)), 1);% & (P < 0.05) ;

% Keep only significant correlations
R_display = nan(size(R));
R_display(sigMask) = R(sigMask);

% Split indices
[~, idxPerf] = ismember(perfMetrics, selectedFields);
[~, idxLight] = ismember(lightFeatures, selectedFields);

% Create submatrix for heatmap: correlations between perf ↔ light features
R_sub = R_display(idxPerf, idxLight);
P_sub = P(idxPerf, idxLight);


figure('Renderer','painters','Position',[500 500 1000 350])
imagesc(R_sub)
% colormap(turbo)
% colorbar
% clim([-1 1])
% cb.Direction = 'reverse';
cmap = turbo;          % default 256‑by‑3 matrix
colormap(flipud(cmap)) % flip it top‑to‑bottom



cb = colorbar;         % optional, but helpful
clim([-1 1])
cb.Ticks = [-1 0 1];   % keep tick order logical
cb.Direction = 'reverse';


xticks(1:length(lightFeatures))
xticklabels(lightFeatures)
yticks(1:length(perfMetrics))
yticklabels(perfMetrics)
xtickangle(45)
% set(gca, 'XAxisLocation', 'top')
%title('Significant Correlations: Light Features vs. Performance Metrics')

% === Add r-values and significance annotations ===
hold on
for i = 1:length(perfMetrics)
    for j = 1:length(lightFeatures)
        r_val = R_sub(i, j);
        p_val = P_sub(i,j);

        if ~isnan(r_val)
            if p_val <= 0.05 && p_val >0.01
                % add superscript “n.s.”
                label = sprintf('%.2f^{\\it*}', r_val); 
            elseif p_val <= 0.01
                label = sprintf('%.2f^{\\it**}', r_val);
            else
                %label = sprintf('%.2f^{\\it n.s.}', r_val); 
                label = sprintf('%.2f', r_val);
            end

            % Plot text in center of each box
            text(j, i, label, 'HorizontalAlignment', 'center', ...
                'Color', [0.3 0.3 0.3], 'FontWeight', 'bold', 'FontSize', 14)
        end
    end
end
hold off

set(gca, 'fontweight','bold','fontsize',14)

