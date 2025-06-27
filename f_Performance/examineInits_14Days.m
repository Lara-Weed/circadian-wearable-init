%% Examine Num Inits
% Lara Weed
% 5 Apr 2025

%% Load Data - Select out 14 days only
basepath = '/Users/lara/Library/CloudStorage/OneDrive-Stanford/Research/Projects/Light Circadian/Data/InitsDays/StH';

rawData = load('/Users/lara/Library/CloudStorage/OneDrive-Stanford/Research/Projects/Light Circadian/Data/Preprocessed/subsetData_20241003.mat');

d = dir(basepath);

fn = {d.name}';
fn = fn(contains(fn,'.mat'));

noSubset = [1:11,6*2.^[1:7],1440];%[1:200];%unique(round(1440./(2.^[0:.25:10])));%[1:11,6*2.^[1:7],1440];
noTimes = 100;

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

    uDays = 14;%unique(T.DayfromDLMO);

    T.DLMO = T.DLMO + hours(1/3); % correction to 6.667 hour offset
    %betweeen CBTmin and DLMO instead of 7 hours

    T.DLMOrel = hour(T.DLMO) + minute(T.DLMO)./60 + second(T.DLMO)./3600 - TrueDLMO_hr;
    T.DLMOrel(T.DLMOrel>12) = T.DLMOrel(T.DLMOrel>12) - 24; 
    T.DLMOrel(T.DLMOrel<-12) = T.DLMOrel(T.DLMOrel<-12) + 24; 

    T2 = T(T.DayfromDLMO == 14,:);

    uStates = unique(T(:,2:3));

    if ~isempty(T2.DLMOrel)

        for k = 1:length(noSubset)
            numInits = noSubset(k);
            fprintf('      Subset: %d\n',numInits)
            for q = 1:noTimes
                EvalPoints = randsample(T2.DLMOrel,numInits,true);
    
                [f1,xf1] = kde(EvalPoints,"EvaluationPoints",linspace(-12,12,97));
                [f2,xf2] = kde(EvalPoints);
                NLL = -log(f1(xf1 == 0)./(sum(f1)) + 10e-6);
                MLE = xf2(find(f2 == max(f2),1));
                RMSE = sqrt(MLE.^2);
                [mu, K, Entropy] = estimateVonMisesParamsWithEntropy(2*pi*EvalPoints./24, 288);
            
                Mu = 24*mu./(2*pi);
    
                tt = [tt;i,14,noSubset(k),q,NLL,MLE,RMSE,Mu,K,Entropy];
            end
        end
    end
end

Metrics = {'NLL','MLE','RMSE','Mu','K','Entropy'};

%% Consolidate and plot

uinits = unique(tt(:,3));

meansRS = nan(length(uinits),6);
STDsRS = nan(length(uinits),6);

meansSW = nan(length(uinits),6);
STDsSW = nan(length(uinits),6);

for i = 1:length(uinits)

    init_ind = tt(:,3) == uinits(i);

    RS_ind = tt(:,1)<57;
    
    urs = unique(tt(RS_ind,1));
    ttrs = [];
    for kkk = 1:length(urs)
        ttrs =[ttrs ;mean(tt(init_ind & RS_ind & tt(:,1)==urs(kkk),5:end))];
    end

    usw = unique(tt(~RS_ind,1));
    ttsw = [];
    for kkk = 1:length(usw)
        ttsw =[ttsw ;mean(tt(init_ind & ~RS_ind & tt(:,1)==usw(kkk),5:end))];
    end

    meansRS(i,:) = mean(ttrs);

    meansSW(i,:) = mean(ttsw);

    STDsRS(i,:) = std(ttrs);

    STDsSW(i,:) = std(ttsw);
end

n_rs = length(unique(tt(tt(:,1)<57,1)));
n_sw = length(unique(tt(tt(:,1)>=57,1)));

% meansRS(isinf(meansRS)) = 10000;
% meansSW(isinf(meansSW)) = 10000;

%%
fig = figure('Renderer','painters','Position',[100 100 700 750]);%figure(Renderer="painters",Position=[500 500 550 500]);
for i = 1:size(meansRS,2)
    ax(i) = subplot(3,2,i);
    hold on

    % RS group (black)
    mean_RS = meansRS(:,i);
    sem_RS = STDsRS(:,i);% ./ sqrt(n_rs);
    upper_RS = mean_RS + sem_RS;
    lower_RS = mean_RS - sem_RS;

    fill([uinits; flipud(uinits)], [lower_RS; flipud(upper_RS)], ...
         'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    plot(uinits, mean_RS, 'k', 'LineWidth', 2)

    % SW group (red)
    mean_SW = meansSW(:,i);
    sem_SW = STDsSW(:,i);% ./ sqrt(n_sw);
    upper_SW = mean_SW + sem_SW;
    lower_SW = mean_SW - sem_SW;

    fill([uinits; flipud(uinits)], [lower_SW; flipud(upper_SW)], ...
         'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    plot(uinits, mean_SW, 'r', 'LineWidth', 2)

    ylabel(Metrics{i})
    if i >= 5
        xlabel('Initializations')
    end
    grid on
    % if i == 6
    %     legend({'', 'RS', '', 'SW'}, 'Location', 'best')
    % end
    if i == 5
        ylim([0 200])
    end
    %set(gca, 'XScale', 'log','fontweight','bold','fontsize',12)
    set(gca,'fontweight','bold','fontsize',15)
    xlim([min(uinits) max(uinits)])
end

linkaxes(ax, 'x');
xlim([1 1440])
%xlim([1 200])

savefig(fig,"Inits1440_V3_6.67hours.fig")
saveas(gcf,'Inits1440_V3_6.67hours.png')
close(fig)

%%
% figure(Renderer="painters",Position=[500 500 600 350])
% for i = 1:size(meansRS,2)
%     ax(i) = subplot(3,2,i);
%     hold on
% 
%     % RS group (black)
%     errorbar(uinits, meansRS(:,i), STDsRS(:,i), ...
%         'k', 'LineWidth', 2)
% 
%     % SW group (red)
%     errorbar(uinits, meansSW(:,i), STDsSW(:,i), ...
%         'r', 'LineWidth', 2)
% 
%     ylabel(['Metric ' num2str(i)])
%     if i >= 5
%         xlabel('Subset size (n)')
%     end
%     grid on
%     if i == 1
%         legend({'RS', 'SW'}, 'Location', 'best')
%     end
%     set(gca, 'fontweight', 'bold', 'fontsize', 12)
%     xlim([min(uinits) max(uinits)])
% end
% 
% linkaxes(ax, 'x');
