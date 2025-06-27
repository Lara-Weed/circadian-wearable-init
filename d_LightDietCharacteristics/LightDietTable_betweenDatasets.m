%% Light Datasets Characteristics Tables
% Lara Weed
% 5 Feb 2025

%% Load Data

% Light data
load('/Users/lara/Library/CloudStorage/OneDrive-Stanford/Research/Projects/Light Circadian/Data/Preprocessed/subsetData_20241003.mat')
%fT = readtable('/Users/lara/Library/CloudStorage/OneDrive-Stanford/Research/Projects/Light Circadian/Data/FPCA/fpca output (1).xlsx');
load('/Users/lara/Library/CloudStorage/OneDrive-Stanford/Research/Projects/Light Circadian/Data/FPCA/TT_forCorr_fPCAResults.mat');

%%

sub = fields(T);

varNames = {'RecordingLength','TAT0days','TAT1000days','TAT0prct','TAT1000prct','IVlight','ISlight','MeanLux','MedianLux','lux95prct','PC1','PC2','PC3','PC4'};

lightData = nan(length(sub),length(varNames));

for i = 1:length(sub)
    
    epoch = median(diff(T.(sub{i}).t));

    % Recording Length
    lightData(i,1) = days(T.(sub{i}).t(end) - T.(sub{i}).t(1));
    
    % TAT>0 lux, days
    lightData(i,2) = days(sum(T.(sub{i}).I>0)*epoch);

    % TAT>1000 lux, days
    lightData(i,3) = days(sum(T.(sub{i}).I>1000)*epoch);

    % TAT>0 lux, %
    lightData(i,4) = lightData(i,2)./lightData(i,1); 

    % TAT>1000 lux, %
    lightData(i,5) = lightData(i,3)./lightData(i,1); 

    % IV light
    [IV,IS] = calcIVIS(T.(sub{i}).I,T.(sub{i}).t);
    lightData(i,6) = IV; 

    % IS Light
    lightData(i,7) = IS;

    % Mean Lux
    lightData(i,8) = mean(T.(sub{i}).I);

    % Median lux
    lightData(i,9) = median(T.(sub{i}).I);

    % 95th percentile lux
    lightData(i,10) = prctile(T.(sub{i}).I,95);

    if i < 88
        % PC1
        lightData(i,11) = TT.PC1(i,14);

        % PC2
        lightData(i,12) = TT.PC2(i,14);

        % PC3
        lightData(i,13) = TT.PC3(i,14);
    
        % PC4
        lightData(i,14) = TT.PC4(i,14);
    end

end

LD = [table(sub),array2table(lightData,"VariableNames",varNames)];

%% Compute Summary data
RS = LD(contains(LD.sub,'RS'),:);

SW = LD(~contains(LD.sub,'RS'),:);


for i = 1:length(varNames)
    H1 = lillietest(RS.(varNames{i}));

    % if H1 == 1
    %     fprintf('%s - %.2f(%.2f)\n',varNames{i},mean(RS.(varNames{i})),std(RS.(varNames{i})))
    % else
         fprintf('%s - %.2f[%.2f-%.2f]\n',varNames{i},nanmedian(RS.(varNames{i})),prctile(RS.(varNames{i}),[25,75]))
    % end

     H2 = lillietest(SW.(varNames{i}));

    % if H2 == 1
    %     fprintf('%s - %.2f(%.2f)\n',varNames{i},mean(SW.(varNames{i})),std(SW.(varNames{i})))
    % else
         fprintf('%s - %.2f[%.2f-%.2f]\n',varNames{i},nanmedian(SW.(varNames{i})),prctile(SW.(varNames{i}),[25,75]))
    % end

    if H1 == 1 && H2 == 1
        % Both normal -> use t-test
        [~, p] = ttest2(RS.(varNames{i}), SW.(varNames{i}));
        test_used = 't-test';
    else
        % At least one non-normal -> use Mann-Whitney U test
        p = ranksum(RS.(varNames{i}), SW.(varNames{i}));
        test_used = 'Mann-Whitney U test';
    end

    fprintf('%s - %s: %.4f \n',varNames{i},test_used,p)
end








