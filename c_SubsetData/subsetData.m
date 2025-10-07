%% Subset data to only include data with the same sampling rate 

%% Initializations
% Lara Weed
% 13 Aug 2024

%% Load Data
% Raw Data
fullset = load('/Users/lara/Library/CloudStorage/OneDrive-Stanford/Research/Projects/Light Circadian/Data/Preprocessed/preprocessed_20240813.mat');

save_path = '/Users/lara/Library/CloudStorage/OneDrive-Stanford/Research/Projects/Light Circadian/Data/Preprocessed';

%%
Subjects = fields(fullset.T);

collectionLength = nan(length(Subjects),1);
timeBeforeDLMO = nan(length(Subjects),1);
epochs = nan(length(Subjects),1);

mv = [];
for i = 1:length(Subjects)
    collectionLength(i) = hours(fullset.T.(Subjects{i}).t(end) - fullset.T.(Subjects{i}).t(1));
    timeBeforeDLMO(i) = hours(fullset.T.(Subjects{i}).DLMO - fullset.T.(Subjects{i}).t(1));
    epochs(i) = seconds(nanmedian(diff(fullset.T.(Subjects{i}).t)));

    % only consider 30s epochs
    if epochs(i) == 30

        T.(Subjects{i}) = fullset.T.(Subjects{i});

        if contains(Subjects{i},'SW')
            fprintf('%s - max: %d\n',Subjects{i},max(T.(Subjects{i}).I))
            mv = [mv;max(T.(Subjects{i}).I)];
        end
        % Regualar schedule only handles light up to 2700
        T.(Subjects{i}).I(T.(Subjects{i}).I>2700) = 2700;
        T.(Subjects{i}).I(isnan(T.(Subjects{i}).I)) = 0; % code doesn't handle nans


    end
end

%save(fullfile(save_path,'subsetData_20241003.mat'),'T')
% 
% median(collectionLength(contains(fields(T),'RS')))./24
% 
% median(collectionLength(contains(fields(T),'SW')))./24
% 
% median(timeBeforeDLMO(contains(fields(T),'RS')))./24
% 
% median(timeBeforeDLMO(contains(fields(T),'SW')))./24
