% Get shiftworker names

%% Load Data
load('subsetData_20241003.mat')

%% 
sub = fields(T);

sub_sw = sub(contains(sub,'RS'));

sub_names = [];
for i = 1:length(sub_sw)
    sub_names = [sub_names;{T.(sub_sw{i}).filename(1:end-4)}];

    
end

