%% Organize Data
% This code oragnaizes two raw datasets into a similar format for analysis.
% The datsets are of actigraphy data including light exposure recorded over
% approximately 1 week followed by an assessment of DLMO. The regular
% schedule dataset was obtained from the Zeitzer Circadian Research lab and
% the shiftwork dataset was obtained from Phillip Cheng.
% Lara Weed
% 14 Jan 2024

%% Ploting
isPlot_RS = 1;
isPlot_SW = 1;


%% Get paths
basepath = '/Users/lara/Library/CloudStorage/OneDrive-Stanford/Research/Projects/Light Circadian/Data';
save_fn = '/Users/lara/Library/CloudStorage/OneDrive-Stanford/Research/Projects/Light Circadian/Data/Organized';

SW = fullfile(basepath,'Raw/Shiftwork_Cheng');
RS = fullfile(basepath,'Raw/RegularSchedule_Zeitzer');

SW_data = fullfile(SW,'Actigraphy Data');
RS_data = fullfile(RS,'Raw');

%% Regular Schedule Data

% Load Meta data
SubInfo = readtable(fullfile(RS, 'Subject Info/subjectinfo.csv'));
fi = SubInfo.DateOfArrival + SubInfo.fi;
O_L = SubInfo.O_L;

% Load Data
dpath = fullfile(RS, 'Raw');
ddpath = dir(dpath);
fn = {ddpath.name}';
fn = fn(3:end);
Sub = {};
numDays = nan(length(fn),1);
missing = [];

for i = 1:length(fn)
    fprintf('%d - %s\n',i,fn{i})
    Subject = {sprintf('RS%02d',i)};

    dt.(Subject{1}).DLMO = fi(strcmp(SubInfo.Filename,fn{i}));
    dt.(Subject{1}).filename = fn{i};

    data = readtable(fullfile(dpath,fn{i}),'FileType','text');
    
    % Correct year timestamp
    t_raw = data.Date + data.Time;
    if year(data.Date)<2000
        t_raw = t_raw + years(2000);
    end
    
    dt.(Subject{1}).t = t_raw ;
    dt.(Subject{1}).I = data.Light;
    dt.(Subject{1}).act = data.Act;

    if isPlot_RS
        figure
        ax(1) = subplot(2,1,1);
        plot(dt.(Subject{1}).t,dt.(Subject{1}).I)
        hold on
        plot(dt.(Subject{1}).DLMO,200,'*','MarkerSize',10)
        title(sprintf('%d - %s\n',i,fn{i}))
        ax(2) = subplot(2,1,2);
        plot(dt.(Subject{1}).t,dt.(Subject{1}).act,'k')
        hold on
        plot(dt.(Subject{1}).DLMO,50,'*','MarkerSize',10)
        linkaxes(ax,'x')
    end
    
end

%% Shift work Data
% Load Meta data
SubInfo = readtable(fullfile(SW, 'DLMO Times.csv'));
visitDay = readtable(fullfile(SW, 'visit Dates.csv'));

%Create one table with participant DLMO
loopSubs = unique(SubInfo.Subject);
DLMO = NaT(length(loopSubs),1);
for q = 1:length(loopSubs)
    currSub = loopSubs{q};
    
    % determine date
    ind_visit = contains(visitDay.SubjectID,currSub(1:2)) & contains(visitDay.SubjectID,currSub(end-2:end));
    subVisitTable = visitDay(ind_visit,:);
    Visit1_Day = day(subVisitTable.Visit1Date);
    Visit1_Month = month(subVisitTable.Visit1Date);
    Visit1_Year = year(subVisitTable.Visit1Date);
    if Visit1_Year < 2000
        Visit1_Year = Visit1_Year + 2000;
    end

    visit1 = datetime(Visit1_Year,Visit1_Month,Visit1_Day);

    Visit2_Day = day(subVisitTable.Visit2Date);
    Visit2_Month = month(subVisitTable.Visit2Date);
    Visit2_Year = year(subVisitTable.Visit2Date);
    if Visit2_Year < 2000
        Visit2_Year = Visit2_Year + 2000;
    end

    visit2 = datetime(Visit2_Year,Visit2_Month,Visit2_Day);

    % Determine Time
    D = day(SubInfo.DLMO_time(q));
    M = month(SubInfo.DLMO_time(q));
    Y = year(SubInfo.DLMO_time(q));

    daySub = datetime(Y,M,D);

    DLMOTiming = SubInfo.DLMO_time(q) - daySub;

    % Create DLMO Datetime on correct date
    if strcmp(currSub,'PH011') || strcmp(currSub,'PH004') 
        DLMO(q) = visit2 + DLMOTiming ;
    else
        DLMO(q) = visit1 + DLMOTiming ;
    end
end

% Load Data
ddpath = dir(SW_data);
fn = {ddpath.name}';
fn = fn(3:end);
Sub = {};
numDays = nan(length(fn),1);
missing = [];

for i = 1:length(fn)
    fprintf('%d - %s\n',i,fn{i})
    Subject = {sprintf('SW%02d',i)};

    dt.(Subject{1}).filename = fn{i};

    data = readtable(fullfile(SW_data,fn{i}));%,'Format',columnFormat);
    
    % Convert time to datetime
    D = day(data.Time);
    M = month(data.Time);
    Y = year(data.Time);

    daySub = datetime(Y(1),M(1),D(1));

    timeDur = data.Time - daySub;

    t = data.Date + timeDur;
    
    dt.(Subject{1}).t = t;
    dt.(Subject{1}).I = data.WhiteLight;

    ind_DLMO = contains(SubInfo.Subject,fn{i}(1:end-4));
    
    dt.(Subject{1}).DLMO = DLMO(ind_DLMO,:);

    %% Plot data
    if isPlot_SW
        figure
        ax(1) = subplot(2,1,1);
        plot(dt.(Subject{1}).t,dt.(Subject{1}).I)
        hold on
        plot(dt.(Subject{1}).DLMO(1),200,'*g','MarkerSize',10)
        title(sprintf('%d - %s\n',i,fn{i}))
    %     ax(2) = subplot(2,1,2);
    %     plot(dt.(Subject{1}).t,dt.(Subject{1}).act,'k')
    %     hold on
    %     plot(dt.(Subject{1}).DLMO,50,'*','MarkerSize',10)
    %     linkaxes(ax,'x')
    end
    
end

save(fullfile(save_fn,'dataOrganized_20241003.mat'),"dt")








