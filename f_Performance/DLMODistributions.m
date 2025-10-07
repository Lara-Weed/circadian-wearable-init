%% Plot the DLMO Timing Distributions 
% 10 Sept 2025
% Lara Weed

%% Load Data
basepath = '/Users/lara/Library/CloudStorage/OneDrive-Stanford/Research/Projects/Light Circadian/Data/InitsDays/stH';

rawData = load('/Users/lara/Library/CloudStorage/OneDrive-Stanford/Research/Projects/Light Circadian/Data/Preprocessed/subsetData_20241003.mat');

subs = fields(rawData.T);

DLMOS = nan(1,length(subs));
for i = 1:length(subs)
    dt = rawData.T.(subs{i}).DLMO;

    DLMOS(i) = hour(dt) + minute(dt)/60 + second(dt)/3600;
end

%%
% Assuming DLMOS is already in hours [0,24) and may contain NaNs
DLMOS = mod(DLMOS, 24);
n = numel(DLMOS);
idx1 = 1:min(56, n);
idx2 = (idx1(end)+1):n;

% Helper: map hours on a 24h circle -> signed hours in (-12, 12]
toRad   = @(h) (h/24)*2*pi;
wrapHrs = @(h) ((mod(toRad(h)+pi, 2*pi) - pi) / (2*pi)) * 24;

signed1 = wrapHrs(DLMOS(idx1(~isnan(DLMOS(idx1)))));
signed2 = wrapHrs(DLMOS(idx2(~isnan(DLMOS(idx2)))));

% 15-minute bin edges, centered at 0
edges_lin = -12:0.25:12;   % hours

ax(1) = subplot(2,1,1);
h1 = histogram(signed1, 'BinEdges', edges_lin,'DisplayStyle','bar', 'LineWidth',1.5,'FaceColor','k');

xlim([-12 12])
xticks(-12:2:12)
ylabel('No. Participants')
title('DLMO Timing Histogram')
grid on
xlabel('Hours (Relative to Midnight)')
set(gca,'FontWeight','bold','FontSize',12)

ax(2) = subplot(2,1,2);
h2 = histogram(signed2, 'BinEdges', edges_lin,'DisplayStyle','bar', 'LineWidth',1.5,'FaceColor','r');
xlim([-12 12])
xticks(-12:2:12)
linkaxes(ax,'xy')
xlabel('Hours (Relative to Midnight)')
ylabel('No. Participants')
grid on
set(gca,'FontWeight','bold','FontSize',12)

% — Optional: shape comparison independent of sample size —
% set([h1 if exist('h2','var'), h2, end], 'Normalization','probability')

%% Load Recording Length Data with 1440 initalizations

basepath = '/Users/lara/Library/CloudStorage/OneDrive-Stanford/Research/Projects/Light Circadian/Data/InitsDays/stH';
rawData = load('/Users/lara/Library/CloudStorage/OneDrive-Stanford/Research/Projects/Light Circadian/Data/Preprocessed/subsetData_20241003.mat');

d = dir(basepath);

fn = {d.name}';
fn = fn(contains(fn,'.mat'));

tt = [];
tt48 = [];
sublist2 =[];
for i = 1:length(fn)

    gg = 13;
    ss = 19;
    
    fprintf('   %s - %d/%d\n',fn{i}(ss:end-gg),i,length(fn))
    
    load(fullfile(basepath,fn{i}))

    sublist2 = [sublist2 ;{fn{i}(ss:end-gg)},{i}];
    TrueDLMO = rawData.T.(fn{i}(ss:end - gg)).DLMO;
    TrueDLMO_hr = hour(TrueDLMO) + minute(TrueDLMO)./60 + second(TrueDLMO)./3600;

    T = extractTinitDays(DT);

    T = rmmissing(T);

    uDays = unique(T.DayfromDLMO);

    T.DLMO = T.DLMO + hours(1/3); % correction to 6.667 hour offset
    %betweeen CBTmin and DLMO instead of 7 hours

    T.DLMOrel = hour(T.DLMO) + minute(T.DLMO)./60 + second(T.DLMO)./3600 - TrueDLMO_hr;
    T.DLMOrel(T.DLMOrel>12) = T.DLMOrel(T.DLMOrel>12) - 24; 
    T.DLMOrel(T.DLMOrel<-12) = T.DLMOrel(T.DLMOrel<-12) + 24; 

    tt1 = nan(length(uDays),9);
    tt1_48 = nan(length(uDays),9);
    parfor j = 1:length(uDays)

        % All 1440 inits
        TDR = T.DLMOrel(T.DayfromDLMO == uDays(j));
        [f1,xf1] = kde(TDR,"EvaluationPoints",linspace(-12,12,97));
        [f2,xf2] = kde(TDR);
        NLL = -log(f1(xf1 == 0)./(sum(f1)) + 10e-6);
        MLE = xf2(f2 == max(f2));
        RMSE = sqrt(MLE.^2);
        [mu, K, Entropy] = estimateVonMisesParamsWithEntropy(2*pi*TDR./24, 288);
        
        Mu = 24*mu./(2*pi);

        tt1(j,:) = [i,uDays(j),length(TDR),NLL,MLE,RMSE,Mu,K,Entropy];

        % 48 inits
        TDR = randsample(TDR,48,true);
        [f1,xf1] = kde(TDR,"EvaluationPoints",linspace(-12,12,97));
        [f2,xf2] = kde(TDR);
        NLL = -log(f1(xf1 == 0)./(sum(f1)) + 10e-6);
        MLE = xf2(f2 == max(f2));
        RMSE = sqrt(MLE.^2);
        [mu, K, Entropy] = estimateVonMisesParamsWithEntropy(2*pi*TDR./24, 288);
        
        Mu = 24*mu./(2*pi);

        tt1_48(j,:) = [i,uDays(j),length(TDR),NLL,MLE,RMSE,Mu,K,Entropy];


    end
    tt = [tt;tt1];
    tt48 = [tt48;tt1_48];
end

% All inits
uSubs = unique(tt(:,1));
TT.NLL = nan(length(uSubs),round(max(tt(:,2))));
TT.MLE = nan(length(uSubs),round(max(tt(:,2))));
TT.RMSE = nan(length(uSubs),round(max(tt(:,2))));
TT.Mu = nan(length(uSubs),round(max(tt(:,2))));
TT.K = nan(length(uSubs),round(max(tt(:,2))));
TT.Entropy = nan(length(uSubs),round(max(tt(:,2))));

for i = 1:length(uSubs)
    ind = tt(:,1) == uSubs(i);
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

% 48 inits
uSubs = unique(tt48(:,1));
TT48.NLL = nan(length(uSubs),round(max(tt48(:,2))));
TT48.MLE = nan(length(uSubs),round(max(tt48(:,2))));
TT48.RMSE = nan(length(uSubs),round(max(tt48(:,2))));
TT48.Mu = nan(length(uSubs),round(max(tt48(:,2))));
TT48.K = nan(length(uSubs),round(max(tt48(:,2))));
TT48.Entropy = nan(length(uSubs),round(max(tt48(:,2))));

for i = 1:length(uSubs)
    ind = tt48(:,1) == uSubs(i);
    for j = 1:round(max(tt48(:,2)))
        subtt = tt48(ind & tt48(:,2)==j,4:9);
        if ~isempty(subtt)
            TT48.NLL(i,j) = subtt(1);
            TT48.MLE(i,j) = subtt(2);
            TT48.RMSE(i,j) = subtt(3);
            TT48.Mu(i,j) = subtt(4);
            TT48.K(i,j) = subtt(5);
            TT48.Entropy(i,j) = subtt(6);
        end
    end

end

%% Compute Lin's Bias for 14 day recording with 1440 initializations
% Match DLMOs and Subs
x = DLMOS(1:end-1)';
y = x+TT.MLE(:,14); % Day 14
y2 = x+TT48.MLE(:,14);

XY = rmmissing([x,y]);
XY2 = rmmissing([x,y2]);

X = XY(:,1);
Y = XY(:,2);

X2 = XY2(:,1);
Y2 = XY2(:,2);

% RS
LB_RS_14days_1440inits = LinsBias(X(1:56),Y(1:56))

% SW
LB_SW_14days_1440inits = LinsBias(X(57:end),Y(57:end))

% RS
LB_RS_14days_48inits = LinsBias(X2(1:56),Y2(1:56))

% SW
LB_SW_14days_48inits = LinsBias(X2(57:end),Y2(57:end))

