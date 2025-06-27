function [IV,IS] = calcIVIS(act,t)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% IS calculation
N = length(act);
p = hours(24)/(t(2)-t(1));
X_bar = nanmean(act);
h_ind = hour(t);%(24*(day(t)-1) + hour(t)) - min((24*(day(t)-1) + hour(t))) + 1;
hi_ind = (24*(day(t)-1) + hour(t)) - ((24*(day(t(1))-1) + hour(t(1)))) + 1;
if sum(hi_ind<0)>0
    hi_ind(hi_ind<0) = hi_ind(hi_ind<0)+ hi_ind(find(hi_ind<0,1)-1) - hi_ind(find(hi_ind<0,1)); 
end

X_i = nan(max(hi_ind),1);
for h = 1:max(hi_ind)
    X_i(h) = nanmean(act(hi_ind == h));   
end

X_h = nan(max(h_ind),1);
for hh = 1:length(unique(h_ind))
    X_h(hh) = nanmean(act(h_ind == hh-1));   
end

topIS =nansum((X_h - X_bar).^2);
bottomIS = nansum((X_i-X_bar).^2);%sum((act-X_bar).^2);

IS = (N*topIS)./(p*bottomIS);

% IV calculation
topIV = nansum(diff(X_i).^2);
bottomIV = nansum((X_i-X_bar).^2);%sum((act-X_bar).^2);

IV = N*topIV./((N-1)*bottomIV);


end

