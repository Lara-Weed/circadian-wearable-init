%% Preprocess Data
% This includes identifying missing data, performing median imputation, and
% and potentially also other things

%% Paths
save_path = '/Users/lara/Library/CloudStorage/OneDrive-Stanford/Research/Projects/Light Circadian/Data/Preprocessed';

is_plot = 0;

%% Load data
T = load('/Users/lara/Library/CloudStorage/OneDrive-Stanford/Research/Projects/Light Circadian/Data/Organized/dataOrganized_20241003.mat');

Subjects = fields(T.dt);

for i = [1:115,117:length(Subjects)]% 116 - timestamp issue
    clear t I act
    fprintf('%s\n',Subjects{i})
    t = T.dt.(Subjects{i}).t;
    I = T.dt.(Subjects{i}).I;

    epoch = median(diff(t));
    
    
    if isfield(T.dt.(Subjects{i}),'act')
        act = T.dt.(Subjects{i}).act;
        
            % Find missing data
            MissingT = findMissing(t,act,epoch);
            if ~isempty(MissingT)
                MissingT = MissingT(MissingT.dur>duration(1,0,0),:);
            
                missing = [missing;sum(MissingT.dur)];
                gaps = [MissingT.start, MissingT.stop];
            
                % Criteria
                min_length = t(end) - t(1)>= duration(3*24,0,0); % at least 3 days of data
                max_missing = (sum(MissingT.dur)./(t(end)-t(1)))<(24/168); % Max ratio of missing data 
            
                if isempty(MissingT) && min_length %&& max_missing
        
                    
        %             for jj = 1:length(MissingT.dur)
        %                 ind = t>=MissingT.start(jj) & t<=MissingT.stop(jj);
        %                 I(ind) = nan;
        %             end
        
                    % Remove missing light data
                    [IImp,timeImp,imputedLight,~]  = medianImputation(I,t,gaps,epoch);
        
                    % Impute Missing activity data
                    [actImp,~,imputedAct,totmissingAct]= medianImputation(act,t,gaps,epoch);
                    
                    % new struct
                    dt.(Subjects{i}).t = timeImp;
                    dt.(Subjects{i}).act = actImp;
                    dt.(Subjects{i}).I = IImp;
                    dt.(Subjects{i}).DLMO = T.dt.(Subjects{i}).DLMO;
                    dt.(Subjects{i}).filename = T.dt.(Subjects{i}).filename;
                    %dt.(Subjects{i}).predicted_CBTmin_act = predicted_CBTmin_act;
                    dt.(Subjects{i}).imputedLight = imputedLight;
                    dt.(Subjects{i}).imputedAct = imputedAct;
                    dt.(Subjects{i}).totmissingAct = totmissingAct;
                end
            else
                dt.(Subjects{i}).t = t;
                dt.(Subjects{i}).act = act;
                dt.(Subjects{i}).I = I;
                dt.(Subjects{i}).DLMO = T.dt.(Subjects{i}).DLMO;
                dt.(Subjects{i}).filename = T.dt.(Subjects{i}).filename;
                %dt.(Subjects{i}).predicted_CBTmin_act = predicted_CBTmin_act;
                dt.(Subjects{i}).imputedLight = I;
                dt.(Subjects{i}).imputedAct = act;
                dt.(Subjects{i}).totmissingAct = zeros(length(I),1);
            end
    
    % %             figure
    % %             axR(1) = subplot(2,1,1);
    % %             plot(timeImp,IImp)
    % %             hold on
    % %             plot(t,I,'k')
    % %             ylabel('Light (lux)')
    % %             ylim([0 3000])
    % % 
    % %             title(sprintf('%s\n',Subjects{i}))
    % %             axR(2) = subplot(2,1,2);
    % %             plot(timeImp,actImp)
    % %             hold on
    % %             plot(t,act)
    % %             ylim([0 400])
    % %             ylabel('Actigraphy (Counts)')
    % %             linkaxes(axR,'x')
        %end

    % No actigraphy data
    % The shiftwork datset (as is) does nto contain any actigraphy data.
    % Additionally, since there is not a reiable way to estimate Initial
    % CBTmin, we will skip this step and initially start from wherever the data starts 
    else 
        
            % Find missing data
            % since we have no actigaphy yet, let's use light exposure but
            % potentially with relaxed criteria
            MissingT = findMissing(t,I,epoch);
            if ~isempty(MissingT)
                MissingT = MissingT(MissingT.dur>duration(1,0,0),:);
            
                missing = [missing;sum(MissingT.dur)];
                gaps = [MissingT.start, MissingT.stop];
            
                % Criteria
                min_length = t(end) - t(1)>= duration(3*24,0,0); % at least 3 days of data
                max_missing = (sum(MissingT.dur)./(t(end)-t(1)))<(24/168); % Max ratio of missing data 
            
                if ~isempty(MissingT) && min_length %&& max_missing
        
                    % Remove missing light data
                    [IImp,timeImp,imputedLight,totmissingLight] = medianImputation(I,t,gaps,epoch);
                
                    % new struct
                    dt.(Subjects{i}).t = timeImp;
                    dt.(Subjects{i}).I = IImp;
                    dt.(Subjects{i}).DLMO = T.dt.(Subjects{i}).DLMO;
                    dt.(Subjects{i}).filename = T.dt.(Subjects{i}).filename;
                    dt.(Subjects{i}).imputedLight = imputedLight;
                    dt.(Subjects{i}).totmissingLight = totmissingLight;
                    %dt.(Subjects{i}).predicted_CBTmin_act = predicted_CBTmin_act;
                end
            else
                % new struct
                dt.(Subjects{i}).t = t;
                dt.(Subjects{i}).I = I;
                dt.(Subjects{i}).DLMO = T.dt.(Subjects{i}).DLMO;
                dt.(Subjects{i}).filename = T.dt.(Subjects{i}).filename;
                dt.(Subjects{i}).imputedLight = I;
                dt.(Subjects{i}).totmissingLight = zeros(length(I),1);
            end
    end
    if is_plot
                figure
    % %             axR(1) = subplot(2,1,1);
    % %             plot(timeImp,IImp)
    % %             hold on
                 plot(dt.(Subjects{i}).t,dt.(Subjects{i}).I,'k')
                 ylabel('Light (lux)')
    % %             ylim([0 3000])
    % % 
                 title(sprintf('%s\n',Subjects{i}))
    % %             axR(2) = subplot(2,1,2);
    % %             plot(timeImp,actImp)
    % %             hold on
    % %             plot(t,act)
    % %             ylim([0 400])
    % %             ylabel('Actigraphy (Counts)')
    % %             linkaxes(axR,'x')
    end
end

clear T

Subjects = fields(dt);

for i = 1:length(Subjects)  
    % Set variables
    t = dt.(Subjects{i}).t;
    clear act
    if isfield(dt.(Subjects{i}),'act')
        act = dt.(Subjects{i}).act;
    end
    I = dt.(Subjects{i}).I;
    %predicted_min = dt.(Subjects{i}).predicted_CBTmin_act;
    fi_i = dt.(Subjects{i}).DLMO + hours(8);

    if (t(1) <= fi_i && t(end)+days(1) >= fi_i && fi_i - t(1) > days(5))
%         fprintf('%d - %s\n',i,Subjects{i})
%         figure
%         plot(t(I>1),log(I(I>1)))
%         hold on
%         plot(fi_i,1,'k*','linewidth',2,'MarkerSize',10)
%         title(Subjects{i})
        T.(Subjects{i}) = dt.(Subjects{i});
    end

end

save(fullfile(save_path,'preprocessed_20241003.mat'),'T')












