function runLightModel_Wearable_ODE15_AMP(sub)

	%% Wearable Light Diets
	% Lara Weed
	% 4 Apr 2025

	%% Add Paths
    addpath('/scratch/groups/jzeitzer/TimeZones/ProcessL')
	addpath('/scratch/groups/jzeitzer/CircadianModeling/New')

    %% Set starting phases
   
    % Evenly spaced phase
    phases = linspace(0,2*pi,1441)';
    phases = phases(1:end-1);


    % Set intialization points (1440)
    % Amplitude = 1
    x_inits = sin(phases);
    xc_inits = cos(phases);

    %% Load Participant Data
    sub = round(sub);

    load('/scratch/groups/jzeitzer/CircadianModeling/New/subsetData_20241003.mat')
    Subjects = fields(T);

    Subject = Subjects{sub};

    DLMO_true = T.(Subjects{sub}).DLMO;
    t_all = T.(Subjects{sub}).t;
    I_all = T.(Subjects{sub}).I;

    [B_hat_all,~,~,~] = processL_stHilaire2007(I_all,t_all);
    

    %% Set Starting Times

     totDays = days(DLMO_true - t_all(1));


     DaysFromDLMO = [1:floor(totDays),totDays];
     startTimes = DLMO_true - days(DaysFromDLMO);


    %% All inits

    all_inits =[];
    for kk = 1:length(startTimes)

        startTime = repmat(startTimes(kk),length(x_inits),1);
        DayFromDLMO = repmat(DaysFromDLMO(kk),length(x_inits),1);
        all_inits = [all_inits; table(DayFromDLMO, startTime, x_inits, xc_inits)];
    end

	
	DT =cell(size(all_inits,1),3);

	parfor pt = 1:size(all_inits,1)
		try
    		fprintf('    %d/%d\n',pt,size(all_inits,1))
     
    		% Process P
            X_init = all_inits.x_inits(pt);
            Xc_init = all_inits.xc_inits(pt);
            Tx = 24.2;

            % Select initialization scheme
            days_ind = t_all>= all_inits.startTime(pt);

            t_pt = t_all(days_ind);
        
            % Prep for ODE
            tVals = hours(t_pt - t_pt(1));
            
            bVals = B_hat_all(days_ind);
        
            tSpan   = tVals;  % numeric start/end
        
            y0 = [X_init;Xc_init];
        
            [~, ySol] = ode15s(@(xx, yy) processPODE(xx, yy, tVals, bVals , Tx), tSpan, y0);

            % Extract States
            X = ySol(:,1);
            Xc = ySol(:,2);
           
            % Compute phase
    		phaseSol = atan2(Xc,X);
        
            % Extract CBTmin estimates
            p_ref = hours(0.97); % St. Hilaire
        
            p_xcx = deg2rad(-170.7); % St. Hilaire
        
            CBTMins = t_pt(phaseSol(1:end-1)<= p_xcx & phaseSol(2:end)>= p_xcx & diff(phaseSol)> 0) + p_ref; % Model Predicts CBTmin

            % Convert to CBTmin (6.667 hours from St. Hilaire et al. 2007)
            DLMOs =  CBTMins - hours(6.667); %CBTmin - 7 hours = DLMO (approximately) 
            
            % Populate Iterations
            tmpRow = cell(1,3);
            tmpRow{1} = [all_inits.DayFromDLMO(pt), all_inits.x_inits(pt), all_inits.xc_inits(pt),24.2];
            tmpRow{2} = all_inits.startTime(pt);
            tmpRow{3} = DLMOs(end);
    
            % Store that row back into DLMOForms
            DT(pt,:) = tmpRow;
		catch
		end
	end

    fprintf('    Saving Data...\n')
	spath = '/scratch/groups/jzeitzer/CircadianModeling/New/Results/StH_AMP';

	sfn1 = sprintf('Results_ODE15_AMP_StH_%s_20250422.mat',Subject);
	save(fullfile(spath,sfn1),"DT","-v7.3")

end
