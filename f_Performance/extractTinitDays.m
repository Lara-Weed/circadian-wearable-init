function T = extractTinitDays(DT)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    % Preallocate arrays
    DayfromDLMO = nan(size(DT, 1), 1);
    X_init = nan(size(DT, 1), 1);
    Xc_init = nan(size(DT, 1), 1);
    Tau = nan(size(DT, 1), 1);
    InitTime = NaT(size(DT, 1), 1);
    DLMO = NaT(size(DT, 1), 1);
    
    % Extract data from the cell array
    for j = 1:size(DT, 1)

        if ~isempty(DT{j,1})
            temp = DT{j,1};     % Get 1x4 numeric array
            DayfromDLMO(j) = temp(1);
            X_init(j) = temp(2);
            Xc_init(j) = temp(3);
            Tau(j) = temp(4);
        end

        if ~isempty(DT{j,2})
            InitTime(j) = DT{j,2};  % Get datetime
        end
        
        if ~isempty(DT{j,3})
            DLMO(j) = DT{j,3};  % Get datetime
        end
    end
    
    % Create table
    T = table(DayfromDLMO, X_init, Xc_init, Tau, InitTime, DLMO);
end