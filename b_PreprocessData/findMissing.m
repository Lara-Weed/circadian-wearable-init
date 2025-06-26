function [MissingT] = findMissing(t,act,epoch)

% Identify Missing Data Phenotype
    dt = diff(t);
    ind_miss = dt> epoch;

    if sum(ind_miss)>0

        start_m = t(logical([ind_miss;0]));
        stop_m = t(logical([0;ind_miss]));
        phenotype_m = repmat('Missing',length(start_m),1);
    
        % Identify Non-wear phenotype
        ind_nonwear = act==0;
    
        da = diff(ind_nonwear);
    
        ind_nwp = da>0;
        ind_nwn = da<0;
    
        start_n = t(logical([0;ind_nwp]));
        stop_n = t(logical([ind_nwn;0]));
    
        if stop_n(1) < start_n(1)
            start_n = [t(1);start_n];
        end
    
        if stop_n(end) < start_n(end)
            stop_n = [stop_n;t(end)];
        end
    
        phenotype_n = repmat('Nonwear',length(start_n),1);
    
        % Missing & Nonwear data table
        start = [start_m;start_n];
        stop = [stop_m;stop_n];
        dur = stop-start;
        phenotype = [phenotype_m;phenotype_n];
    else
        % Missing & Nonwear data table
        start = [];
        stop = [];
        dur = [];
        phenotype = [];
    end


    MissingT = table(start,stop,dur,phenotype);
    MissingT = sortrows(MissingT);

end