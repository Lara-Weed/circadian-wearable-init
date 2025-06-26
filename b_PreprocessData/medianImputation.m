function [medianimputed,t_new,imputed,totmissing] = medianImputation(act,t,gaps,epoch)
    
    t_new = [t(1):epoch:t(end)]';
    act_new = nan(length(t),1);

    for i = 1:length(t)
        ind_new = find(t_new == t(i));
        act_new(ind_new) = act(i);
    end

  
    medianimputed = act_new;
    tHrs = hour(t_new)+minute(t_new)/60 + second(t_new)/3600;
    for q = 1:size(gaps,1)
        
        ind = t_new> gaps(q,1) & t_new<= gaps(q,2);

        gapTime = t_new(ind);
        
        gapHrs = hour(gapTime)+minute(gapTime)/60 + second(gapTime)/3600;
        
        impval = nan(length(gapHrs),1);
        for qq = 1:length(gapHrs)
            all_vals = act_new(tHrs==gapHrs(qq));
            vals = all_vals(all_vals~=0 & ~isnan(all_vals));
            impval(qq) = median(vals);
        end
        
        medianimputed(ind) = impval;
    end

    imputed = act_new == medianimputed & ~isnan(medianimputed);
    totmissing = ~isnan(medianimputed);
end

