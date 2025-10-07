function LB = LinsBias(x,y)
% Computes Lin's bias from measured (x) and predicted (y) data 
    p = corrcoef(x,y);
    cc = p(1,2);
    LB = 2*cc*std(x)*std(y)/(var(x) + var(y) + (mean(x) - mean(y))^2);
end